#include "zinc/io/qe_log.h"
#include "zinc/core/constants.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cctype>
#include <optional>

namespace zinc {
namespace io {

bool QELogDriver::identify(const std::string& first_line) const {
  std::string lower = first_line;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  if (lower.find('<') != std::string::npos) return false;
  return lower.find("program pwscf") != std::string::npos ||
         lower.find("quantum espresso") != std::string::npos ||
         lower.find("pwscf") != std::string::npos ||
         lower.find("espresso") != std::string::npos;
}

// -- helpers --

static bool try_parse_double(const std::string& s, double& out) {
  try { out = std::stod(s); return true; } catch (...) { return false; }
}

static bool try_parse_int(const std::string& s, int& out) {
  try { out = std::stoi(s); return true; } catch (...) { return false; }
}

static core::Vector3d parse_vec3_from_line(const std::string& line) {
  std::istringstream ss(line);
  double x = 0, y = 0, z = 0;
  ss >> x >> y >> z;
  return core::Vector3d(x, y, z);
}

static std::pair<core::LengthUnit, double> parse_cell_unit(const std::string& header, double alat_bohr) {
  std::string lower = header;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

  auto lp = header.find('(');
  auto rp = header.rfind(')');
  if (lp != std::string::npos && rp != std::string::npos && rp > lp) {
    std::string u = header.substr(lp + 1, rp - lp - 1);
    std::string ul = u;
    std::transform(ul.begin(), ul.end(), ul.begin(), ::tolower);
    if (ul.find("angstrom") != std::string::npos) return {core::LengthUnit::Angstrom, 1.0};
    if (ul.find("bohr") != std::string::npos || ul.find("a.u") != std::string::npos) return {core::LengthUnit::Bohr, 1.0};
    if (ul.find("alat=") != std::string::npos) {
      size_t eq = ul.find('=');
      std::string val = ul.substr(eq + 1);
      size_t vs = val.find_first_not_of(" \t");
      if (vs != std::string::npos) val = val.substr(vs);
      double alat_val;
      if (std::istringstream(val) >> alat_val) return {core::LengthUnit::Alat, alat_val};
      return {core::LengthUnit::Alat, alat_bohr};
    }
    if (ul.find("alat") != std::string::npos) return {core::LengthUnit::Alat, alat_bohr};
  }

  // check second space-separated token
  std::istringstream ss(header);
  std::string tok;
  ss >> tok >> tok;
  {
    std::string tl = tok;
    std::transform(tl.begin(), tl.end(), tl.begin(), ::tolower);
    if (tl.find("angstrom") != std::string::npos) return {core::LengthUnit::Angstrom, 1.0};
    if (tl.find("bohr") != std::string::npos || tl.find("a.u") != std::string::npos) return {core::LengthUnit::Bohr, 1.0};
    if (tl.find("alat") != std::string::npos) return {core::LengthUnit::Alat, alat_bohr};
  }

  if (lower.find("alat") != std::string::npos) return {core::LengthUnit::Alat, alat_bohr};
  return {core::LengthUnit::Bohr, 1.0};
}

static core::LengthUnit parse_position_unit_cart(const std::string& header) {
  std::string lower = header;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

  auto lp = header.find('(');
  auto rp = header.rfind(')');
  if (lp != std::string::npos && rp != std::string::npos && rp > lp) {
    std::string u = header.substr(lp + 1, rp - lp - 1);
    std::string ul = u;
    std::transform(ul.begin(), ul.end(), ul.begin(), ::tolower);
    if (ul.find("angstrom") != std::string::npos) return core::LengthUnit::Angstrom;
    if (ul.find("bohr") != std::string::npos || ul.find("a.u") != std::string::npos) return core::LengthUnit::Bohr;
    if (ul.find("alat") != std::string::npos) return core::LengthUnit::Alat;
  }

  std::istringstream ss(header);
  std::string tok;
  ss >> tok >> tok;
  {
    std::string tl = tok;
    std::transform(tl.begin(), tl.end(), tl.begin(), ::tolower);
    if (tl.find("angstrom") != std::string::npos) return core::LengthUnit::Angstrom;
    if (tl.find("bohr") != std::string::npos || tl.find("a.u") != std::string::npos) return core::LengthUnit::Bohr;
    if (tl.find("alat") != std::string::npos) return core::LengthUnit::Alat;
  }

  if (lower.find("alat") != std::string::npos) return core::LengthUnit::Alat;
  return core::LengthUnit::Bohr;
}

static bool is_position_fractional(const std::string& header) {
  std::string lower = header;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  return lower.find("crystal") != std::string::npos;
}

static std::optional<double> extract_last_number(const std::string& line) {
  std::istringstream ss(line);
  std::string tok;
  double last = 0;
  bool has = false;
  while (ss >> tok) {
    double v;
    if (try_parse_double(tok, v)) { last = v; has = true; }
  }
  if (has) return last;
  return std::nullopt;
}

struct AtomLineParsed {
  std::string element;
  double coords[3];
  std::array<int, 3> flags;
  bool has_flags = false;
};

static AtomLineParsed parse_atom_line(const std::string& line) {
  AtomLineParsed al;
  al.coords[0] = al.coords[1] = al.coords[2] = 0;
  al.flags = {1, 1, 1};

  std::istringstream ss(line);
  ss >> al.element;

  std::vector<std::string> tokens;
  std::string t;
  while (ss >> t) tokens.push_back(t);

  std::vector<double> nums;
  for (const auto& tok : tokens) {
    double v;
    if (try_parse_double(tok, v)) nums.push_back(v);
  }
  if (nums.size() >= 3) {
    al.coords[0] = nums[0];
    al.coords[1] = nums[1];
    al.coords[2] = nums[2];
  }

  // check for trailing 0/1 flags
  std::vector<int> iflags;
  for (size_t k = tokens.size(); k > 0 && iflags.size() < 3; --k) {
    int v;
    if (try_parse_int(tokens[k-1], v) && (v == 0 || v == 1)) {
      iflags.push_back(v);
    } else {
      break;
    }
  }
  if (iflags.size() == 3) {
    al.flags = {iflags[2], iflags[1], iflags[0]};
    al.has_flags = true;
  }

  return al;
}

static core::Vector3d parse_axes_line(const std::string& line) {
  auto eq = line.find('=');
  auto lp = (eq != std::string::npos) ? line.find('(', eq) : line.find('(');
  auto rp = (lp != std::string::npos) ? line.find(')', lp + 1) : std::string::npos;
  if (lp != std::string::npos && rp != std::string::npos && rp > lp + 1) {
    return parse_vec3_from_line(line.substr(lp + 1, rp - lp - 1));
  }
  // fallback: last 3 numbers
  std::vector<double> nums;
  std::istringstream ss(line);
  std::string tok;
  while (ss >> tok) {
    double v;
    if (try_parse_double(tok, v)) nums.push_back(v);
  }
  if (nums.size() >= 3) {
    return core::Vector3d(nums[nums.size()-3], nums[nums.size()-2], nums[nums.size()-1]);
  }
  return core::Vector3d(0, 0, 0);
}

static std::array<double, 3> parse_force_line(const std::string& line) {
  std::vector<double> nums;
  std::istringstream ss(line);
  std::string tok;
  while (ss >> tok) {
    double v;
    if (try_parse_double(tok, v)) nums.push_back(v);
  }
  if (nums.size() >= 3) {
    return {nums[nums.size()-3], nums[nums.size()-2], nums[nums.size()-1]};
  }
  return {0, 0, 0};
}

static std::string frames_to_extxyz(const std::vector<core::Structure>& frames) {
  std::ostringstream out;
  out << std::setprecision(12);
  for (size_t idx = 0; idx < frames.size(); ++idx) {
    const auto& frame = frames[idx];
    size_t nat = frame.atoms.size();
    const auto& mat = frame.lattice.matrix();
    bool has_force = std::any_of(frame.atoms.begin(), frame.atoms.end(),
                                  [](const auto& a) { return a.force.has_value(); });
    out << nat << "\nLattice=\"";
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        out << std::fixed << mat(i, j) << (i == 2 && j == 2 ? "" : " ");
    out << "\" Properties=species:S:1:pos:R:3" << (has_force ? ":forces:R:3" : "")
        << " frame=" << idx << "\n";

    out << std::setprecision(9);
    for (const auto& atom : frame.atoms) {
      out << atom.element << " " << atom.position.x() << " " << atom.position.y() << " " << atom.position.z();
      if (has_force && atom.force) {
        out << " " << atom.force->x() << " " << atom.force->y() << " " << atom.force->z();
      }
      out << "\n";
    }
    out << std::setprecision(12);
  }
  return out.str();
}

static void format_double(char* buf, size_t bufsz, const char* fmt, double val) {
  snprintf(buf, bufsz, fmt, val);
}

core::Structure QELogDriver::read(const std::filesystem::path& path) const {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + path.string());
  }

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(file, line)) {
    lines.push_back(line);
  }

  double alat_bohr = 1.0;
  int nat = 0;
  double total_energy_ry = 0.0;
  bool has_total_energy = false;
  double energy_conv_thr_ry = 0.0;
  double force_conv_thr_ry_bohr = 0.0;
  double press_conv_thr = 0.0;
  bool scf_converged = false;
  bool relax_converged = false;
  bool job_done = false;

  // first pass: global info
  for (const auto& l : lines) {
    std::string lower = l;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower.find("number of atoms/cell") != std::string::npos) {
      std::istringstream ss(l);
      std::string tok;
      while (ss >> tok) {
        int n;
        if (try_parse_int(tok, n)) nat = n;
      }
    }
    if (lower.find("lattice parameter (alat)") != std::string::npos) {
      std::istringstream ss(l);
      std::string tok;
      while (ss >> tok) {
        if (tok == "=") {
          ss >> tok;
          double v;
          if (try_parse_double(tok, v)) alat_bohr = v;
          break;
        }
      }
    }
    if (lower.find("convergence has been achieved") != std::string::npos) scf_converged = true;
    if (lower.find("bfgs converged") != std::string::npos ||
        lower.find("reached required accuracy") != std::string::npos) relax_converged = true;
    if (lower.find("job done") != std::string::npos) job_done = true;

    if (lower.find("total energy") != std::string::npos &&
        lower.find("thresh") == std::string::npos &&
        l.find('=') != std::string::npos) {
      auto parsed = extract_last_number(l);
      if (parsed.has_value()) {
        total_energy_ry = *parsed;
        has_total_energy = true;
      }
    }
    if (lower.find("energy convergence thresh") != std::string::npos ||
        lower.find("etot_conv_thr") != std::string::npos) {
      auto parsed = extract_last_number(l);
      if (parsed.has_value()) energy_conv_thr_ry = *parsed;
    }
    if (lower.find("force convergence thresh") != std::string::npos ||
        lower.find("forc_conv_thr") != std::string::npos) {
      auto parsed = extract_last_number(l);
      if (parsed.has_value()) force_conv_thr_ry_bohr = *parsed;
    }
    if (lower.find("press_conv_thr") != std::string::npos) {
      auto parsed = extract_last_number(l);
      if (parsed.has_value()) press_conv_thr = *parsed;
    }
  }

  if (nat == 0) {
    throw std::runtime_error("Missing number of atoms in log");
  }

  // second pass: frames
  std::vector<core::Structure> frames;
  core::Lattice last_lattice;
  bool has_lattice = false;

  size_t i = 0;
  while (i < lines.size()) {
    std::string li = lines[i];
    size_t ns = li.find_first_not_of(" \t\r");
    std::string lt = (ns != std::string::npos) ? li.substr(ns) : li;
    std::string lower = lt;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    // crystal axes banner
    if (lower.find("crystal axes") != std::string::npos && lower.find("cart") != std::string::npos) {
      if (i + 3 >= lines.size()) break;
      core::LengthUnit unit = core::LengthUnit::Alat;
      double s = core::to_angstrom(unit, 1.0, alat_bohr);
      last_lattice = core::Lattice(
        parse_axes_line(lines[i+1]) * s,
        parse_axes_line(lines[i+2]) * s,
        parse_axes_line(lines[i+3]) * s);
      has_lattice = true;
      i += 4;
      continue;
    }

    // CELL_PARAMETERS
    if (lower.find("cell_parameters") != std::string::npos) {
      if (i + 3 >= lines.size()) break;
      auto [unit, alat_scale] = parse_cell_unit(lt, alat_bohr);
      double s = core::to_angstrom(unit, 1.0, alat_scale);
      last_lattice = core::Lattice(
        parse_vec3_from_line(lines[i+1]) * s,
        parse_vec3_from_line(lines[i+2]) * s,
        parse_vec3_from_line(lines[i+3]) * s);
      has_lattice = true;
      i += 4;
      continue;
    }

    // ATOMIC_POSITIONS
    if (lower.find("atomic_positions") != std::string::npos) {
      if (!has_lattice) {
        throw std::runtime_error("CELL_PARAMETERS not found before ATOMIC_POSITIONS");
      }
      bool fractional = is_position_fractional(lt);
      core::LengthUnit pos_unit = fractional ? core::LengthUnit::Angstrom
                                             : parse_position_unit_cart(lt);

      if (i + 1 + nat > lines.size()) break;

      core::Structure structure(last_lattice);
      for (int k = 0; k < nat && (i + 1 + k) < static_cast<int>(lines.size()); ++k) {
        auto al = parse_atom_line(lines[i + 1 + k]);
        core::Vector3d pos(al.coords[0], al.coords[1], al.coords[2]);

        if (fractional) {
          pos = last_lattice.frac_to_cart(pos);
        } else {
          double s = core::to_angstrom(pos_unit, 1.0, alat_bohr);
          pos *= s;
        }

        std::optional<std::array<int, 3>> constr;
        if (al.has_flags) constr = al.flags;
        structure.add_atom(al.element, {pos.x(), pos.y(), pos.z()}, constr);
      }
      frames.push_back(std::move(structure));
      i += nat + 1;
      continue;
    }

    // Forces block
    if (lower.find("forces acting on atoms") != std::string::npos && lower.find("ry/au") != std::string::npos) {
      std::vector<std::array<double, 3>> forces;
      size_t j = i + 1;
      while (j < lines.size() && forces.size() < static_cast<size_t>(nat)) {
        std::string lj = lines[j];
        size_t ns2 = lj.find_first_not_of(" \t\r");
        std::string ljt = (ns2 != std::string::npos) ? lj.substr(ns2) : lj;
        std::string ljl = ljt;
        std::transform(ljl.begin(), ljl.end(), ljl.begin(), ::tolower);

        if (ljl.find("the ionic contribution") != std::string::npos) break;
        if (ljl.find("atom") != std::string::npos && ljl.find("force") != std::string::npos) {
          forces.push_back(parse_force_line(ljt));
          j++;
          continue;
        }
        if (ljt.empty()) { j++; continue; }
        break;
      }

      if (!forces.empty() && !frames.empty()) {
        auto& last = frames.back();
        for (size_t k = 0; k < forces.size() && k < last.atoms.size(); ++k) {
          auto& f = forces[k];
          last.atoms[k].force = core::Vector3d(
            f[0] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM,
            f[1] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM,
            f[2] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM);
        }
      }
      i = j;
      continue;
    }

    i++;
  }

  if (frames.empty()) {
    throw std::runtime_error("No ATOMIC_POSITIONS found in log");
  }

  auto frames_xyz = frames;
  auto structure = std::move(frames.back());
  frames.pop_back();

  structure.metadata["frames"] = std::to_string(frames_xyz.size());
  if (frames_xyz.size() > 1) {
    structure.metadata["trajectory_extxyz"] = frames_to_extxyz(frames_xyz);
  }

  if (has_total_energy) {
    char buf[64];
    format_double(buf, sizeof(buf), "%.12f", total_energy_ry);
    structure.metadata["qe_total_energy_ry"] = buf;
    format_double(buf, sizeof(buf), "%.12f", total_energy_ry * core::RY_TO_EV);
    structure.metadata["qe_total_energy_ev"] = buf;
  }
  if (energy_conv_thr_ry > 0) {
    char buf[64];
    format_double(buf, sizeof(buf), "%.6e", energy_conv_thr_ry);
    structure.metadata["qe_energy_threshold_ry"] = buf;
    format_double(buf, sizeof(buf), "%.6e", energy_conv_thr_ry * core::RY_TO_EV);
    structure.metadata["qe_energy_threshold_ev"] = buf;
  }
  if (force_conv_thr_ry_bohr > 0) {
    char buf[64];
    format_double(buf, sizeof(buf), "%.6e", force_conv_thr_ry_bohr);
    structure.metadata["qe_force_threshold_ry_bohr"] = buf;
    format_double(buf, sizeof(buf), "%.6e",
      force_conv_thr_ry_bohr * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM);
    structure.metadata["qe_force_threshold_ev_ang"] = buf;
  }
  if (press_conv_thr > 0) {
    char buf[64];
    format_double(buf, sizeof(buf), "%.6e", press_conv_thr);
    structure.metadata["qe_press_threshold"] = buf;
  }

  // max force
  {
    double max_evang = 0.0;
    for (const auto& atom : structure.atoms) {
      if (atom.force) {
        double mag = std::sqrt(
          atom.force->x() * atom.force->x() +
          atom.force->y() * atom.force->y() +
          atom.force->z() * atom.force->z());
        if (mag > max_evang) max_evang = mag;
      }
    }
    if (max_evang > 0.0) {
      char buf[64];
      format_double(buf, sizeof(buf), "%.6e", max_evang);
      structure.metadata["qe_max_force_ev_ang"] = buf;
      format_double(buf, sizeof(buf), "%.6e", max_evang / core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM);
      structure.metadata["qe_max_force_ry_bohr"] = buf;
    }
  }

  std::string status = job_done ? "completed" :
                       relax_converged ? "relax converged" :
                       scf_converged ? "scf converged" : "not converged";
  structure.metadata["qe_status"] = status;
  structure.metadata["qe_job_done"] = job_done ? "true" : "false";
  structure.metadata["qe_scf_converged"] = scf_converged ? "true" : "false";
  structure.metadata["qe_relax_converged"] = relax_converged ? "true" : "false";
  structure.metadata["calculation"] = "scf";

  return structure;
}

void QELogDriver::write(const core::Structure& structure, const std::filesystem::path& path) const {
  std::ofstream out(path);
  if (!out.is_open()) throw std::runtime_error("Cannot open output file: " + path.string());
  out << std::fixed << std::setprecision(10);
  std::string calc_name = "pwscf";
  auto it = structure.metadata.find("calculation");
  if (it != structure.metadata.end()) calc_name = it->second;

  out << "     PROGRAM: " << calc_name << "\n";
  out << "     number of atoms/cell          " << structure.natoms() << "\n";
  const auto& mat = structure.lattice.matrix();
  out << "     crystal axes (cartesian coord. - cartesian axes)\n";
  for (int i = 0; i < 3; ++i) {
    out << "a(" << i+1 << ") = ( "
        << std::setw(18) << mat(i,0)/core::BOHR_TO_ANGSTROM
        << std::setw(18) << mat(i,1)/core::BOHR_TO_ANGSTROM
        << std::setw(18) << mat(i,2)/core::BOHR_TO_ANGSTROM << " )\n";
  }
  out << "     ATOMIC_POSITIONS (angstrom)\n";
  for (const auto& atom : structure.atoms) {
    out << " " << std::left << std::setw(4) << atom.element << std::right
        << std::setw(18) << atom.position.x()
        << std::setw(18) << atom.position.y()
        << std::setw(18) << atom.position.z() << "\n";
  }
}

} // namespace io
} // namespace zinc
