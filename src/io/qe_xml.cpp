#include "zinc/io/qe_xml.h"
#include "zinc/core/constants.h"
#include <pugixml.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cctype>
#include <iomanip>

namespace zinc {
namespace io {

bool QEXmlDriver::identify(const std::string& first_line) const {
  std::string lower = first_line;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  if (lower.find('<') == std::string::npos) return false;
  return lower.find("qes:espresso") != std::string::npos ||
         lower.find("xmlns:qes") != std::string::npos ||
         lower.find("qexsd") != std::string::npos ||
         lower.find("units=\"hartree atomic units\"") != std::string::npos ||
         lower.find("<atomic_structure") != std::string::npos;
}

// -- helpers --

static std::string local_name(const std::string& name) {
  auto pos = name.find(':');
  return (pos != std::string::npos) ? name.substr(pos + 1) : name;
}

static pugi::xml_node child_by_local_name(pugi::xml_node node, const char* name) {
  for (auto child : node.children()) {
    if (local_name(child.name()) == name) return child;
  }
  return {};
}

static void collect_elements(pugi::xml_node node, std::vector<pugi::xml_node>& nodes) {
  for (auto child : node.children()) {
    if (child.type() == pugi::node_element) {
      nodes.push_back(child);
      collect_elements(child, nodes);
    }
  }
}

static std::string attr_val(pugi::xml_node node, const char* attr) {
  auto a = node.attribute(attr);
  return a ? a.value() : "";
}

static double attr_double(pugi::xml_node node, const char* attr, double def = 0.0) {
  auto s = attr_val(node, attr);
  if (s.empty()) return def;
  try { return std::stod(s); } catch (...) { return def; }
}

static int attr_int(pugi::xml_node node, const char* attr, int def = 0) {
  auto s = attr_val(node, attr);
  if (s.empty()) return def;
  try { return std::stoi(s); } catch (...) { return def; }
}

static core::Vector3d parse_vec3_str(const std::string& text) {
  std::istringstream ss(text);
  double x = 0, y = 0, z = 0;
  ss >> x >> y >> z;
  return core::Vector3d(x, y, z);
}

static bool parse_vec3_arr(const std::string& text, double out[3]) {
  std::istringstream ss(text);
  return static_cast<bool>(ss >> out[0] >> out[1] >> out[2]);
}

enum class EnergyUnit { Hartree, Rydberg };
enum class ForceUnit { HartreePerBohr, RydbergPerBohr };
enum class PosType { Cartesian, Fractional };

struct PendingAtom {
  std::string name;
  double coords[3];
  PosType unit;
};

// Convert cell vector text to angstroms
static core::Vector3d cell_to_angstrom(const std::string& text, const std::string& unit, double alat) {
  auto v = parse_vec3_str(text);
  double s = 1.0;
  std::string u = unit;
  std::transform(u.begin(), u.end(), u.begin(), ::tolower);
  if (u == "bohr" || u == "a.u.") s = core::BOHR_TO_ANGSTROM;
  else if (u == "alat") s = alat * core::BOHR_TO_ANGSTROM;
  return v * s;
}

// Convert atom coordinates to angstroms
static core::Vector3d pos_to_angstrom(const double coords[3], const std::string& unit,
                                       double alat, PosType ptype,
                                       const core::Lattice& lattice) {
  core::Vector3d v(coords[0], coords[1], coords[2]);
  if (ptype == PosType::Fractional) {
    return lattice.frac_to_cart(v);
  }
  std::string u = unit;
  std::transform(u.begin(), u.end(), u.begin(), ::tolower);
  double s = 1.0;
  if (u == "bohr" || u == "a.u.") s = core::BOHR_TO_ANGSTROM;
  else if (u == "alat") s = alat * core::BOHR_TO_ANGSTROM;
  return v * s;
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

core::Structure QEXmlDriver::read(const std::filesystem::path& path) const {
  pugi::xml_document doc;
  if (!doc.load_file(path.c_str())) {
    throw std::runtime_error("Failed to parse XML: " + path.string());
  }

  auto root = doc.child("qes:espresso");
  if (!root) {
    for (auto c : doc.children()) {
      if (local_name(c.name()) == "espresso") { root = c; break; }
    }
  }
  if (!root) root = doc.first_child();

  std::string calc_mode = "scf";
  double energy_conv_thr_ry = 0.0;
  double force_conv_thr_ry_bohr = 0.0;
  double press_conv_thr = 0.0;
  double last_etot_ry = 0.0;
  bool has_etot = false;

  std::vector<core::Structure> frames;
  double global_alat = 1.0;
  std::vector<pugi::xml_node> elements;
  collect_elements(root, elements);

  // parse thresholds and global info
  for (auto node : elements) {
    std::string ln = local_name(node.name());

    if (ln == "calculation") {
      std::string val = node.child_value();
      if (!val.empty()) {
        std::string lowered = val;
        std::transform(lowered.begin(), lowered.end(), lowered.begin(), ::tolower);
        if (!lowered.empty()) calc_mode = lowered;
      }
    }
    if (ln == "etot_conv_thr") {
      try { energy_conv_thr_ry = std::stod(node.child_value()) * 2.0; } catch (...) {}
    }
    if (ln == "forc_conv_thr") {
      try { force_conv_thr_ry_bohr = std::stod(node.child_value()); } catch (...) {}
    }
    if (ln == "press_conv_thr") {
      try { press_conv_thr = std::stod(node.child_value()); } catch (...) {}
    }
    if (ln == "etot") {
      std::string vu = attr_val(node, "units");
      std::string vu_l = vu;
      std::transform(vu_l.begin(), vu_l.end(), vu_l.begin(), ::tolower);
      EnergyUnit eu = (vu_l.find("ryd") != std::string::npos || vu_l.find("ry") != std::string::npos)
                      ? EnergyUnit::Rydberg : EnergyUnit::Hartree;
      try {
        double val = std::stod(node.child_value());
        last_etot_ry = (eu == EnergyUnit::Hartree) ? val * 2.0 : val;
        has_etot = true;
      } catch (...) {}
    }
  }

  // accumulate forces and constraints; associate with next atomic_structure
  std::vector<std::array<double, 3>> pending_forces;
  std::vector<std::array<int, 3>> pending_constraints;

  for (auto node : elements) {
    std::string ln = local_name(node.name());

    if (ln == "atomic_structure") {
      double alat = attr_double(node, "alat", 1.0);
      if (global_alat == 1.0) global_alat = alat;
      int nat = attr_int(node, "nat", 0);

      // cell
      auto cell = child_by_local_name(node, "cell");
      std::string cell_unit_str = cell ? attr_val(cell, "units") : "bohr";
      if (cell_unit_str.empty()) cell_unit_str = "bohr";
      core::Vector3d v1(1,0,0), v2(0,1,0), v3(0,0,1);
      if (cell) {
        auto a1 = child_by_local_name(cell, "a1"); if (a1) v1 = cell_to_angstrom(a1.child_value(), cell_unit_str, alat);
        auto a2 = child_by_local_name(cell, "a2"); if (a2) v2 = cell_to_angstrom(a2.child_value(), cell_unit_str, alat);
        auto a3 = child_by_local_name(cell, "a3"); if (a3) v3 = cell_to_angstrom(a3.child_value(), cell_unit_str, alat);
      }

      // positions
      auto apos = child_by_local_name(node, "atomic_positions");
      std::string pos_unit_str = apos ? attr_val(apos, "units") : "bohr";
      if (pos_unit_str.empty()) pos_unit_str = "bohr";
      PosType ptype = PosType::Cartesian;
      {
        std::string pl = pos_unit_str;
        std::transform(pl.begin(), pl.end(), pl.begin(), ::tolower);
        if (pl == "crystal" || pl == "crystal_sg") ptype = PosType::Fractional;
      }

      core::Lattice lattice(v1, v2, v3);
      core::Structure structure(lattice);

      // atoms
      int atom_idx = 0;
      if (!apos) {
        throw std::runtime_error("atomic_structure missing atomic_positions in XML");
      }
      for (auto atom : apos.children()) {
        if (local_name(atom.name()) != "atom") continue;
        if (nat > 0 && atom_idx >= nat) break;
        std::string at_name = attr_val(atom, "name");

        // position from attribute or text
        std::string pos_str = attr_val(atom, "position");
        double coords[3] = {0,0,0};
        bool has_pos = false;
        if (!pos_str.empty()) {
          has_pos = parse_vec3_arr(pos_str, coords);
        }
        if (!has_pos) {
          pos_str = atom.child_value();
          has_pos = parse_vec3_arr(pos_str, coords);
        }
        if (!has_pos) {
          // try the atom itself (empty element with position attr)
          continue;
        }

        core::Vector3d cart = pos_to_angstrom(coords, pos_unit_str, alat, ptype, lattice);

        std::optional<std::array<int, 3>> constr;
        if (atom_idx < static_cast<int>(pending_constraints.size())) {
          constr = pending_constraints[atom_idx];
        }
        structure.add_atom(at_name, {cart.x(), cart.y(), cart.z()}, constr);
        atom_idx++;
      }

      // attach pending forces
      if (!pending_forces.empty()) {
        for (size_t i = 0; i < pending_forces.size() && i < structure.atoms.size(); ++i) {
          auto& f = pending_forces[i];
          structure.atoms[i].force = core::Vector3d(
            f[0] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM,
            f[1] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM,
            f[2] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM);
        }
      }

      frames.push_back(std::move(structure));
      pending_forces.clear();
      pending_constraints.clear();
    }

    if (ln == "forces") {
      pending_forces.clear();
      std::string fu_s = attr_val(node, "units");
      std::string fu_l = fu_s;
      std::transform(fu_l.begin(), fu_l.end(), fu_l.begin(), ::tolower);
      ForceUnit fu = (fu_l.find("hartree") != std::string::npos)
                     ? ForceUnit::HartreePerBohr : ForceUnit::RydbergPerBohr;

      int nat_forces = 0;
      std::string dims = attr_val(node, "dims");
      if (!dims.empty()) {
        std::istringstream dss(dims);
        int v;
        while (dss >> v) nat_forces = v;
      }

      std::string ftext = node.child_value();
      std::istringstream fss(ftext);
      std::vector<double> fbuf;
      double fv;
      while (fss >> fv) fbuf.push_back(fv);

      size_t nf = (nat_forces > 0) ? static_cast<size_t>(nat_forces) : fbuf.size() / 3;
      for (size_t i = 0; i + 2 < fbuf.size() && i / 3 < nf; i += 3) {
        double fx = fbuf[i], fy = fbuf[i+1], fz = fbuf[i+2];
        if (fu == ForceUnit::HartreePerBohr) { fx *= 2.0; fy *= 2.0; fz *= 2.0; }
        pending_forces.push_back({fx, fy, fz});
      }

      // attach to previous frame if it exists
      if (!frames.empty() && !pending_forces.empty()) {
        auto& last = frames.back();
        for (size_t i = 0; i < pending_forces.size() && i < last.atoms.size(); ++i) {
          auto& f = pending_forces[i];
          last.atoms[i].force = core::Vector3d(
            f[0] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM,
            f[1] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM,
            f[2] * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM);
        }
        pending_forces.clear();
      }
    }

    if (ln == "free_positions") {
      pending_constraints.clear();
      std::string fptext = node.child_value();
      std::istringstream fpss(fptext);
      std::vector<int> fbuf;
      int v;
      while (fpss >> v) fbuf.push_back(v);

      for (size_t i = 0; i + 2 < fbuf.size(); i += 3) {
        pending_constraints.push_back({fbuf[i], fbuf[i+1], fbuf[i+2]});
      }
    }
  }

  if (frames.empty()) {
    throw std::runtime_error("No atomic_structure found in XML");
  }

  auto frames_xyz = frames;
  auto structure = std::move(frames.back());
  frames.pop_back();

  bool is_relax = (calc_mode == "relax" || calc_mode == "vc-relax");
  structure.metadata["calculation"] = calc_mode;
  structure.metadata["frames"] = std::to_string(frames_xyz.size());

  // build extxyz trajectory
  if (is_relax) {
    size_t start = 0;
    if (frames_xyz.size() > 1) {
      bool first_no_force = true;
      for (const auto& a : frames_xyz[0].atoms) {
        if (a.force.has_value()) { first_no_force = false; break; }
      }
      if (first_no_force) start = 1;
    }
    if (start < frames_xyz.size()) {
      std::vector<core::Structure> traj(frames_xyz.begin() + start, frames_xyz.end());
      structure.metadata["trajectory_extxyz"] = frames_to_extxyz(traj);
    }
  } else if (frames_xyz.size() > 1) {
    structure.metadata["trajectory_extxyz"] = frames_to_extxyz(frames_xyz);
  }

  // energy metadata
  if (has_etot) {
    char buf[64];
    format_double(buf, sizeof(buf), "%.12f", last_etot_ry);
    structure.metadata["qe_total_energy_ry"] = buf;
    format_double(buf, sizeof(buf), "%.12f", last_etot_ry * core::RY_TO_EV);
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
    format_double(buf, sizeof(buf), "%.6e", force_conv_thr_ry_bohr * core::FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM);
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

  return structure;
}

void QEXmlDriver::write(const core::Structure& structure, const std::filesystem::path& path) const {
  std::ofstream out(path);
  if (!out.is_open()) throw std::runtime_error("Cannot open output file: " + path.string());
  out << std::fixed << std::setprecision(10);
  out << "<?xml version=\"1.0\"?>\n<qes:espresso>\n";
  out << "  <atomic_structure nat=\"" << structure.natoms() << "\">\n";
  out << "    <cell units=\"angstrom\">\n";
  const auto& mat = structure.lattice.matrix();
  for (int i = 0; i < 3; ++i)
    out << "      <a" << i+1 << ">"
        << std::setw(18) << mat(i,0)
        << std::setw(18) << mat(i,1)
        << std::setw(18) << mat(i,2)
        << "</a" << i+1 << ">\n";
  out << "    </cell>\n    <atomic_positions units=\"angstrom\">\n";
  for (const auto& atom : structure.atoms)
    out << "      <atom name=\"" << atom.element << "\">"
        << std::setw(18) << atom.position.x()
        << std::setw(18) << atom.position.y()
        << std::setw(18) << atom.position.z()
        << "</atom>\n";
  out << "    </atomic_positions>\n  </atomic_structure>\n</qes:espresso>\n";
}

} // namespace io
} // namespace zinc
