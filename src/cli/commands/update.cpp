#include <iostream>
#include <chrono>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include "zinc/io/driver.h"
#include "zinc/core/structure.h"
#include "zinc/core/constants.h"

namespace fs = std::filesystem;

static double parse_fortran(const std::string& s) {
  std::string f = s; std::replace(f.begin(), f.end(), 'D', 'E'); std::replace(f.begin(), f.end(), 'd', 'e');
  try { return std::stod(f); } catch (...) { return 0.0; }
}

static std::string json_escape(std::string s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (char c : s) {
    switch (c) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: out.push_back(c); break;
    }
  }
  return out;
}

// -- atom mass database (essential elements) --
static std::unordered_map<std::string, double> atom_masses = {
  {"H", 1.008}, {"He", 4.0026}, {"Li", 6.94}, {"Be", 9.0122}, {"B", 10.81}, {"C", 12.011},
  {"N", 14.007}, {"O", 15.999}, {"F", 18.998}, {"Ne", 20.180}, {"Na", 22.990}, {"Mg", 24.305},
  {"Al", 26.982}, {"Si", 28.085}, {"P", 30.974}, {"S", 32.06}, {"Cl", 35.45}, {"Ar", 39.948},
  {"K", 39.098}, {"Ca", 40.078}, {"Sc", 44.956}, {"Ti", 47.867}, {"V", 50.942}, {"Cr", 51.996},
  {"Mn", 54.938}, {"Fe", 55.845}, {"Co", 58.933}, {"Ni", 58.693}, {"Cu", 63.546}, {"Zn", 65.38},
  {"Ga", 69.723}, {"Ge", 72.630}, {"As", 74.922}, {"Se", 78.971}, {"Br", 79.904}, {"Kr", 83.798},
  {"Rb", 85.468}, {"Sr", 87.62}, {"Y", 88.906}, {"Zr", 91.224}, {"Nb", 92.906}, {"Mo", 95.95},
  {"Tc", 98.0}, {"Ru", 101.07}, {"Rh", 102.91}, {"Pd", 106.42}, {"Ag", 107.87}, {"Cd", 112.41},
  {"In", 114.82}, {"Sn", 118.71}, {"Sb", 121.76}, {"Te", 127.60}, {"I", 126.90}, {"Xe", 131.29},
  {"Cs", 132.91}, {"Ba", 137.33}, {"La", 138.91}, {"Ce", 140.12}, {"Pr", 140.91}, {"Nd", 144.24},
  {"Pm", 145.0}, {"Sm", 150.36}, {"Eu", 151.96}, {"Gd", 157.25}, {"Tb", 158.93}, {"Dy", 162.50},
  {"Ho", 164.93}, {"Er", 167.26}, {"Tm", 168.93}, {"Yb", 173.05}, {"Lu", 174.97}, {"Hf", 178.49},
  {"Ta", 180.95}, {"W", 183.84}, {"Re", 186.21}, {"Os", 190.23}, {"Ir", 192.22}, {"Pt", 195.08},
  {"Au", 196.97}, {"Hg", 200.59}, {"Tl", 204.38}, {"Pb", 207.2}, {"Bi", 208.98}, {"Po", 209.0},
  {"At", 210.0}, {"Rn", 222.0}, {"Fr", 223.0}, {"Ra", 226.0}, {"Ac", 227.0}, {"Th", 232.04},
  {"Pa", 231.04}, {"U", 238.03}, {"Np", 237.0}, {"Pu", 244.0}, {"Am", 243.0}, {"Cm", 247.0}
};

static double get_atom_mass(const std::string& element) {
  auto it = atom_masses.find(element);
  return (it != atom_masses.end()) ? it->second : 0.0;
}

// common pseudopotential filename suffixes
static std::unordered_map<std::string, std::string> atom_pseudos = {
  {"H",  "H.pbe-rrkjus_psl.1.0.0.UPF"},
  {"Li", "Li.pbe-s-rrkjus_psl.1.0.0.UPF"},
  {"C",  "C.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"N",  "N.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"O",  "O.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"F",  "F.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"Na", "Na.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Mg", "Mg.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Al", "Al.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"Si", "Si.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"P",  "P.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"S",  "S.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"Cl", "Cl.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"K",  "K.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ca", "Ca.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Sc", "Sc.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ti", "Ti.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"V",  "V.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Cr", "Cr.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Mn", "Mn.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Fe", "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Co", "Co.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ni", "Ni.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Cu", "Cu.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Zn", "Zn.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Ga", "Ga.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Ge", "Ge.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"As", "As.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"Br", "Br.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"Rb", "Rb.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Sr", "Sr.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Y",  "Y.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Zr", "Zr.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Nb", "Nb.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Mo", "Mo.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ru", "Ru.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Rh", "Rh.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Pd", "Pd.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ag", "Ag.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Cd", "Cd.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Sn", "Sn.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"I",  "I.pbe-n-rrkjus_psl.1.0.0.UPF"},
  {"Cs", "Cs.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ba", "Ba.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Ce", "Ce.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"W",  "W.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Pt", "Pt.pbe-spn-rrkjus_psl.1.0.0.UPF"},
  {"Au", "Au.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Pb", "Pb.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"Bi", "Bi.pbe-dn-rrkjus_psl.1.0.0.UPF"},
  {"U",  "U.pbe-spn-rrkjus_psl.1.0.0.UPF"}
};

static std::string get_atom_pseudo(const std::string& element) {
  auto it = atom_pseudos.find(element);
  return (it != atom_pseudos.end()) ? it->second : "";
}

static std::string render_cell(const zinc::core::Structure& s, const std::string& unit, double alat) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(10);
  const auto& mat = s.lattice.matrix();
  double scale = (unit == "bohr") ? 1.0/zinc::core::BOHR_TO_ANGSTROM : (unit == "alat") ? 1.0/alat : 1.0;
  for (int i = 0; i < 3; ++i)
    oss << std::setw(18) << mat(i,0)*scale
        << std::setw(18) << mat(i,1)*scale
        << std::setw(18) << mat(i,2)*scale << "\n";
  return oss.str();
}

static bool has_fixed_atom(const zinc::core::Structure& s) {
  return std::any_of(s.atoms.begin(), s.atoms.end(), [](const auto& atom) {
    return atom.constraints && ((*atom.constraints)[0] == 0 ||
                                (*atom.constraints)[1] == 0 ||
                                (*atom.constraints)[2] == 0);
  });
}

static std::array<int, 3> atom_constraints_or_free(const zinc::core::Atom& atom) {
  return atom.constraints.value_or(std::array<int, 3>{1, 1, 1});
}

static std::string render_pos(const zinc::core::Structure& s, const std::string& unit, double alat) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(10);
  const bool write_constraints = has_fixed_atom(s);
  for (const auto& atom : s.atoms) {
    if (unit == "crystal" || unit == "crystal_sg") {
      auto f = s.lattice.cart_to_frac(atom.position);
      if (f) {
        oss << std::left << std::setw(4) << atom.element << std::right
            << std::setw(18) << f->x()
            << std::setw(18) << f->y()
            << std::setw(18) << f->z();
        if (write_constraints) {
          const auto c = atom_constraints_or_free(atom);
          oss << "   " << c[0] << "   " << c[1] << "   " << c[2];
        }
        oss << "\n";
      }
    } else {
      double scale = (unit == "bohr") ? 1.0/zinc::core::BOHR_TO_ANGSTROM : (unit == "alat") ? 1.0/alat : 1.0;
      oss << std::left << std::setw(4) << atom.element << std::right
          << std::setw(18) << atom.position.x()*scale
          << std::setw(18) << atom.position.y()*scale
          << std::setw(18) << atom.position.z()*scale;
      if (write_constraints) {
        const auto c = atom_constraints_or_free(atom);
        oss << "   " << c[0] << "   " << c[1] << "   " << c[2];
      }
      oss << "\n";
    }
  }
  return oss.str();
}

// Parse lattice from CELL_PARAMETERS block in QE input content
static std::optional<zinc::core::Lattice> parse_lattice_from_qe_input(const std::string& content) {
  std::vector<std::string> lines;
  std::istringstream css(content);
  std::string line;
  while (std::getline(css, line)) lines.push_back(line);

  double celldm1 = 0.0;
  for (const auto& l : lines) {
    std::string lower = l;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("celldm(1)") != std::string::npos) {
      std::istringstream ss(l);
      std::string tok;
      while (ss >> tok) {
        if (tok.find('=') != std::string::npos) {
          auto eq = tok.find('=');
          std::string right = tok.substr(eq + 1);
          // remove trailing comma
          if (!right.empty() && right.back() == ',') right.pop_back();
          try { celldm1 = std::stod(right); break; } catch (...) {}
        } else {
          try { celldm1 = std::stod(tok); break; } catch (...) {}
        }
      }
      break;
    }
  }

  for (size_t i = 0; i < lines.size(); ++i) {
    std::string lt = lines[i];
    size_t ns = lt.find_first_not_of(" \t\r");
    if (ns != std::string::npos) lt = lt.substr(ns);
    std::string lower = lt;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower.find("cell_parameters") == 0) {
      if (i + 3 >= lines.size()) return std::nullopt;

      // determine unit
      zinc::core::LengthUnit unit = zinc::core::LengthUnit::Bohr;
      {
        std::string ul = lower;
        if (ul.find("angstrom") != std::string::npos) unit = zinc::core::LengthUnit::Angstrom;
        else if (ul.find("alat") != std::string::npos) unit = zinc::core::LengthUnit::Alat;
      }

      auto parse_v = [](const std::string& l) -> zinc::core::Vector3d {
        std::istringstream ss(l);
        double x = 0, y = 0, z = 0;
        ss >> x >> y >> z;
        return zinc::core::Vector3d(x, y, z);
      };

      auto v1 = parse_v(lines[i+1]);
      auto v2 = parse_v(lines[i+2]);
      auto v3 = parse_v(lines[i+3]);

      double s = zinc::core::to_angstrom(unit, 1.0,
        (unit == zinc::core::LengthUnit::Alat) ? celldm1 : 1.0);
      return zinc::core::Lattice(v1 * s, v2 * s, v3 * s);
    }
  }
  return std::nullopt;
}

// Parse existing ATOMIC_SPECIES entries: element -> (mass, pseudo)
static std::unordered_map<std::string, std::pair<double, std::string>>
parse_atomic_species(const std::string& content) {
  std::unordered_map<std::string, std::pair<double, std::string>> map;
  std::vector<std::string> lines;
  std::istringstream css(content);
  std::string line;
  while (std::getline(css, line)) lines.push_back(line);

  bool in_block = false;
  for (size_t i = 0; i < lines.size(); ++i) {
    std::string lt = lines[i];
    size_t ns = lt.find_first_not_of(" \t\r");
    if (ns != std::string::npos) lt = lt.substr(ns);
    std::string lower = lt;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower.find("atomic_species") == 0) {
      in_block = true;
      continue;
    }
    if (!in_block) continue;

    if (lt.empty()) break;
    if (lt[0] == '&' || lt[0] == '/' ||
        lower.find("cell_parameters") == 0 ||
        lower.find("atomic_positions") == 0 ||
        lower.find("k_points") == 0) break;

    std::istringstream ss(lt);
    std::string el, mass_str, pseudo;
    if (ss >> el >> mass_str) {
      double mass = parse_fortran(mass_str);
      ss >> pseudo;
      map[el] = {mass > 0 ? mass : get_atom_mass(el), pseudo};
    }
  }
  return map;
}

// Replace ATOMIC_SPECIES block
static std::string replace_atomic_species(const std::string& content,
    const std::set<std::string>& elements, bool fill_pseudo) {
  auto existing = parse_atomic_species(content);

  // build new block
  std::ostringstream new_block;
  new_block << "ATOMIC_SPECIES\n";

  for (const auto& el : elements) {
    double mass;
    std::string pseudo;
    auto it = existing.find(el);
    if (it != existing.end()) {
      mass = it->second.first;
      pseudo = it->second.second;
    } else {
      mass = get_atom_mass(el);
      if (fill_pseudo) pseudo = get_atom_pseudo(el);
    }
    new_block << el << "  " << std::setprecision(6) << mass;
    if (!pseudo.empty()) new_block << "  " << pseudo;
    new_block << "\n";
  }

  std::string block_str = new_block.str();

  // find and replace ATOMIC_SPECIES block
  std::vector<std::string> lines;
  std::istringstream css(content);
  std::string line;
  while (std::getline(css, line)) lines.push_back(line);

  std::ostringstream out;
  bool replaced = false;
  size_t i = 0;
  while (i < lines.size()) {
    std::string lt = lines[i];
    size_t ns = lt.find_first_not_of(" \t\r");
    std::string lt2 = (ns != std::string::npos) ? lt.substr(ns) : lt;
    std::string lower = lt2;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    if (lower.find("atomic_species") == 0) {
      out << block_str;
      replaced = true;
      i++;
      // skip old block lines
      while (i < lines.size()) {
        std::string nl = lines[i];
        ns = nl.find_first_not_of(" \t\r");
        std::string nt = (ns != std::string::npos) ? nl.substr(ns) : nl;
        if (nt.empty()) break;
        std::string nlower = nt;
        std::transform(nlower.begin(), nlower.end(), nlower.begin(), ::tolower);
        if (nt[0] == '&' || nt[0] == '/' ||
            nlower.find("cell_parameters") == 0 ||
            nlower.find("atomic_positions") == 0 ||
            nlower.find("k_points") == 0) break;
        i++;
      }
      continue;
    }
    out << lines[i] << "\n";
    i++;
  }

  if (!replaced) {
    // Insert before ATOMIC_POSITIONS if present, else before CELL_PARAMETERS, else append.
    std::string out_str = out.str();
    css.clear();
    css.str(out_str);
    lines.clear();
    while (std::getline(css, line)) lines.push_back(line);

    out.str("");
    out.clear();
    bool inserted = false;
    for (const auto& ol : lines) {
      std::string tl = ol;
      size_t tl_ns = tl.find_first_not_of(" \t\r");
      if (tl_ns != std::string::npos) tl = tl.substr(tl_ns);
      std::string tlower = tl;
      std::transform(tlower.begin(), tlower.end(), tlower.begin(), ::tolower);
      if (!inserted &&
          (tlower.find("atomic_positions") == 0 || tlower.find("cell_parameters") == 0)) {
        out << block_str;
        inserted = true;
      }
      out << ol << "\n";
    }
    if (!inserted) {
      out << "\n" << block_str;
    }
  }

  return out.str();
}

void run_update(zinc::io::DriverManager& manager, const std::string& source, const std::string& target,
                const fs::path* output, bool in_place, const std::string& pos_unit, const std::string& cell_unit,
                bool fill_pseudo, const std::string& report_format) {
  // validate
  static const std::set<std::string> pos_units = {"angstrom", "bohr", "alat", "crystal", "crystal_sg"};
  static const std::set<std::string> cell_units = {"angstrom", "bohr", "alat"};
  if (source == "-" && target == "-")
    throw std::runtime_error("Cannot read both source and target from stdin");
  if (in_place && output)
    throw std::invalid_argument("--in-place cannot be used together with --output");
  if (in_place && target == "-")
    throw std::invalid_argument("--in-place requires a real target file (target cannot be '-')");
  if (pos_units.find(pos_unit) == pos_units.end())
    throw std::invalid_argument("unsupported --pos-unit: " + pos_unit);
  if (cell_units.find(cell_unit) == cell_units.end())
    throw std::invalid_argument("unsupported --cell-unit: " + cell_unit);

  // load source structure (handle stdin via temp file)
  fs::path source_path;
  if (source == "-") {
    std::string buf((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
    source_path = fs::temp_directory_path() / ("zinc-stdin-" + std::to_string(
      std::chrono::steady_clock::now().time_since_epoch().count()) + ".qe");
    std::ofstream tmp(source_path);
    if (!tmp.is_open()) throw std::runtime_error("Cannot create temporary input file: " + source_path.string());
    tmp << buf;
  } else {
    source_path = source;
  }
  auto struct_ = manager.load_auto(source_path);
  if (source == "-") std::filesystem::remove(source_path);

  // load target content
  std::string content;
  if (target == "-") {
    content = std::string((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
  } else {
    std::ifstream tgt(target);
    if (!tgt.is_open()) throw std::runtime_error("Cannot open target: " + target);
    content = std::string((std::istreambuf_iterator<char>(tgt)), std::istreambuf_iterator<char>());
  }

  // inherit lattice if source missing
  if (struct_.metadata["lattice_missing"] == "true") {
    auto inherited = parse_lattice_from_qe_input(content);
    if (inherited) {
      struct_.lattice = *inherited;
    } else {
      throw std::runtime_error("Source lattice missing and target provides no CELL_PARAMETERS block to inherit");
    }
  }

  double alat = struct_.lattice.matrix().row(0).norm();
  if (alat < 1e-12) alat = 1.0;

  std::string cell_header = "CELL_PARAMETERS (" + cell_unit + ")\n";
  std::string pos_header = "ATOMIC_POSITIONS (" + pos_unit + ")\n";

  // Replace CELL_PARAMETERS block
  {
    std::vector<std::string> lines;
    std::istringstream css(content);
    std::string line;
    while (std::getline(css, line)) lines.push_back(line);

    std::ostringstream out;
    bool replaced = false;
    size_t i = 0;
    while (i < lines.size()) {
      std::string lt = lines[i];
      size_t ns = lt.find_first_not_of(" \t\r");
      std::string lt2 = (ns != std::string::npos) ? lt.substr(ns) : lt;
      std::string lower = lt2;
      std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

      if (lower.find("cell_parameters") == 0) {
        out << cell_header << render_cell(struct_, cell_unit, alat);
        replaced = true;
        i += 4;
        continue;
      }
      out << lines[i] << "\n";
      i++;
    }
    if (!replaced) {
      out << cell_header << render_cell(struct_, cell_unit, alat);
    }
    content = out.str();
  }

  // Replace ATOMIC_POSITIONS block
  {
    std::vector<std::string> lines;
    std::istringstream css(content);
    std::string line;
    while (std::getline(css, line)) lines.push_back(line);

    std::ostringstream out;
    bool replaced = false;
    int nat = static_cast<int>(struct_.natoms());
    size_t i = 0;
    while (i < lines.size()) {
      std::string lt = lines[i];
      size_t ns = lt.find_first_not_of(" \t\r");
      std::string lt2 = (ns != std::string::npos) ? lt.substr(ns) : lt;
      std::string lower = lt2;
      std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

      if (lower.find("atomic_positions") == 0) {
        out << pos_header << render_pos(struct_, pos_unit, alat);
        replaced = true;
        i++;
        int skipped = 0;
        while (i < lines.size() && skipped < nat) {
          std::string nl = lines[i];
          size_t ns2 = nl.find_first_not_of(" \t\r");
          std::string nt = (ns2 != std::string::npos) ? nl.substr(ns2) : nl;
          if (nt.empty()) { i++; continue; }
          std::string nlower = nt;
          std::transform(nlower.begin(), nlower.end(), nlower.begin(), ::tolower);
          if (nt[0] == '&' || nt[0] == '/' ||
              nlower.find("cell_parameters") == 0 ||
              nlower.find("atomic_species") == 0 ||
              nlower.find("k_points") == 0) break;
          if (std::isalpha(static_cast<unsigned char>(nt[0]))) skipped++;
          i++;
        }
        continue;
      }
      out << lines[i] << "\n";
      i++;
    }
    if (!replaced) {
      out << "\n" << pos_header << render_pos(struct_, pos_unit, alat);
    }
    content = out.str();
  }

  // Replace ATOMIC_SPECIES block
  std::set<std::string> elems;
  for (const auto& atom : struct_.atoms) elems.insert(atom.element);
  content = replace_atomic_species(content, elems, fill_pseudo);

  // Update nat/ntyp
  {
    auto p = content.find("&SYSTEM");
    if (p != std::string::npos) {
      auto e = content.find('/', p);
      if (e != std::string::npos) {
        std::string block = content.substr(p, e - p + 1);
        auto replace_val = [](std::string& b, const std::string& key, int val) {
          auto kp = b.find(key);
          if (kp == std::string::npos) return;
          auto eq = b.find('=', kp);
          if (eq == std::string::npos) return;
          auto vs = b.find_first_of("0123456789", eq);
          if (vs == std::string::npos) return;
          auto ve = b.find_first_not_of("0123456789.", vs);
          if (ve == std::string::npos) ve = b.length();
          b.replace(vs, ve - vs, std::to_string(val));
        };
        replace_val(block, "nat", static_cast<int>(struct_.natoms()));
        replace_val(block, "ntyp", static_cast<int>(elems.size()));
        content = content.substr(0, p) + block + content.substr(e + 1);
      }
    }
  }

  // write output
  std::string written_to;
  std::string mode;
  if (in_place) {
    std::ofstream ofs(target);
    if (!ofs.is_open()) throw std::runtime_error("Cannot open output file: " + target);
    ofs << content;
    written_to = target;
    mode = "in-place";
  } else if (output) {
    std::ofstream ofs(*output);
    if (!ofs.is_open()) throw std::runtime_error("Cannot open output file: " + output->string());
    ofs << content;
    written_to = output->string();
    mode = "file";
  } else {
    fmt::print("{}", content);
    return;
  }

  if (report_format == "json" || report_format == "jsonl") {
    const auto json = fmt::format(
        "{{\"schema_version\":\"1.0\",\"command\":\"update\",\"source\":\"{}\",\"target\":\"{}\",\"output\":\"{}\",\"mode\":\"{}\",\"natoms\":{},\"ntyp\":{},\"status\":\"ok\"}}",
        json_escape(source),
        json_escape(target),
        json_escape(written_to),
        json_escape(mode),
        struct_.natoms(),
        elems.size());
    if (report_format == "jsonl") {
      fmt::println("{}", json);
    } else {
      fmt::print("{}", json);
    }
  } else if (in_place) {
    fmt::println("Updated {} in-place", target);
  } else {
    fmt::println("Written to {}", written_to);
  }
}
