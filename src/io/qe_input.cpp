#include "zinc/io/qe_input.h"
#include "zinc/core/constants.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace zinc {
namespace io {

namespace {

enum class QEUnit { Angstrom, Bohr, Alat, Crystal };

static std::string lower_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

static std::string trim_copy(const std::string& s) {
  const auto begin = s.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return "";
  const auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(begin, end - begin + 1);
}

static std::string clean_line(const std::string& line) {
  auto end = line.find_first_of("!#");
  return trim_copy(line.substr(0, end));
}

static double parse_number(std::string token) {
  std::replace(token.begin(), token.end(), 'D', 'E');
  std::replace(token.begin(), token.end(), 'd', 'e');
  auto slash = token.find('/');
  if (slash != std::string::npos) {
    const double left = parse_number(token.substr(0, slash));
    const double right = parse_number(token.substr(slash + 1));
    if (std::abs(right) < 1e-12) return 0.0;
    return left / right;
  }
  return std::stod(token);
}

static std::optional<double> parse_assignment_number(const std::string& content, const std::string& key) {
  const std::regex re(key + R"(\s*=\s*['"]?([+-]?\d*\.?\d+(?:[eEdD][+-]?\d+)?))",
                      std::regex_constants::icase);
  std::smatch m;
  if (!std::regex_search(content, m, re)) return std::nullopt;
  try { return parse_number(m[1].str()); } catch (...) { return std::nullopt; }
}

static std::optional<int> parse_assignment_int(const std::string& content, const std::string& key) {
  const std::regex re(key + R"(\s*=\s*(-?\d+))", std::regex_constants::icase);
  std::smatch m;
  if (!std::regex_search(content, m, re)) return std::nullopt;
  try { return std::stoi(m[1].str()); } catch (...) { return std::nullopt; }
}

static std::optional<double> extract_alat_angstrom(const std::vector<std::string>& lines) {
  const std::string content = [&]() {
    std::ostringstream oss;
    for (const auto& line : lines) oss << clean_line(line) << "\n";
    return oss.str();
  }();

  if (auto v = parse_assignment_number(content, R"(celldm\s*\(\s*1\s*\))")) {
    return *v * core::BOHR_TO_ANGSTROM;
  }
  if (auto v = parse_assignment_number(content, R"(\bA\b)")) {
    return *v;
  }
  return std::nullopt;
}

static QEUnit detect_cell_unit(const std::string& header) {
  auto lower = lower_copy(header);
  if (lower.find("angstrom") != std::string::npos) return QEUnit::Angstrom;
  if (lower.find("alat") != std::string::npos) return QEUnit::Alat;
  return QEUnit::Bohr;
}

static QEUnit detect_position_unit(const std::string& header) {
  auto lower = lower_copy(header);
  if (lower.find("crystal") != std::string::npos) return QEUnit::Crystal;
  if (lower.find("bohr") != std::string::npos) return QEUnit::Bohr;
  if (lower.find("alat") != std::string::npos) return QEUnit::Alat;
  return QEUnit::Angstrom;
}

static double unit_scale_to_angstrom(QEUnit unit, double alat_angstrom) {
  if (unit == QEUnit::Angstrom || unit == QEUnit::Crystal) return 1.0;
  if (unit == QEUnit::Bohr) return core::BOHR_TO_ANGSTROM;
  return alat_angstrom;
}

static std::string label_to_symbol(const std::string& label) {
  std::string letters;
  for (char ch : label) {
    if (std::isalpha(static_cast<unsigned char>(ch))) {
      letters.push_back(ch);
    } else if (!letters.empty()) {
      break;
    }
  }
  if (letters.empty()) return "X";
  letters[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0])));
  if (letters.size() > 1) {
    letters[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1])));
    letters.resize(2);
  }
  return letters;
}

static std::optional<std::array<double, 3>> parse_vec3(const std::string& line) {
  std::istringstream ss(line);
  std::string sx, sy, sz;
  if (!(ss >> sx >> sy >> sz)) return std::nullopt;
  try {
    return std::array<double, 3>{parse_number(sx), parse_number(sy), parse_number(sz)};
  } catch (...) {
    return std::nullopt;
  }
}

static std::optional<core::Lattice> parse_cell(
    const std::vector<std::string>& lines,
    std::size_t start,
    double alat_angstrom) {
  const auto unit = detect_cell_unit(lines[start]);
  const double scale = unit_scale_to_angstrom(unit, alat_angstrom);
  std::vector<core::Vector3d> vecs;

  for (std::size_t i = start + 1; i < lines.size() && vecs.size() < 3; ++i) {
    const auto cleaned = clean_line(lines[i]);
    if (cleaned.empty()) continue;
    auto parsed = parse_vec3(cleaned);
    if (!parsed) break;
    vecs.emplace_back((*parsed)[0] * scale, (*parsed)[1] * scale, (*parsed)[2] * scale);
  }

  if (vecs.size() != 3) return std::nullopt;
  return core::Lattice(vecs[0], vecs[1], vecs[2]);
}

static bool is_block_start(const std::string& cleaned) {
  const auto lower = lower_copy(cleaned);
  return cleaned.empty() ||
         cleaned[0] == '&' ||
         cleaned[0] == '/' ||
         lower.find("atomic_species") == 0 ||
         lower.find("atomic_positions") == 0 ||
         lower.find("cell_parameters") == 0 ||
         lower.find("k_points") == 0;
}

} // namespace

bool QEInputDriver::identify(const std::string& first_line) const {
  const auto lower = lower_copy(first_line);
  if (lower.find("program pwscf") != std::string::npos ||
      lower.find("quantum espresso") != std::string::npos ||
      lower.find("<") != std::string::npos) {
    return false;
  }
  return lower.find("&control") != std::string::npos ||
         lower.find("&system") != std::string::npos ||
         (lower.find("atomic_species") != std::string::npos &&
          lower.find("atomic_positions") != std::string::npos);
}

core::Structure QEInputDriver::read(const std::filesystem::path& path) const {
  std::ifstream file(path);
  if (!file.is_open()) throw std::runtime_error("Cannot open: " + path.string());

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(file, line)) lines.push_back(line);

  const std::string content = [&]() {
    std::ostringstream oss;
    for (const auto& l : lines) oss << clean_line(l) << "\n";
    return oss.str();
  }();

  if (auto ibrav = parse_assignment_int(content, R"(\bibrav\b)")) {
    if (*ibrav != 0) {
      throw std::runtime_error("QE input with ibrav != 0 is not supported; use explicit CELL_PARAMETERS");
    }
  }

  const auto nat = parse_assignment_int(content, R"(\bnat\b)");
  const double alat_angstrom = extract_alat_angstrom(lines).value_or(1.0);

  std::optional<core::Lattice> lattice;
  std::size_t positions_start = lines.size();
  for (std::size_t i = 0; i < lines.size(); ++i) {
    const auto cleaned = clean_line(lines[i]);
    const auto lower = lower_copy(cleaned);
    if (lower.find("cell_parameters") == 0) {
      lattice = parse_cell(lines, i, alat_angstrom);
    } else if (lower.find("atomic_positions") == 0) {
      positions_start = i;
    }
  }

  if (positions_start == lines.size()) {
    throw std::runtime_error("Invalid QE input: missing ATOMIC_POSITIONS");
  }

  core::Structure structure(lattice.value_or(core::Lattice()));
  structure.metadata["detected_format_detail"] = "qe-input";
  structure.metadata["lattice_missing"] = lattice ? "false" : "true";

  const auto pos_unit = detect_position_unit(lines[positions_start]);
  const double pos_scale = unit_scale_to_angstrom(pos_unit, alat_angstrom);
  int parsed_atoms = 0;

  for (std::size_t i = positions_start + 1; i < lines.size(); ++i) {
    const auto cleaned = clean_line(lines[i]);
    if (is_block_start(cleaned)) {
      if (parsed_atoms > 0) break;
      continue;
    }

    std::istringstream ss(cleaned);
    std::string label, sx, sy, sz;
    if (!(ss >> label >> sx >> sy >> sz)) {
      if (parsed_atoms > 0) break;
      continue;
    }

    double x = 0.0, y = 0.0, z = 0.0;
    try {
      x = parse_number(sx);
      y = parse_number(sy);
      z = parse_number(sz);
    } catch (...) {
      if (parsed_atoms > 0) break;
      continue;
    }

    std::optional<std::array<int, 3>> constraints;
    std::string c1, c2, c3;
    if (ss >> c1 >> c2 >> c3) {
      try {
        constraints = std::array<int, 3>{std::stoi(c1), std::stoi(c2), std::stoi(c3)};
      } catch (...) {}
    }

    core::Vector3d cart(x, y, z);
    if (pos_unit == QEUnit::Crystal) {
      if (!lattice) {
        throw std::runtime_error("QE input uses crystal coordinates but no CELL_PARAMETERS were found");
      }
      cart = structure.lattice.frac_to_cart(cart);
    } else {
      cart *= pos_scale;
    }

    structure.add_atom(label_to_symbol(label), {cart.x(), cart.y(), cart.z()}, constraints);
    parsed_atoms++;
    if (nat && parsed_atoms >= *nat) break;
  }

  if (structure.natoms() == 0) {
    throw std::runtime_error("Invalid QE input: no atoms parsed from ATOMIC_POSITIONS");
  }

  return structure;
}

void QEInputDriver::write(const core::Structure&, const std::filesystem::path&) const {
  throw std::runtime_error("QE input writing is not supported");
}

} // namespace io
} // namespace zinc
