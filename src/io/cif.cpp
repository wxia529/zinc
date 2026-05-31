#include "zinc/io/cif.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cctype>
#include <optional>

namespace zinc {
namespace io {

bool CifDriver::identify(const std::string& first_line) const {
  std::string lower = first_line;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  return lower.find("data_") != std::string::npos ||
         lower.find("_cell_length_a") != std::string::npos ||
         (lower.find("loop_") != std::string::npos && lower.find("_atom_site_") != std::string::npos);
}

// extract number after a cif tag (case-insensitive)
static double extract_val_ci(const std::string& content, const std::string& key) {
  std::string lower = content;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  std::string lkey = key;
  std::transform(lkey.begin(), lkey.end(), lkey.begin(), ::tolower);

  auto pos = lower.find(lkey);
  if (pos == std::string::npos) return 0.0;
  pos += lkey.length();
  while (pos < content.length() && (content[pos] == ' ' || content[pos] == '\t')) pos++;

  // read number including scientific notation
  std::string num;
  while (pos < content.length() &&
         (std::isdigit(static_cast<unsigned char>(content[pos])) ||
          content[pos] == '.' || content[pos] == '-' ||
          content[pos] == '+' || content[pos] == 'e' || content[pos] == 'E')) {
    num += content[pos++];
  }
  try { return num.empty() ? 0.0 : std::stod(num); } catch (...) { return 0.0; }
}

// extract element symbol from label like "C1", "Fe2", "H_1"
static std::string extract_element_from_label(const std::string& label) {
  // strip surrounding quotes first
  std::string s = label;
  while (!s.empty() && (s.front() == '\'' || s.front() == '"')) s.erase(0,1);
  while (!s.empty() && (s.back() == '\'' || s.back() == '"')) s.pop_back();
  if (s.empty()) return "";
  if (std::isupper(static_cast<unsigned char>(s[0]))) {
    std::string el;
    el += s[0];
    if (s.size() > 1 && std::islower(static_cast<unsigned char>(s[1])))
      el += s[1];
    return el;
  }
  return "";
}

// CIF number parse: skip '?' and '.', strip (esd) suffix
static std::optional<double> parse_number_cif(const std::string& raw) {
  std::string s = raw;
  while (!s.empty() && (s.front() == '\'' || s.front() == '"')) s.erase(0,1);
  while (!s.empty() && (s.back() == '\'' || s.back() == '"')) s.pop_back();
  if (s == "?" || s == "." || s.empty()) return std::nullopt;
  auto paren = s.find('(');
  if (paren != std::string::npos) s = s.substr(0, paren);
  try { return std::stod(s); } catch (...) { return std::nullopt; }
}

static std::string strip_quotes(const std::string& s) {
  if (s.size() >= 2 &&
      ((s.front() == '\'' && s.back() == '\'') || (s.front() == '"' && s.back() == '"')))
    return s.substr(1, s.size() - 2);
  return s;
}

static std::vector<std::string> tokenize_line(const std::string& line) {
  std::vector<std::string> tokens;
  std::string current;
  char quote = 0;

  for (size_t i = 0; i < line.size(); ++i) {
    char ch = line[i];
    if (quote) {
      if (ch == quote) { quote = 0; }
      else { current += ch; }
      continue;
    }
    if (ch == '\'' || ch == '"') { quote = ch; continue; }
    if (std::isspace(static_cast<unsigned char>(ch))) {
      if (!current.empty()) { tokens.push_back(current); current.clear(); }
    } else {
      current += ch;
    }
  }
  if (!current.empty()) tokens.push_back(current);
  return tokens;
}

core::Structure CifDriver::read(const std::filesystem::path& path) const {
  std::ifstream file(path);
  std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

  double a = extract_val_ci(content, "_cell_length_a");
  double b = extract_val_ci(content, "_cell_length_b");
  double c = extract_val_ci(content, "_cell_length_c");
  double alpha = extract_val_ci(content, "_cell_angle_alpha");
  double beta = extract_val_ci(content, "_cell_angle_beta");
  double gamma = extract_val_ci(content, "_cell_angle_gamma");
  if (a < 1e-6) a = 1.0; if (b < 1e-6) b = 1.0; if (c < 1e-6) c = 1.0;

  core::Lattice lattice = core::Lattice::from_cell_parameters({a, b, c}, {alpha, beta, gamma});
  core::Structure structure(lattice);

  // split into lines for loop parsing
  std::vector<std::string> lines;
  std::istringstream css(content);
  std::string line;
  while (std::getline(css, line)) lines.push_back(line);

  // parse loose _atom_site_* columns
  std::unordered_map<std::string, std::vector<std::string>> loose_cols;

  size_t li = 0;
  while (li < lines.size()) {
    std::string raw = lines[li];
    // trim
    size_t ns = raw.find_first_not_of(" \t\r");
    std::string trimmed = (ns != std::string::npos) ? raw.substr(ns) : "";
    if (trimmed.empty() || trimmed[0] == '#') { li++; continue; }

    std::string lower = trimmed;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    // loop_ block
    if (lower == "loop_") {
      // read headers
      std::vector<std::string> headers;
      size_t hi = li + 1;
      while (hi < lines.size()) {
        std::string hraw = lines[hi];
        size_t hns = hraw.find_first_not_of(" \t\r");
        std::string ht = (hns != std::string::npos) ? hraw.substr(hns) : "";
        if (ht.empty() || ht[0] == '#') { hi++; continue; }
        if (ht[0] != '_') break;
        std::string hl = ht;
        std::transform(hl.begin(), hl.end(), hl.begin(), ::tolower);
        size_t he = ht.find_first_of(" \t");
        if (he != std::string::npos) ht = ht.substr(0, he);
        headers.push_back(hl);
        hi++;
      }

      // read data rows
      size_t ri = hi;
      std::vector<std::string> row_buffer;
      while (ri < lines.size()) {
        std::string rraw = lines[ri];
        // handle ';' multiline value
        if (rraw.size() > 0 && rraw[0] == ';') {
          std::string multiline = rraw.substr(1);
          ri++;
          while (ri < lines.size() && lines[ri] != ";") {
            multiline += "\n" + lines[ri];
            ri++;
          }
          ri++; // skip closing ';'
          row_buffer.push_back(multiline);
          continue;
        }

        size_t rns = rraw.find_first_not_of(" \t\r");
        if (rns == std::string::npos) { ri++; continue; }
        std::string rt = rraw.substr(rns);
        if (rt[0] == '#') { ri++; continue; }
        // strip inline comment
        {
          bool in_quote = false;
          size_t cp = 0;
          for (; cp < rt.size(); ++cp) {
            if (rt[cp] == '\'' || rt[cp] == '"') in_quote = !in_quote;
            else if (!in_quote && rt[cp] == '#') break;
          }
          if (cp < rt.size()) rt = rt.substr(0, cp);
          // trim trailing space after stripping comment
          while (!rt.empty() && (rt.back() == ' ' || rt.back() == '\t')) rt.pop_back();
        }
        if (rt.empty()) { ri++; continue; }

        std::string rl = rt;
        std::transform(rl.begin(), rl.end(), rl.begin(), ::tolower);
        if (rt[0] == '_' || rl == "loop_" || rl.find("data_") == 0) break;

        auto tokens = tokenize_line(rt);
        for (auto& t : tokens) row_buffer.push_back(t);

        while (!headers.empty() && row_buffer.size() >= headers.size()) {
          std::vector<std::string> row(row_buffer.begin(), row_buffer.begin() + headers.size());
          row_buffer.erase(row_buffer.begin(), row_buffer.begin() + headers.size());

          std::string sym;
          std::optional<double> ofx, ofy, ofz, ocx, ocy, ocz;

          for (size_t hi2 = 0; hi2 < headers.size() && hi2 < row.size(); ++hi2) {
            const auto& h = headers[hi2];
            if (h == "_atom_site_type_symbol") { sym = strip_quotes(row[hi2]); if (sym == "?" || sym == ".") sym = ""; }
            else if (h == "_atom_site_label" && sym.empty()) sym = extract_element_from_label(row[hi2]);
            else if (h == "_atom_site_fract_x") ofx = parse_number_cif(row[hi2]);
            else if (h == "_atom_site_fract_y") ofy = parse_number_cif(row[hi2]);
            else if (h == "_atom_site_fract_z") ofz = parse_number_cif(row[hi2]);
            else if (h == "_atom_site_cartn_x") ocx = parse_number_cif(row[hi2]);
            else if (h == "_atom_site_cartn_y") ocy = parse_number_cif(row[hi2]);
            else if (h == "_atom_site_cartn_z") ocz = parse_number_cif(row[hi2]);
          }

          if (sym.empty()) sym = "X";
          if (ofx && ofy && ofz) {
            core::Vector3d cart = lattice.frac_to_cart(core::Vector3d(*ofx, *ofy, *ofz));
            structure.add_atom(sym, {cart.x(), cart.y(), cart.z()});
          } else if (ocx && ocy && ocz) {
            structure.add_atom(sym, {*ocx, *ocy, *ocz});
          }
        }
        ri++;
      }
      li = ri;
      continue;
    }

    // loose _atom_site_* column
    if (lower.find("_atom_site_") == 0) {
      size_t sp = trimmed.find_first_of(" \t");
      std::string key;
      std::string val;
      if (sp != std::string::npos) {
        key = trimmed.substr(0, sp);
        val = trimmed.substr(sp);
        // trim val
        size_t vs = val.find_first_not_of(" \t");
        if (vs != std::string::npos) val = val.substr(vs);
        size_t ve = val.find_last_not_of(" \t\r");
        if (ve != std::string::npos) val = val.substr(0, ve + 1);
      }
      std::string kl = key;
      std::transform(kl.begin(), kl.end(), kl.begin(), ::tolower);
      if (!kl.empty() && !val.empty()) {
        loose_cols[kl].push_back(val);
      }
    }

    li++;
  }

  // build atoms from loose columns if no atoms found from loop_ blocks
  if (structure.natoms() == 0 && !loose_cols.empty()) {
    // determine max length
    size_t len = 0;
    for (const auto& [k, v] : loose_cols) {
      if (v.size() > len) len = v.size();
    }

    for (size_t idx = 0; idx < len; ++idx) {
      std::string sym;
      if (loose_cols.count("_atom_site_type_symbol") && idx < loose_cols["_atom_site_type_symbol"].size())
        sym = strip_quotes(loose_cols["_atom_site_type_symbol"][idx]);
      if (sym.empty() && loose_cols.count("_atom_site_label") && idx < loose_cols["_atom_site_label"].size())
        sym = extract_element_from_label(loose_cols["_atom_site_label"][idx]);
      if (sym.empty()) sym = "X";

      auto getv = [&](const std::string& k) -> std::optional<double> {
        if (loose_cols.count(k) && idx < loose_cols[k].size()) {
          return parse_number_cif(loose_cols[k][idx]);
        }
        return std::nullopt;
      };

      auto ofx = getv("_atom_site_fract_x");
      auto ofy = getv("_atom_site_fract_y");
      auto ofz = getv("_atom_site_fract_z");
      auto ocx = getv("_atom_site_cartn_x");
      auto ocy = getv("_atom_site_cartn_y");
      auto ocz = getv("_atom_site_cartn_z");

      if (ofx && ofy && ofz) {
        core::Vector3d cart = lattice.frac_to_cart(core::Vector3d(*ofx, *ofy, *ofz));
        structure.add_atom(sym, {cart.x(), cart.y(), cart.z()});
      } else if (ocx && ocy && ocz) {
        structure.add_atom(sym, {*ocx, *ocy, *ocz});
      }
    }
  }

  return structure;
}

void CifDriver::write(const core::Structure& structure, const std::filesystem::path& path) const {
  std::ofstream out(path);
  if (!out.is_open()) throw std::runtime_error("Cannot open output file: " + path.string());
  auto params = structure.lattice.cell_parameters();
  if (!params) throw std::runtime_error("Cannot write CIF: lattice is degenerate");
  auto [lengths, angles] = *params;

  out << "data_zinc\n" << std::fixed << std::setprecision(10);
  out << "_cell_length_a    " << std::setw(18) << lengths[0] << "\n";
  out << "_cell_length_b    " << std::setw(18) << lengths[1] << "\n";
  out << "_cell_length_c    " << std::setw(18) << lengths[2] << "\n";
  out << "_cell_angle_alpha " << std::setw(18) << angles[0] << "\n";
  out << "_cell_angle_beta  " << std::setw(18) << angles[1] << "\n";
  out << "_cell_angle_gamma " << std::setw(18) << angles[2] << "\n";
  out << "_space_group_name_P1   P1\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n";

  for (size_t i = 0; i < structure.atoms.size(); ++i) {
    auto frac = structure.lattice.cart_to_frac(structure.atoms[i].position);
    if (frac) {
      out << std::left << std::setw(8) << (structure.atoms[i].element + std::to_string(i + 1))
          << std::setw(4) << structure.atoms[i].element << std::right
          << std::setw(18) << frac->x()
          << std::setw(18) << frac->y()
          << std::setw(18) << frac->z() << "\n";
    }
  }
}

} // namespace io
} // namespace zinc
