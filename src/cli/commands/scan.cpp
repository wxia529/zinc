#include <fmt/core.h>
#include <fmt/ranges.h>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <set>
#include <string>
#include <algorithm>
#include <cctype>
#include "zinc/io/driver.h"
#include "zinc/core/structure.h"

namespace fs = std::filesystem;

static std::string get_composition(const zinc::core::Structure& s) {
  std::unordered_map<std::string, int> counts;
  for (const auto& atom : s.atoms) {
    counts[atom.element]++;
  }
  std::set<std::string> elems;
  for (const auto& [e, _] : counts) elems.insert(e);

  std::string res;
  for (const auto& e : elems) {
    res += e + std::to_string(counts[e]);
  }
  return res;
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

static bool is_structure_file(const fs::path& path) {
  std::string name = path.filename().string();
  std::string ext = path.extension().string();
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  return name == "poscar" || name == "contcar" ||
      ext == ".cif" || ext == ".xyz" || ext == ".extxyz" ||
      ext == ".vasp" || ext == ".poscar" || ext == ".contcar" ||
      ext == ".xml" || ext == ".out" || ext == ".log" ||
      ext == ".in" || ext == ".pwi";
}

int run_scan(
    zinc::io::DriverManager& manager,
    const std::filesystem::path& path,
    const std::string& filter_element,
    const std::string& output_format) {
  if (output_format == "text") {
    fmt::println("{:<40} {:<10} {:>5}  {}", "File", "Type", "Atoms", "Composition");
    fmt::println("{:-<40} {:-<10} {:->5}  {:-<15}", "", "", "", "");
  }

  std::size_t ok_count = 0;
  std::size_t fail_count = 0;

  for (const auto& entry : fs::recursive_directory_iterator(path)) {
    if (!entry.is_regular_file()) continue;
    auto p = entry.path();
    if (!is_structure_file(p)) continue;
    try {
      auto loaded = manager.load_with_detection(p);
      auto& struct_ = loaded.structure;
      if (!filter_element.empty()) {
        bool has_elem = false;
        for (const auto& atom : struct_.atoms) {
          if (atom.element == filter_element) { has_elem = true; break; }
        }
        if (!has_elem) continue;
      }
      std::string display_path = fs::relative(p, path).string();
      ok_count++;
      if (output_format == "jsonl") {
        fmt::println(
            "{{\"schema_version\":\"1.0\",\"command\":\"scan\",\"input\":\"{}\",\"detected_format\":\"{}\",\"natoms\":{},\"cell_volume\":{},\"status\":\"ok\"}}",
            json_escape(display_path),
            json_escape(loaded.detected_format),
            struct_.natoms(),
            struct_.lattice.volume());
      } else {
        if (display_path.length() > 40) display_path = "..." + display_path.substr(display_path.length() - 37);
        fmt::println("{:<40} {:<10} {:>5}  {}", display_path, loaded.detected_format, struct_.natoms(), get_composition(struct_));
      }
    } catch (const std::exception& e) {
      fail_count++;
      std::string display_path = fs::relative(p, path).string();
      fmt::println(stderr, "[scan] parse failed: {}", display_path);
      fmt::println(stderr, "reason: {}", e.what());
      fmt::println(stderr, "hint: try --from qe-log or verify file completeness");
      if (output_format == "jsonl") {
        fmt::println(
            "{{\"schema_version\":\"1.0\",\"command\":\"scan\",\"input\":\"{}\",\"status\":\"failed\",\"error\":{{\"code\":\"PARSE_ERROR\",\"message\":\"{}\"}}}}",
            json_escape(display_path),
            json_escape(e.what()));
      }
    }
  }

  if (ok_count > 0 && fail_count > 0) {
    return 2;
  }
  if (ok_count == 0 && fail_count > 0) {
    return 1;
  }
  return 0;
}
