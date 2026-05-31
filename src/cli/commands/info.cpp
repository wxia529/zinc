#include <fmt/core.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include "zinc/io/driver.h"

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

static std::string build_inspect_json(
    const std::filesystem::path& file,
    const zinc::io::DriverManager::LoadResult& loaded) {
  const auto& structure = loaded.structure;
  std::ostringstream meta;
  bool first = true;
  for (const auto& kv : structure.metadata) {
    if (kv.first == "trajectory_extxyz") {
      continue;
    }
    if (!first) {
      meta << ",";
    }
    first = false;
    meta << "\"" << json_escape(kv.first) << "\":\"" << json_escape(kv.second) << "\"";
  }

  std::ostringstream out;
  out << std::setprecision(16);
  out << "{"
      << "\"schema_version\":\"1.0\","
      << "\"command\":\"info\","
      << "\"input\":\"" << json_escape(file.string()) << "\","
      << "\"detected_format\":\"" << json_escape(loaded.detected_format) << "\","
      << "\"natoms\":" << structure.natoms() << ","
      << "\"cell_volume\":" << structure.lattice.volume() << ","
      << "\"metadata\":{" << meta.str() << "}"
      << "}";
  return out.str();
}

int run_info(
    zinc::io::DriverManager& manager,
    const std::filesystem::path& file,
    const std::string& from_format,
    const std::string& output_format,
    const std::filesystem::path* output) {
  const auto loaded = manager.load_with_detection(file, from_format);
  const auto& structure = loaded.structure;

  std::string rendered;
  if (output_format == "json" || output_format == "jsonl") {
    rendered = build_inspect_json(file, loaded);
    if (output_format == "jsonl") {
      rendered.push_back('\n');
    }
  } else {
    std::ostringstream oss;
    oss << "File: " << file.string() << "\n";
    oss << "--------------------------------\n";
    oss << "Detected Format: " << loaded.detected_format << "\n";
    oss << "Number of Atoms: " << structure.natoms() << "\n";
    oss << std::setprecision(8) << "Cell Volume    : " << structure.lattice.volume() << " A^3\n";
    rendered = oss.str();
  }

  if (output) {
    std::ofstream ofs(*output);
    if (!ofs.is_open()) {
      throw std::runtime_error("Cannot open output file: " + output->string());
    }
    ofs << rendered;
  } else {
    fmt::print("{}", rendered);
  }
  return 0;
}
