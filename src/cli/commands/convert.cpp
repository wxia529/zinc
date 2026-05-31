#include <fmt/core.h>
#include <filesystem>
#include <fstream>
#include <optional>
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

static std::size_t count_extxyz_frames(const std::string& text) {
  std::size_t count = 0;
  std::size_t pos = 0;
  while ((pos = text.find("Lattice=\"", pos)) != std::string::npos) {
    ++count;
    pos += 9;
  }
  return count;
}

int run_convert(
    zinc::io::DriverManager& manager,
    const std::filesystem::path& input,
    const std::filesystem::path& output,
    const std::string& from_format,
    const std::string& to_format,
    const std::string& output_format) {
  const auto loaded = manager.load_with_detection(input, from_format);
  std::size_t frame_count = 1;
  if (to_format == "extxyz") {
    auto traj = loaded.structure.metadata.find("trajectory_extxyz");
    if (traj != loaded.structure.metadata.end() && !traj->second.empty()) {
      std::ofstream ofs(output);
      if (!ofs.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output.string());
      }
      ofs << traj->second;
      frame_count = count_extxyz_frames(traj->second);
    } else {
      manager.write_as(loaded.structure, output, "xyz");
    }
  } else {
    manager.write_as(loaded.structure, output, to_format);
  }

  if (output_format == "json" || output_format == "jsonl") {
    const auto json = fmt::format(
        "{{\"schema_version\":\"1.0\",\"command\":\"convert\",\"input\":\"{}\",\"output\":\"{}\",\"detected_format\":\"{}\",\"target_format\":\"{}\",\"natoms\":{},\"frames\":{},\"status\":\"ok\"}}",
        json_escape(input.string()),
        json_escape(output.string()),
        json_escape(loaded.detected_format),
        json_escape(to_format),
        loaded.structure.natoms(),
        frame_count);
    if (output_format == "jsonl") {
      fmt::println("{}", json);
    } else {
      fmt::print("{}", json);
    }
  } else {
    if (to_format == "extxyz" && frame_count > 1) {
      fmt::println("converted {} -> {} (extxyz, {} frames)", input.string(), output.string(), frame_count);
    } else {
      fmt::println("converted {} -> {} ({})", input.string(), output.string(), to_format);
    }
  }
  return 0;
}
