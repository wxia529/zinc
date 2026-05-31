#include <CLI/CLI.hpp>
#include <fmt/core.h>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include "zinc/io/driver.h"

int run_scan(
    zinc::io::DriverManager& manager,
    const std::filesystem::path& path,
    const std::string& filter_element,
    const std::string& output_format);
int run_info(
    zinc::io::DriverManager& manager,
    const std::filesystem::path& file,
    const std::string& from_format,
    const std::string& output_format,
    const std::filesystem::path* output);
void run_update(
    zinc::io::DriverManager& manager,
    const std::string& source,
    const std::string& target,
    const std::filesystem::path* output,
    bool in_place,
    const std::string& pos_unit,
    const std::string& cell_unit,
    bool fill_pseudo,
    const std::string& report_format);
int run_convert(
    zinc::io::DriverManager& manager,
    const std::filesystem::path& input,
    const std::filesystem::path& output,
    const std::string& from_format,
    const std::string& to_format,
    const std::string& output_format);

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

static std::string normalize_format(std::string fmt) {
  for (auto& ch : fmt) {
    if (ch >= 'A' && ch <= 'Z') {
      ch = static_cast<char>(ch - 'A' + 'a');
    } else if (ch == '_') {
      ch = '-';
    }
  }
  if (fmt == "vasp") return "poscar";
  if (fmt == "extxyz" || fmt == "ext-xyz") return "extxyz";
  if (fmt == "qe-xml" || fmt == "qexml") return "qe-xml";
  if (fmt == "qe-log" || fmt == "qelog") return "qe-log";
  if (fmt == "qe-input" || fmt == "qeinput" || fmt == "pwi") return "qe-input";
  return fmt;
}

static std::optional<std::string> infer_format_from_path(const std::filesystem::path& path) {
  auto ext = path.extension().string();
  for (auto& ch : ext) {
    if (ch >= 'A' && ch <= 'Z') ch = static_cast<char>(ch - 'A' + 'a');
  }
  if (ext == ".cif") return "cif";
  if (ext == ".xyz") return "xyz";
  if (ext == ".extxyz") return "extxyz";
  if (ext == ".vasp" || ext == ".poscar" || ext == ".contcar") return "poscar";
  if (ext == ".xml") return "qe-xml";
  if (ext == ".out" || ext == ".log") return "qe-log";
  if (ext == ".in" || ext == ".pwi") return "qe-input";
  return std::nullopt;
}

static void validate_from_conflict(const std::filesystem::path& input, const std::string& from_format) {
  if (from_format.empty()) return;
  auto inferred = infer_format_from_path(input);
  if (inferred && normalize_format(from_format) != normalize_format(*inferred)) {
    throw std::invalid_argument(
        "format conflict: --from " + normalize_format(from_format) +
        " does not match input extension format " + *inferred);
  }
}

static void add_report_format_option(CLI::App* cmd, std::string& report_format) {
  cmd->add_option("--report-format,--format", report_format, "Command report format: text, json, jsonl")
      ->default_val("text")
      ->check(CLI::IsMember({"text", "json", "jsonl"}));
}

static std::string resolve_to_format(const std::filesystem::path& output, const std::string& to_format) {
  auto inferred = infer_format_from_path(output);
  if (!to_format.empty() && inferred && normalize_format(to_format) != normalize_format(*inferred)) {
    throw std::invalid_argument(
        "format conflict: --to " + normalize_format(to_format) +
        " does not match output extension format " + *inferred);
  }
  if (!to_format.empty()) return normalize_format(to_format);
  if (inferred) return *inferred;
  throw std::invalid_argument("missing output format: provide --to or an output filename with known extension");
}

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "inspect") {
    fmt::println(stderr, "usage error: inspect command has been removed; use info instead");
    return 64;
  }

  CLI::App app{"zinc CLI - crystal structure tooling"};

  auto* info_cmd = app.add_subcommand("info", "Show structure information");
  std::string info_file;
  std::string info_from;
  std::string info_output;
  std::string info_report_format = "text";
  info_cmd->add_option("file", info_file, "File path")->required();
  info_cmd->add_option("--from", info_from, "Input structure format");
  info_cmd->add_option("-o,--output", info_output, "Output report file");
  add_report_format_option(info_cmd, info_report_format);

  auto* convert_cmd = app.add_subcommand("convert", "Convert one structure file");
  std::string convert_input;
  std::string convert_output;
  std::string convert_from;
  std::string convert_to;
  std::string convert_report_format = "text";
  convert_cmd->add_option("file", convert_input, "Input file")->required();
  convert_cmd->add_option("-o,--output", convert_output, "Output structure file")->required();
  convert_cmd->add_option("--from", convert_from, "Input structure format");
  convert_cmd->add_option("--to", convert_to, "Output structure format");
  add_report_format_option(convert_cmd, convert_report_format);

  auto* scan_cmd = app.add_subcommand("scan", "Scan directory for structure files");
  std::string scan_path = ".";
  std::string filter_elem;
  std::string scan_report_format = "text";
  scan_cmd->add_option("path", scan_path, "Target directory")->default_val(".");
  scan_cmd->add_option("-e,--element", filter_elem, "Filter by element");
  add_report_format_option(scan_cmd, scan_report_format);

  auto* update_cmd = app.add_subcommand("update", "Update QE input with new structure");
  std::string src, tgt, update_output;
  bool in_place = false;
  std::string pos_unit = "angstrom", cell_unit = "angstrom";
  bool fill_pseudo = false;
  std::string update_report_format = "text";
  update_cmd->add_option("source", src, "Source structure file")->required();
  update_cmd->add_option("target", tgt, "Target QE input file")->required();
  update_cmd->add_option("-o,--output", update_output, "Output QE input file");
  update_cmd->add_flag("--in-place", in_place, "Overwrite target");
  update_cmd->add_option("--pos-unit", pos_unit, "Position unit")->default_val("angstrom");
  update_cmd->add_option("--cell-unit", cell_unit, "Cell unit")->default_val("angstrom");
  update_cmd->add_flag("--fill-pseudo", fill_pseudo, "Fill pseudopotentials");
  add_report_format_option(update_cmd, update_report_format);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    int code = app.exit(e);
    return (code == 0) ? 0 : 64;
  }

  zinc::io::DriverManager manager;
  std::string active_report_format = "text";
  try {
    if (*scan_cmd) {
      active_report_format = scan_report_format;
      return run_scan(manager, scan_path, filter_elem, scan_report_format);
    }

    if (*info_cmd) {
      active_report_format = info_report_format;
      validate_from_conflict(info_file, info_from);
      const std::filesystem::path out_path = info_output;
      const std::filesystem::path* out_ptr = info_output.empty() ? nullptr : &out_path;
      return run_info(manager, info_file, info_from, info_report_format, out_ptr);
    }

    if (*convert_cmd) {
      active_report_format = convert_report_format;
      validate_from_conflict(convert_input, convert_from);
      const auto to_resolved = resolve_to_format(convert_output, convert_to);
      return run_convert(manager, convert_input, convert_output, convert_from, to_resolved, convert_report_format);
    }

    if (*update_cmd) {
      active_report_format = update_report_format;
      const std::filesystem::path up_fs_path = update_output;
      const std::filesystem::path* out_ptr = update_output.empty() ? nullptr : &up_fs_path;
      run_update(manager, src, tgt, out_ptr, in_place, pos_unit, cell_unit, fill_pseudo, update_report_format);
      return 0;
    }

    fmt::print("{}", app.help());
    return 0;
  } catch (const std::invalid_argument& e) {
    fmt::println(stderr, "usage error: {}", e.what());
    return 64;
  } catch (const std::exception& e) {
    if (active_report_format == "json" || active_report_format == "jsonl") {
      const auto err = fmt::format(
          "{{\"error\":{{\"code\":\"EXECUTION_ERROR\",\"message\":\"{}\"}}}}",
          json_escape(e.what()));
      fmt::println(stderr, "{}", err);
    } else {
      fmt::println(stderr, "error: {}", e.what());
    }
    return 1;
  }
}
