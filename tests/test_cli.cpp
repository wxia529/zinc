#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

namespace {

namespace fs = std::filesystem;

struct CommandResult {
  int exit_code;
  std::string stdout_text;
  std::string stderr_text;
};

static fs::path find_cli() {
  fs::path p = fs::current_path();
  for (int i = 0; i < 8; ++i) {
    const std::vector<fs::path> candidates = {
      p / "zinc",
      p / "build" / "zinc",
      p / "build" / "Debug" / "zinc",
      p / "build" / "Release" / "zinc",
    };
    for (const auto& candidate : candidates) {
      if (fs::exists(candidate)) {
        return candidate;
      }
    }
    if (!p.has_parent_path()) {
      break;
    }
    p = p.parent_path();
  }
  throw std::runtime_error("Cannot locate zinc executable");
}

static fs::path find_repo_root() {
  fs::path p = fs::current_path();
  for (int i = 0; i < 8; ++i) {
    if (fs::exists(p / "Example")) {
      return p;
    }
    if (!p.has_parent_path()) {
      break;
    }
    p = p.parent_path();
  }
  throw std::runtime_error("Cannot locate repository root");
}

static std::string read_file(const fs::path& path) {
  std::ifstream ifs(path);
  std::ostringstream oss;
  oss << ifs.rdbuf();
  return oss.str();
}

static CommandResult run_cli(const std::string& args) {
  const auto stamp = std::to_string(std::rand());
  const auto out_path = fs::temp_directory_path() / ("zinc-test-out-" + stamp + ".txt");
  const auto err_path = fs::temp_directory_path() / ("zinc-test-err-" + stamp + ".txt");
  const auto exe = find_cli();
  std::string cmd = "\"" + exe.string() + "\" " + args + " >\"" + out_path.string() + "\" 2>\"" + err_path.string() + "\"";
  int rc = std::system(cmd.c_str());
  int code = (rc == -1) ? 1 : (rc >> 8);
  CommandResult r{code, read_file(out_path), read_file(err_path)};
  fs::remove(out_path);
  fs::remove(err_path);
  return r;
}

static fs::path find_example_file(const std::string& name) {
  auto candidate = find_repo_root() / "Example" / name;
  if (fs::exists(candidate)) {
    return candidate;
  }
  throw std::runtime_error("Cannot locate Example/" + name);
}

}  // namespace

TEST(CLITest, HelpListsCanonicalCommands) {
  auto r = run_cli("--help");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("info"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("convert"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("scan"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("update"), std::string::npos);
  EXPECT_EQ(r.stdout_text.find("inspect"), std::string::npos);
}

TEST(CLITest, InvalidUsageReturns64) {
  auto r = run_cli("info");
  EXPECT_EQ(r.exit_code, 64);
  EXPECT_FALSE(r.stderr_text.empty());
}

TEST(CLITest, InfoJsonContainsSchemaVersion) {
  auto file = find_example_file("yvo4.cif");
  auto r = run_cli("info \"" + file.string() + "\" --report-format json");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("\"schema_version\""), std::string::npos);
  EXPECT_NE(r.stdout_text.find("\"command\""), std::string::npos);
  EXPECT_NE(r.stdout_text.find("\"detected_format\""), std::string::npos);
}

TEST(CLITest, FormatAliasStillWorksForReports) {
  auto file = find_example_file("yvo4.cif");
  auto r = run_cli("info \"" + file.string() + "\" --format json");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("\"command\":\"info\""), std::string::npos);
}

TEST(CLITest, ScanReturnsPartialSuccessWhenSomeFilesFail) {
  auto dir = fs::temp_directory_path() / ("zinc-scan-" + std::to_string(std::rand()));
  fs::create_directories(dir);
  auto ok = dir / "ok.xyz";
  auto bad = dir / "bad.xyz";
  {
    std::ofstream ofs(ok);
    ofs << "1\ncomment\nH 0 0 0\n";
  }
  {
    std::ofstream ofs(bad);
    ofs << "not-a-valid-structure";
  }
  auto r = run_cli("scan \"" + dir.string() + "\" --format jsonl");
  EXPECT_EQ(r.exit_code, 2);
  EXPECT_NE(r.stdout_text.find("\"status\":\"ok\""), std::string::npos);
  EXPECT_NE(r.stderr_text.find("parse failed"), std::string::npos);
  fs::remove_all(dir);
}

TEST(CLITest, ScanIgnoresUnrelatedFiles) {
  auto dir = fs::temp_directory_path() / ("zinc-scan-" + std::to_string(std::rand()));
  fs::create_directories(dir);
  {
    std::ofstream ofs(dir / ".gitignore");
    ofs << "build\n";
  }
  {
    std::ofstream ofs(dir / "notes.txt");
    ofs << "not a structure\n";
  }
  {
    std::ofstream ofs(dir / "ok.xyz");
    ofs << "1\ncomment\nH 0 0 0\n";
  }

  auto r = run_cli("scan \"" + dir.string() + "\" --format jsonl");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("\"input\":\"ok.xyz\""), std::string::npos);
  EXPECT_EQ(r.stdout_text.find(".gitignore"), std::string::npos);
  EXPECT_EQ(r.stdout_text.find("notes.txt"), std::string::npos);
  EXPECT_TRUE(r.stderr_text.empty());
  fs::remove_all(dir);
}

TEST(CLITest, ConvertRejectsToAndOutputExtensionConflict) {
  auto file = find_example_file("yvo4.cif");
  auto out = fs::temp_directory_path() / ("zinc-conflict-" + std::to_string(std::rand()) + ".xyz");
  auto r = run_cli("convert \"" + file.string() + "\" -o \"" + out.string() + "\" --to poscar");
  EXPECT_EQ(r.exit_code, 64);
  EXPECT_NE(r.stderr_text.find("format conflict"), std::string::npos);
}

TEST(CLITest, InfoRejectsFromAndInputExtensionConflict) {
  auto file = find_example_file("yvo4.cif");
  auto r = run_cli("info \"" + file.string() + "\" --from xyz");
  EXPECT_EQ(r.exit_code, 64);
  EXPECT_NE(r.stderr_text.find("format conflict"), std::string::npos);
}

TEST(CLITest, InfoAcceptsQEInputExtension) {
  auto file = find_example_file("pwscf.in");
  auto r = run_cli("info \"" + file.string() + "\" --report-format json");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("\"detected_format\":\"qe-input\""), std::string::npos);
  EXPECT_NE(r.stdout_text.find("\"natoms\":12"), std::string::npos);
}

TEST(CLITest, InspectCommandRemoved) {
  auto r = run_cli("inspect --help");
  EXPECT_EQ(r.exit_code, 64);
}

TEST(CLITest, InfoHelpShowsOnlyInfoOptions) {
  auto r = run_cli("info --help");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("--from"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("--report-format"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("--format"), std::string::npos);
  EXPECT_EQ(r.stdout_text.find("--to"), std::string::npos);
}

TEST(CLITest, ConvertHelpSeparatesStructureAndReportFormats) {
  auto r = run_cli("convert --help");
  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("Input structure format"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("Output structure format"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("Command report format"), std::string::npos);
}

TEST(CLITest, ConvertQEXmlOptimizationToExtxyzWritesTrajectory) {
  auto source = find_example_file("data-file-schema.xml");
  auto out = fs::temp_directory_path() / ("zinc-trajectory-" + std::to_string(std::rand()) + ".extxyz");
  auto r = run_cli("convert \"" + source.string() + "\" -o \"" + out.string() +
                   "\" --report-format json");

  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("\"target_format\":\"extxyz\""), std::string::npos);
  EXPECT_NE(r.stdout_text.find("\"frames\":"), std::string::npos);

  auto text = read_file(out);
  EXPECT_EQ(text.rfind("12\nLattice=", 0), 0);
  auto second = text.find("\n12\nLattice=", 1);
  EXPECT_NE(second, std::string::npos);
  EXPECT_NE(text.find("Properties=species:S:1:pos:R:3:forces:R:3"), std::string::npos);
  fs::remove(out);
}

TEST(CLITest, InfoJsonOmitsInternalTrajectoryText) {
  auto source = find_example_file("data-file-schema.xml");
  auto r = run_cli("info \"" + source.string() + "\" --report-format json");

  EXPECT_EQ(r.exit_code, 0);
  EXPECT_EQ(r.stdout_text.find("trajectory_extxyz"), std::string::npos);
  EXPECT_NE(r.stdout_text.find("\"frames\":\"11\""), std::string::npos);
}

TEST(CLITest, ConvertRejectsExtxyzExtensionConflict) {
  auto source = find_example_file("data-file-schema.xml");
  auto out = fs::temp_directory_path() / ("zinc-conflict-" + std::to_string(std::rand()) + ".extxyz");
  auto r = run_cli("convert \"" + source.string() + "\" -o \"" + out.string() + "\" --to xyz");
  EXPECT_EQ(r.exit_code, 64);
  EXPECT_NE(r.stderr_text.find("format conflict"), std::string::npos);
}

TEST(CLITest, ConvertReportsOutputOpenFailure) {
  auto source = find_example_file("yvo4.vasp");
  auto dir = fs::temp_directory_path() / ("zinc-missing-" + std::to_string(std::rand()));
  auto out = dir / "out.xyz";
  auto r = run_cli("convert \"" + source.string() + "\" -o \"" + out.string() + "\"");

  EXPECT_EQ(r.exit_code, 1);
  EXPECT_NE(r.stderr_text.find("Cannot open output file"), std::string::npos);
}

TEST(CLITest, InfoReportsOutputOpenFailure) {
  auto source = find_example_file("yvo4.vasp");
  auto dir = fs::temp_directory_path() / ("zinc-missing-" + std::to_string(std::rand()));
  auto out = dir / "report.txt";
  auto r = run_cli("info \"" + source.string() + "\" -o \"" + out.string() + "\"");

  EXPECT_EQ(r.exit_code, 1);
  EXPECT_NE(r.stderr_text.find("Cannot open output file"), std::string::npos);
}

TEST(CLITest, UpdateRejectsUnsupportedUnits) {
  auto source = find_example_file("yvo4.vasp");
  auto target = find_example_file("pwscf.in");
  auto r = run_cli("update \"" + source.string() + "\" \"" + target.string() +
                   "\" --pos-unit nonsense");

  EXPECT_EQ(r.exit_code, 64);
  EXPECT_NE(r.stderr_text.find("unsupported --pos-unit"), std::string::npos);
}

TEST(CLITest, UpdateAcceptsQEXmlDataFileSchema) {
  auto source = find_example_file("data-file-schema.xml");
  auto target = find_example_file("pwscf.in");
  auto out = fs::temp_directory_path() / ("zinc-updated-" + std::to_string(std::rand()) + ".in");
  auto r = run_cli("update \"" + source.string() + "\" \"" + target.string() +
                   "\" -o \"" + out.string() + "\" --report-format json");

  EXPECT_EQ(r.exit_code, 0);
  EXPECT_NE(r.stdout_text.find("\"status\":\"ok\""), std::string::npos);
  auto text = read_file(out);
  EXPECT_NE(text.find("ATOMIC_POSITIONS (angstrom)"), std::string::npos);
  EXPECT_NE(text.find("     -2.2705355746"), std::string::npos);
  EXPECT_NE(text.find("      4.2529271402"), std::string::npos);
  fs::remove(out);
}
