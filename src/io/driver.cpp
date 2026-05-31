#include "zinc/io/driver.h"
#include "zinc/io/qe_xml.h"
#include "zinc/io/qe_log.h"
#include "zinc/io/qe_input.h"
#include "zinc/io/poscar.h"
#include "zinc/io/cif.h"
#include "zinc/io/xyz.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace zinc {
namespace io {

static std::string normalize_format(std::string fmt) {
  for (auto& ch : fmt) {
    if (ch >= 'A' && ch <= 'Z') {
      ch = static_cast<char>(ch - 'A' + 'a');
    } else if (ch == '_') {
      ch = '-';
    }
  }
  if (fmt == "vasp") {
    return "poscar";
  }
  if (fmt == "qe-xml" || fmt == "qexml") {
    return "qe-xml";
  }
  if (fmt == "qe-log" || fmt == "qelog") {
    return "qe-log";
  }
  if (fmt == "qe-input" || fmt == "qeinput" || fmt == "pwi") {
    return "qe-input";
  }
  return fmt;
}

static std::string normalize_driver_name(std::string name) {
  name = normalize_format(std::move(name));
  if (name == "qe-xml") {
    return "qe-xml";
  }
  if (name == "qe-log") {
    return "qe-log";
  }
  if (name == "qe-input") {
    return "qe-input";
  }
  return name;
}

DriverManager::DriverManager() {
  register_driver(std::make_unique<QELogDriver>());
  register_driver(std::make_unique<QEXmlDriver>());
  register_driver(std::make_unique<QEInputDriver>());
  register_driver(std::make_unique<CifDriver>());
  register_driver(std::make_unique<XyzDriver>());
  register_driver(std::make_unique<PoscarDriver>());
}

void DriverManager::register_driver(std::unique_ptr<StructureDriver> driver) {
  drivers_.push_back(std::move(driver));
}

core::Structure DriverManager::load_as(const std::filesystem::path& path, const std::string& format) const {
  const std::string want = normalize_format(format);
  for (const auto& driver : drivers_) {
    if (normalize_driver_name(driver->name()) == want) {
      return driver->read(path);
    }
  }
  throw std::runtime_error("Unsupported input format: " + format);
}

core::Structure DriverManager::load_auto(const std::filesystem::path& path) const {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }

  std::ifstream file(path, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + path.string());
  }

  char buffer[1024] = {};
  file.read(buffer, sizeof(buffer));
  std::streamsize n = file.gcount();
  file.close();
  std::string header(buffer, static_cast<size_t>(n));

  for (const auto& driver : drivers_) {
    if (driver->identify(header)) {
      return driver->read(path);
    }
  }

  throw std::runtime_error("Unsupported file format: " + path.string());
}

DriverManager::LoadResult DriverManager::load_with_detection(
    const std::filesystem::path& path,
    const std::string& from_format) const {
  if (!from_format.empty()) {
    return {load_as(path, from_format), normalize_format(from_format)};
  }

  if (!std::filesystem::exists(path)) {
    throw std::runtime_error("File not found: " + path.string());
  }

  std::ifstream file(path, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + path.string());
  }

  char buffer[1024] = {};
  file.read(buffer, sizeof(buffer));
  std::streamsize n = file.gcount();
  file.close();
  std::string header(buffer, static_cast<size_t>(n));

  for (const auto& driver : drivers_) {
    if (driver->identify(header)) {
      return {driver->read(path), normalize_driver_name(driver->name())};
    }
  }
  throw std::runtime_error("Unsupported file format: " + path.string());
}

void DriverManager::write_as(const core::Structure& structure,
                             const std::filesystem::path& path,
                             const std::string& format) const {
  const std::string want = normalize_format(format);
  for (const auto& driver : drivers_) {
    if (normalize_driver_name(driver->name()) == want) {
      driver->write(structure, path);
      return;
    }
  }
  throw std::runtime_error("Unsupported output format: " + format);
}

} // namespace io
} // namespace zinc
