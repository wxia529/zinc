#pragma once

#include "zinc/core/structure.h"
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace zinc {
namespace io {

enum class IoError {
  NotFound,
  ParseError,
  MissingField,
  Unsupported
};

class StructureDriver {
public:
  virtual ~StructureDriver() = default;
  virtual const char* name() const = 0;
  virtual bool identify(const std::string& first_line) const = 0;
  virtual core::Structure read(const std::filesystem::path& path) const = 0;
  virtual void write(const core::Structure& structure, const std::filesystem::path& path) const = 0;
};

class DriverManager {
public:
  struct LoadResult {
    core::Structure structure;
    std::string detected_format;
  };

  DriverManager();
  void register_driver(std::unique_ptr<StructureDriver> driver);
  core::Structure load_as(const std::filesystem::path& path, const std::string& format) const;
  core::Structure load_auto(const std::filesystem::path& path) const;
  LoadResult load_with_detection(const std::filesystem::path& path, const std::string& from_format = "") const;
  void write_as(const core::Structure& structure, const std::filesystem::path& path, const std::string& format) const;

private:
  std::vector<std::unique_ptr<StructureDriver>> drivers_;
};

} // namespace io
} // namespace zinc
