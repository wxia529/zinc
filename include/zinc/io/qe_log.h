#pragma once

#include "zinc/io/driver.h"

namespace zinc {
namespace io {

class QELogDriver : public StructureDriver {
public:
  const char* name() const override { return "qe_log"; }
  bool identify(const std::string& first_line) const override;
  core::Structure read(const std::filesystem::path& path) const override;
  void write(const core::Structure& structure, const std::filesystem::path& path) const override;
};

} // namespace io
} // namespace zinc
