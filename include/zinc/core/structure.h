#pragma once

#include "zinc/core/lattice.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <optional>
#include <array>

namespace zinc {
namespace core {

struct Atom {
  std::string element;
  Vector3d position;
  std::optional<Vector3d> force;
  std::optional<std::array<int, 3>> constraints;
  std::unordered_map<std::string, std::string> properties;
};

struct Structure {
  Lattice lattice;
  std::vector<Atom> atoms;
  std::unordered_map<std::string, std::string> metadata;

  Structure() = default;
  explicit Structure(const Lattice& lat);

  void add_atom(const std::string& element, const std::array<double, 3>& position);
  void add_atom(const std::string& element, const std::array<double, 3>& position,
                const std::optional<std::array<int, 3>>& constraints);

  size_t natoms() const { return atoms.size(); }
};

} // namespace core
} // namespace zinc
