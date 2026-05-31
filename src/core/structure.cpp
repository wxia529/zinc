#include "zinc/core/structure.h"

namespace zinc {
namespace core {

Structure::Structure(const Lattice& lat) : lattice(lat) {}

void Structure::add_atom(const std::string& element, const std::array<double, 3>& position) {
  Atom atom;
  atom.element = element;
  atom.position = Vector3d(position[0], position[1], position[2]);
  atoms.push_back(std::move(atom));
}

void Structure::add_atom(const std::string& element, const std::array<double, 3>& position,
                         const std::optional<std::array<int, 3>>& constraints) {
  Atom atom;
  atom.element = element;
  atom.position = Vector3d(position[0], position[1], position[2]);
  atom.constraints = constraints;
  atoms.push_back(std::move(atom));
}

} // namespace core
} // namespace zinc
