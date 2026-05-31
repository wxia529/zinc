#include "zinc/io/xyz.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>

namespace zinc {
namespace io {

bool XyzDriver::identify(const std::string& first_line) const {
  std::istringstream ss(first_line);
  int nat;
  return (ss >> nat) && nat > 0;
}

core::Structure XyzDriver::read(const std::filesystem::path& path) const {
  std::ifstream file(path);
  if (!file.is_open()) throw std::runtime_error("Cannot open: " + path.string());

  std::string line;
  std::getline(file, line);
  int nat = std::stoi(line);
  if (nat <= 0) throw std::runtime_error("Invalid XYZ atom count");
  std::getline(file, line);

  core::Lattice lattice;
  bool lattice_provided = false;
  if (line.find("Lattice=") != std::string::npos) {
    size_t s = line.find('"'), e = line.find('"', s+1);
    if (s == std::string::npos || e == std::string::npos) {
      throw std::runtime_error("Invalid XYZ Lattice field");
    }
    std::istringstream ss(line.substr(s+1, e-s-1));
    double v[9];
    for (double& value : v) {
      if (!(ss >> value)) throw std::runtime_error("Invalid XYZ Lattice field");
    }
    lattice = core::Lattice(core::Vector3d(v[0],v[1],v[2]), core::Vector3d(v[3],v[4],v[5]), core::Vector3d(v[6],v[7],v[8]));
    lattice_provided = true;
  }

  core::Structure structure(lattice);
  structure.metadata["lattice_missing"] = lattice_provided ? "false" : "true";
  for (int i = 0; i < nat; ++i) {
    if (!std::getline(file, line)) throw std::runtime_error("Invalid XYZ: missing atom line");
    std::istringstream ss(line);
    std::string elem; double x, y, z, fx=0, fy=0, fz=0;
    if (!(ss >> elem >> x >> y >> z)) throw std::runtime_error("Invalid XYZ atom line");
    bool has_f = false;
    if (ss >> fx) {
      if (!(ss >> fy >> fz)) throw std::runtime_error("Invalid XYZ force columns");
      has_f = true;
    }
    if (has_f) {
      core::Atom atom; atom.element = elem; atom.position = core::Vector3d(x,y,z); atom.force = core::Vector3d(fx,fy,fz);
      structure.atoms.push_back(std::move(atom));
    } else {
      structure.add_atom(elem, {x, y, z});
    }
  }
  return structure;
}

void XyzDriver::write(const core::Structure& structure, const std::filesystem::path& path) const {
  std::ofstream out(path);
  if (!out.is_open()) throw std::runtime_error("Cannot open output file: " + path.string());
  out << std::fixed << std::setprecision(10);
  const auto& mat = structure.lattice.matrix();
  out << structure.natoms() << "\nLattice=\""
      << std::setw(18) << mat(0,0) << " " << std::setw(18) << mat(0,1) << " " << std::setw(18) << mat(0,2) << " "
      << std::setw(18) << mat(1,0) << " " << std::setw(18) << mat(1,1) << " " << std::setw(18) << mat(1,2) << " "
      << std::setw(18) << mat(2,0) << " " << std::setw(18) << mat(2,1) << " " << std::setw(18) << mat(2,2) << "\" ";
  bool hf = std::any_of(structure.atoms.begin(), structure.atoms.end(),
      [](const auto& a) { return a.force.has_value(); });
  out << "Properties=species:S:1:pos:R:3" << (hf ? ":forces:R:3" : "") << "\n";
  for (const auto& atom : structure.atoms) {
    out << std::left << std::setw(4) << atom.element << std::right
        << std::setw(18) << atom.position.x()
        << std::setw(18) << atom.position.y()
        << std::setw(18) << atom.position.z();
    if (hf && atom.force) {
      out << std::setw(18) << atom.force->x()
          << std::setw(18) << atom.force->y()
          << std::setw(18) << atom.force->z();
    }
    out << "\n";
  }
}

} // namespace io
} // namespace zinc
