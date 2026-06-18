#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <sstream>
#include <array>
#include "zinc/io/driver.h"
#include "zinc/io/poscar.h"
#include "zinc/io/xyz.h"
#include "zinc/io/cif.h"

using namespace zinc::io;

static std::filesystem::path find_example_file(const std::string& name) {
  std::filesystem::path p = std::filesystem::current_path();
  for (int i = 0; i < 8; ++i) {
    auto candidate = p / "Example" / name;
    if (std::filesystem::exists(candidate)) return candidate;
    if (!p.has_parent_path()) break;
    p = p.parent_path();
  }
  throw std::runtime_error("Cannot locate Example/" + name);
}

static std::string read_file(const std::filesystem::path& path) {
  std::ifstream ifs(path);
  std::ostringstream oss;
  oss << ifs.rdbuf();
  return oss.str();
}

TEST(IOTest, PoscarRoundTrip) {
  std::filesystem::path path = "/tmp/test_poscar.vasp";
  zinc::core::Lattice lat(zinc::core::Vector3d(3.0,0,0), zinc::core::Vector3d(0,3.0,0), zinc::core::Vector3d(0,0,3.0));
  zinc::core::Structure s(lat);
  s.add_atom("Si", {0.0, 0.0, 0.0});
  s.add_atom("Si", {1.5, 1.5, 1.5});

  PoscarDriver driver;
  driver.write(s, path);
  auto s2 = driver.read(path);

  EXPECT_EQ(s2.natoms(), 2);
  EXPECT_DOUBLE_EQ(s2.lattice.volume(), 27.0);
  std::filesystem::remove(path);
}

TEST(IOTest, PoscarWriteUsesFixedTenDecimalPlaces) {
  std::filesystem::path path = "/tmp/test_poscar_format.vasp";
  zinc::core::Lattice lat(zinc::core::Vector3d(3.0,0,0), zinc::core::Vector3d(0,3.0,0), zinc::core::Vector3d(0,0,3.0));
  zinc::core::Structure s(lat);
  s.add_atom("Si", {0.125, -1.5, 2.0});

  PoscarDriver().write(s, path);
  auto text = read_file(path);

  EXPECT_NE(text.find("      3.0000000000"), std::string::npos);
  EXPECT_NE(text.find("      0.1250000000"), std::string::npos);
  EXPECT_NE(text.find("     -1.5000000000"), std::string::npos);
  std::filesystem::remove(path);
}

TEST(IOTest, PoscarOmitsSelectiveDynamicsWhenAllAtomsFree) {
  auto path = std::filesystem::temp_directory_path() / "zinc-poscar-free.vasp";
  zinc::core::Lattice lat(zinc::core::Vector3d(3.0,0,0), zinc::core::Vector3d(0,3.0,0), zinc::core::Vector3d(0,0,3.0));
  zinc::core::Structure s(lat);
  s.add_atom("Si", {0.0, 0.0, 0.0}, std::array<int, 3>{1, 1, 1});
  s.add_atom("Si", {1.5, 1.5, 1.5});

  PoscarDriver().write(s, path);
  auto text = read_file(path);

  EXPECT_EQ(text.find("Selective dynamics"), std::string::npos);
  EXPECT_EQ(text.find("  T  T  T"), std::string::npos);
  std::filesystem::remove(path);
}

TEST(IOTest, PoscarWritesSelectiveDynamicsWhenAnyAtomFixed) {
  auto path = std::filesystem::temp_directory_path() / "zinc-poscar-fixed.vasp";
  zinc::core::Lattice lat(zinc::core::Vector3d(3.0,0,0), zinc::core::Vector3d(0,3.0,0), zinc::core::Vector3d(0,0,3.0));
  zinc::core::Structure s(lat);
  s.add_atom("Si", {0.0, 0.0, 0.0}, std::array<int, 3>{0, 1, 1});
  s.add_atom("Si", {1.5, 1.5, 1.5});

  PoscarDriver().write(s, path);
  auto text = read_file(path);

  EXPECT_NE(text.find("Selective dynamics"), std::string::npos);
  EXPECT_NE(text.find("  F  T  T"), std::string::npos);
  EXPECT_NE(text.find("  T  T  T"), std::string::npos);
  std::filesystem::remove(path);
}

TEST(IOTest, XyzRoundTrip) {
  std::filesystem::path path = "/tmp/test_xyz.xyz";
  zinc::core::Lattice lat(zinc::core::Vector3d(5.0,0,0), zinc::core::Vector3d(0,5.0,0), zinc::core::Vector3d(0,0,5.0));
  zinc::core::Structure s(lat);
  s.add_atom("C", {0.0, 0.0, 0.0});

  XyzDriver driver;
  driver.write(s, path);
  auto s2 = driver.read(path);

  EXPECT_EQ(s2.natoms(), 1);
  EXPECT_DOUBLE_EQ(s2.lattice.volume(), 125.0);
  std::filesystem::remove(path);
}

TEST(IOTest, XyzWriteUsesFixedTenDecimalPlaces) {
  std::filesystem::path path = "/tmp/test_xyz_format.xyz";
  zinc::core::Lattice lat(zinc::core::Vector3d(5.0,0,0), zinc::core::Vector3d(0,5.0,0), zinc::core::Vector3d(0,0,5.0));
  zinc::core::Structure s(lat);
  s.add_atom("C", {0.25, 1.0, -2.5});

  XyzDriver().write(s, path);
  auto text = read_file(path);

  EXPECT_NE(text.find("      5.0000000000"), std::string::npos);
  EXPECT_NE(text.find("C         0.2500000000"), std::string::npos);
  EXPECT_NE(text.find("     -2.5000000000"), std::string::npos);
  std::filesystem::remove(path);
}

TEST(IOTest, XyzRejectsMalformedLattice) {
  auto path = std::filesystem::temp_directory_path() / "zinc-bad-lattice.xyz";
  {
    std::ofstream ofs(path);
    ofs << "1\nLattice=\"1 0\"\nH 0 0 0\n";
  }

  EXPECT_THROW(XyzDriver().read(path), std::runtime_error);
  std::filesystem::remove(path);
}

TEST(IOTest, XyzRejectsMalformedAtom) {
  auto path = std::filesystem::temp_directory_path() / "zinc-bad-atom.xyz";
  {
    std::ofstream ofs(path);
    ofs << "1\ncomment\nH bad 0 0\n";
  }

  EXPECT_THROW(XyzDriver().read(path), std::runtime_error);
  std::filesystem::remove(path);
}

TEST(IOTest, PoscarToXyzKeepsLattice) {
  auto path = std::filesystem::temp_directory_path() / "zinc-yvo4.xyz";
  DriverManager mgr;
  auto s = mgr.load_auto(find_example_file("yvo4.vasp"));

  mgr.write_as(s, path, "xyz");
  auto s2 = mgr.load_auto(path);
  const auto& m = s2.lattice.matrix();

  EXPECT_NEAR(m(0, 0), 3.5933092523, 1e-10);
  EXPECT_NEAR(m(0, 1), 3.5933092523, 1e-10);
  EXPECT_NEAR(m(0, 2), 3.1556474122, 1e-10);
  EXPECT_NEAR(m(1, 0), -3.5933092523, 1e-10);
  EXPECT_NEAR(m(2, 2), -3.1556474122, 1e-10);
  std::filesystem::remove(path);
}

TEST(IOTest, DriverManagerAutoDetect) {
  std::filesystem::path path = "/tmp/test_auto.xyz";
  zinc::core::Structure s;
  s.add_atom("H", {0.0, 0.0, 0.0});

  XyzDriver().write(s, path);

  DriverManager mgr;
  auto s2 = mgr.load_auto(path);
  EXPECT_EQ(s2.natoms(), 1);
  std::filesystem::remove(path);
}

TEST(IOTest, QELogParsesLatticeAndEnergy) {
  DriverManager mgr;
  auto s = mgr.load_auto(find_example_file("pwscf.out"));

  EXPECT_EQ(s.natoms(), 12);
  EXPECT_NEAR(s.lattice.volume(), 162.981254024, 1e-6);

  auto it = s.metadata.find("qe_total_energy_ry");
  ASSERT_NE(it, s.metadata.end());
  EXPECT_NEAR(std::stod(it->second), -800.0241048, 1e-9);
}

TEST(IOTest, QEInputParsesRelaxInput) {
  DriverManager mgr;
  auto s = mgr.load_auto(find_example_file("pwscf.in"));

  EXPECT_EQ(s.natoms(), 12);
  EXPECT_NEAR(s.lattice.volume(), 160.263235223, 1e-6);
  ASSERT_FALSE(s.atoms.empty());
  EXPECT_EQ(s.atoms[0].element, "O");
  EXPECT_NEAR(s.atoms[0].position.x(), 0.0, 1e-9);
  EXPECT_NEAR(s.atoms[0].position.y(), -2.250239923, 1e-9);
}

TEST(IOTest, QEInputParsesCrystalCoordinates) {
  auto path = std::filesystem::temp_directory_path() / "zinc-qe-crystal.pwi";
  {
    std::ofstream ofs(path);
    ofs << "&SYSTEM\n  ibrav = 0\n  nat = 1\n  ntyp = 1\n/\n"
        << "ATOMIC_SPECIES\nSi 28.085 Si.UPF\n"
        << "CELL_PARAMETERS angstrom\n"
        << "2.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 6.0\n"
        << "ATOMIC_POSITIONS crystal\n"
        << "Si 0.5 0.25 0.75\n";
  }

  DriverManager mgr;
  auto s = mgr.load_auto(path);

  EXPECT_EQ(s.natoms(), 1);
  EXPECT_NEAR(s.atoms[0].position.x(), 1.0, 1e-12);
  EXPECT_NEAR(s.atoms[0].position.y(), 1.0, 1e-12);
  EXPECT_NEAR(s.atoms[0].position.z(), 4.5, 1e-12);
  std::filesystem::remove(path);
}

TEST(IOTest, QEXmlParsesNestedAtomicStructure) {
  DriverManager mgr;
  auto s = mgr.load_auto(find_example_file("data-file-schema.xml"));

  EXPECT_EQ(s.natoms(), 12);
  EXPECT_NEAR(s.lattice.volume(), 162.981254064, 1e-6);
  ASSERT_FALSE(s.atoms.empty());
  EXPECT_NEAR(s.atoms[0].position.x(), 0.0, 1e-9);
  EXPECT_NEAR(s.atoms[0].position.y(), -2.270535575, 1e-9);
  auto it = s.metadata.find("frames");
  ASSERT_NE(it, s.metadata.end());
  EXPECT_GT(std::stoi(it->second), 1);
}

TEST(IOTest, PoscarKeepsLatticeVectorOrder) {
  DriverManager mgr;
  auto s = mgr.load_auto(find_example_file("yvo4.vasp"));
  const auto& m = s.lattice.matrix();

  EXPECT_NEAR(m(0, 0), 3.5933092523, 1e-10);
  EXPECT_NEAR(m(0, 1), 3.5933092523, 1e-10);
  EXPECT_NEAR(m(0, 2), 3.1556474122, 1e-10);

  EXPECT_NEAR(m(1, 0), -3.5933092523, 1e-10);
  EXPECT_NEAR(m(1, 1), -3.5933092523, 1e-10);
  EXPECT_NEAR(m(1, 2), 3.1556474122, 1e-10);

  EXPECT_NEAR(m(2, 0), 3.5933092523, 1e-10);
  EXPECT_NEAR(m(2, 1), -3.5933092523, 1e-10);
  EXPECT_NEAR(m(2, 2), -3.1556474122, 1e-10);
}
