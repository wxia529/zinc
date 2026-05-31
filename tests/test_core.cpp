#include <gtest/gtest.h>
#include "zinc/core/constants.h"
#include "zinc/core/lattice.h"
#include "zinc/core/structure.h"

using namespace zinc::core;

TEST(ConstantsTest, BohrToAngstrom) {
  EXPECT_DOUBLE_EQ(to_angstrom(LengthUnit::Bohr, 1.0), BOHR_TO_ANGSTROM);
  EXPECT_DOUBLE_EQ(to_angstrom(LengthUnit::Angstrom, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(to_angstrom(LengthUnit::Alat, 1.0, 10.0), 10.0 * BOHR_TO_ANGSTROM);
}

TEST(LatticeTest, Volume) {
  Lattice lat(Vector3d(2,0,0), Vector3d(0,3,0), Vector3d(0,0,4));
  EXPECT_DOUBLE_EQ(lat.volume(), 24.0);
}

TEST(LatticeTest, FracCartRoundTrip) {
  Lattice lat(Vector3d(3,0,0), Vector3d(0,3,0), Vector3d(0,0,3));
  Vector3d frac(0.5, 0.25, 0.75);
  Vector3d cart = lat.frac_to_cart(frac);
  auto back = lat.cart_to_frac(cart);
  ASSERT_TRUE(back.has_value());
  EXPECT_NEAR(back->x(), frac.x(), 1e-9);
  EXPECT_NEAR(back->y(), frac.y(), 1e-9);
  EXPECT_NEAR(back->z(), frac.z(), 1e-9);
}

TEST(StructureTest, AddAtom) {
  Structure s;
  s.add_atom("C", {0.0, 0.0, 0.0});
  s.add_atom("H", {1.0, 0.0, 0.0});
  EXPECT_EQ(s.natoms(), 2);
}
