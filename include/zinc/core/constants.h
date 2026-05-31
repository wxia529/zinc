#pragma once

#include <cmath>
#include <string>
#include <variant>

namespace zinc {
namespace core {

inline constexpr double PI = 3.14159265358979323846;
inline constexpr double BOHR_TO_ANGSTROM = 0.529177210903;
inline constexpr double BOHR_TO_NM = BOHR_TO_ANGSTROM / 10.0;
inline constexpr double ANGSTROM_TO_NM = 0.1;
inline constexpr double RY_TO_EV = 13.605693122994;
inline constexpr double HARTREE_TO_EV = 27.211386245988;
inline constexpr double EV_TO_KCAL_MOL = 23.0605419454;
inline constexpr double EV_TO_CM1 = 8065.5440050769;
inline constexpr double FORCE_CONSTANT_RY_BOHR_TO_EV_ANGSTROM = RY_TO_EV / BOHR_TO_ANGSTROM;
inline constexpr double FORCE_CONSTANT_HARTREE_BOHR_TO_EV_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM;

enum class LengthUnit {
  Angstrom,
  Bohr,
  Alat
};

inline double to_angstrom(LengthUnit unit, double value, double alat_scale = 1.0) {
  switch (unit) {
    case LengthUnit::Angstrom: return value;
    case LengthUnit::Bohr: return value * BOHR_TO_ANGSTROM;
    case LengthUnit::Alat: return value * alat_scale * BOHR_TO_ANGSTROM;
  }
  return value;
}

} // namespace core
} // namespace zinc
