#pragma once

#include <Eigen/Dense>
#include <array>
#include <optional>

namespace zinc {
namespace core {

using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix3d;

class Lattice {
public:
  Lattice();
  Lattice(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);
  static Lattice from_cell_parameters(const std::array<double, 3>& lengths,
                                      const std::array<double, 3>& angles_deg);

  double volume() const;
  Vector3d frac_to_cart(const Vector3d& frac) const;
  std::optional<Vector3d> cart_to_frac(const Vector3d& cart) const;
  std::optional<std::pair<std::array<double, 3>, std::array<double, 3>>> cell_parameters() const;

  const Matrix3d& matrix() const { return matrix_; }

private:
  Matrix3d matrix_;
  static double angle_deg(const Vector3d& u, const Vector3d& v);
};

} // namespace core
} // namespace zinc
