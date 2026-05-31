#include "zinc/core/lattice.h"
#include "zinc/core/constants.h"
#include <cmath>
#include <stdexcept>

namespace zinc {
namespace core {

Lattice::Lattice() : matrix_(Matrix3d::Identity()) {}

Lattice::Lattice(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3) {
  matrix_.row(0) = v1;
  matrix_.row(1) = v2;
  matrix_.row(2) = v3;
}

Lattice Lattice::from_cell_parameters(const std::array<double, 3>& lengths,
                                      const std::array<double, 3>& angles_deg) {
  double a = lengths[0], b = lengths[1], c = lengths[2];
  double alpha = angles_deg[0] * PI / 180.0;
  double beta = angles_deg[1] * PI / 180.0;
  double gamma = angles_deg[2] * PI / 180.0;

  double sin_gamma = std::sin(gamma);
  double cos_alpha = std::cos(alpha);
  double cos_beta = std::cos(beta);
  double cos_gamma = std::cos(gamma);

  Vector3d v1(a, 0.0, 0.0);
  Vector3d v2(b * cos_gamma, b * sin_gamma, 0.0);
  double v3_x = c * cos_beta;
  double v3_y = c * (cos_alpha - cos_beta * cos_gamma) / std::max(sin_gamma, 1e-12);
  double v3_z_sq = c * c - v3_x * v3_x - v3_y * v3_y;
  double v3_z = (v3_z_sq > 0.0) ? std::sqrt(v3_z_sq) : 0.0;
  Vector3d v3(v3_x, v3_y, v3_z);

  return Lattice(v1, v2, v3);
}

double Lattice::volume() const {
  return std::abs(matrix_.determinant());
}

Vector3d Lattice::frac_to_cart(const Vector3d& frac) const {
  return matrix_.transpose() * frac;
}

std::optional<Vector3d> Lattice::cart_to_frac(const Vector3d& cart) const {
  Matrix3d inv = matrix_.transpose().inverse();
  if (std::abs(inv.determinant()) < 1e-10) return std::nullopt;
  return inv * cart;
}

std::optional<std::pair<std::array<double, 3>, std::array<double, 3>>> Lattice::cell_parameters() const {
  Vector3d r0 = matrix_.row(0), r1 = matrix_.row(1), r2 = matrix_.row(2);
  double a = r0.norm(), b = r1.norm(), c = r2.norm();
  if (a < 1e-10 || b < 1e-10 || c < 1e-10) return std::nullopt;
  return {{{a, b, c}, {angle_deg(r1, r2), angle_deg(r0, r2), angle_deg(r0, r1)}}};
}

double Lattice::angle_deg(const Vector3d& u, const Vector3d& v) {
  double denom = u.norm() * v.norm();
  if (denom < 1e-12) return 0.0;
  double cos_theta = u.dot(v) / denom;
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
  return std::acos(cos_theta) * 180.0 / PI;
}

} // namespace core
} // namespace zinc
