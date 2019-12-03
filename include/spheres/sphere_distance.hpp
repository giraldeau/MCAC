#ifndef INCLUDE_SPHERES_SPHERE_DISTANCE_HPP_
#define INCLUDE_SPHERES_SPHERE_DISTANCE_HPP_
#include "spheres/sphere.hpp"
#include <array>


namespace MCAC {
double distance(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
double distance_2(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
double relative_distance(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
double relative_distance_2(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
double distance(const std::array<double, 3> &point_1,
                const std::array<double, 3> &point_2,
                double box_size) noexcept;
double relative_distance(const std::array<double, 3> &point_1, const std::array<double, 3> &point_2) noexcept;
double distance_2(const std::array<double, 3> &point_1, const std::array<double, 3> &point_2,
                  double box_size) noexcept;
double relative_distance_2(const std::array<double, 3> &point_1, const std::array<double, 3> &point_2) noexcept;
bool contact(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
}  // namespace MCAC
#endif //INCLUDE_SPHERES_SPHERE_DISTANCE_HPP_
