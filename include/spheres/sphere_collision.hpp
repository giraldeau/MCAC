#ifndef INCLUDE_SPHERES_SPHERE_COLLISION_HPP
#define INCLUDE_SPHERES_SPHERE_COLLISION_HPP
#include "spheres/sphere.hpp"
#include "spheres/sphere_list.hpp"
#include <array>


namespace mcac {
std::pair<bool, double> sphere_collision_r(const Sphere &sphere_1, const Sphere &sphere_2,
                                           const std::array<std::array<double, 3>, 3> &rot_mat) noexcept;
std::pair<bool, double> sphere_collision(const Sphere &sphere_1,
                                         const Sphere &sphere_2,
                                         const std::array<double, 3> &displacement) noexcept;
std::vector<double> sphere_collisions(const Sphere &sphere,
                                      const SphereList &list,
                                      const std::array<double, 3> &displacement) noexcept;
std::array<std::array<double, 3>, 3> get_rot_mat(const std::array<double, 3> &displacement) noexcept;
}  // namespace mcac
#endif //INCLUDE_SPHERES_SPHERE_COLLISION_HPP
