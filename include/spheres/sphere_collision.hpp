#ifndef INCLUDE_SPHERES_SPHERE_COLLISION_HPP_
#define INCLUDE_SPHERES_SPHERE_COLLISION_HPP_
#include "spheres/sphere.hpp"
#include "spheres/sphere_list.hpp"
#include <array>


namespace MCAC {
std::pair<bool, double> collision_r(const Sphere &sphere_1, const Sphere &sphere_2,
                                    std::array<std::array<double, 3>, 3> rot_mat) noexcept;
std::pair<bool, double> collision(const Sphere &sphere_1,
                                  const Sphere &sphere_2,
                                  std::array<double, 3> displacement) noexcept;
std::vector<double> collisions(const Sphere &sphere, const SphereList &list, std::array<double, 3> displacement) noexcept;
std::array<std::array<double, 3>, 3> get_rot_mat(std::array<double, 3> displacement) noexcept;
}  // namespace MCAC
#endif //INCLUDE_SPHERES_SPHERE_COLLISION_HPP_
