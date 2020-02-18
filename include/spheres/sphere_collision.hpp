/*
 * MCAC
 * Copyright (C) 2020 CORIA
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
