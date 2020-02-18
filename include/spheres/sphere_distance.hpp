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
#ifndef INCLUDE_SPHERES_SPHERE_DISTANCE_HPP
#define INCLUDE_SPHERES_SPHERE_DISTANCE_HPP
#include "spheres/sphere.hpp"
#include <array>


namespace mcac {
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
}  // namespace mcac
#endif //INCLUDE_SPHERES_SPHERE_DISTANCE_HPP
