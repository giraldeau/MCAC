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
#ifndef INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP
#define INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP
#include "spheres/sphere.hpp"
#include <array>


namespace mcac {
class Intersection {
public:
    double volume_1;
    double volume_2;
    double surface_1;
    double surface_2;
    double dist;
    Intersection(const Sphere &sphere_1, const Sphere &sphere_2, double dist) noexcept;
    Intersection(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
};  // Intersection
}  // namespace mcac
#endif //INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP
