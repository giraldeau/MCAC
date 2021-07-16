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

/*

Sphere.h and Sphere.cpp defines the data storage.

 * Sphere *
 This is an object representing a sphere (!) with its volume and surface
 Data can be stored in an external Aggregat for vectorization purposes

 Beyond managing its properties, it can compute
 - its distance to a point or an another sphere
 - its intersection with an another sphere (volume and surface)
 - detect a collision with an another sphere

 * Aggregat *
 This is container for an aggregat which is an enhanced list of spheres
 Data can be shared between multiple Aggregat

*/

#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include <cmath>


namespace mcac {
/* #############################################################################################################
 * ################################# distance between a sphere and a point #####################################
 * #############################################################################################################*/
[[gnu::pure]]  double distance(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return distance(sphere_1.get_position(), sphere_2.get_position(), sphere_1.physicalmodel->box_lenght);
}
[[gnu::pure]]  double distance_2(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return distance_2(sphere_1.get_position(), sphere_2.get_position(), sphere_1.physicalmodel->box_lenght);
}
[[gnu::pure]]  double relative_distance(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return relative_distance(sphere_1.get_relative_position(), sphere_2.get_relative_position());
}
[[gnu::pure]]  double relative_distance_2(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return relative_distance_2(sphere_1.get_relative_position(), sphere_2.get_relative_position());
}
[[gnu::pure]]  double distance(const std::array<double, 3> &point_1,
                               const std::array<double, 3> &point_2,
                               double box_size) noexcept {
    return std::sqrt(distance_2(point_1, point_2, box_size));
}
[[gnu::pure]]  double relative_distance(const std::array<double, 3> &point_1,
                                        const std::array<double, 3> &point_2) noexcept {
    return std::sqrt(relative_distance_2(point_1, point_2));
}
[[gnu::pure]]  double distance_2(const std::array<double, 3> &point_1, const std::array<double, 3> &point_2,
                                 double box_size) noexcept {
    std::array<double, 3> diff = point_1 - point_2;
    double dx(periodic_distance(diff[0], box_size));
    double dy(periodic_distance(diff[1], box_size));
    double dz(periodic_distance(diff[2], box_size));
    return std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2);
}
[[gnu::pure]]  double relative_distance_2(const std::array<double, 3> &point_1,
                                          const std::array<double, 3> &point_2) noexcept {
    std::array<double, 3> diff = point_1 - point_2;
    return std::pow(diff[0], 2) + std::pow(diff[1], 2) + std::pow(diff[2], 2);
}
[[gnu::pure]]  bool contact(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    //$ Compute signed distance for contact between two spheres
    double distance = distance_2(sphere_1, sphere_2);

    //$ Compute minimum distance for contact
    double dist_contact = std::pow(sphere_1.get_radius() + sphere_2.get_radius(), 2);

    // 1e-28 is for rounding error (1e-14 ^ 2)
    return (distance - dist_contact <= _CONTACT_EPSILON);
}
}  // namespace mcac

