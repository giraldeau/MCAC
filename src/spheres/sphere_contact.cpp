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



#include "spheres/sphere_contact.hpp"
#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include <cmath>


namespace mcac {
[[gnu::pure]] SphereContactInfo distance_to_contact(const Sphere &sphere_1,
                                                    const Sphere &sphere_2,
                                                    const std::array<double, 3> &displacement_vector,
                                                    const double displacement_distance) noexcept {
    double dist_contact = sphere_1.get_radius() + sphere_2.get_radius();
    double dist_contact_2 = std::pow(dist_contact, 2);
    double box_lenght = sphere_1.physicalmodel->box_lenght;
    std::array<double, 3> pos1 = sphere_1.get_position();
    std::array<double, 3> pos2 = sphere_2.get_position();

    if (distance_2(pos1, pos2, box_lenght) <= dist_contact_2) {
        // already in contact
        return SphereContactInfo(0.);
    }

    std::array<double, 3> total_displacement = displacement_vector * displacement_distance;

    std::array<int, 3> nper{0, 0, 0};
    std::array<double, 3> zone_dimension{0., 0., 0.};

    // apply periodicity to sphere_2
    for (size_t l = 0; l < 3; ++l) {
        double base = std::min(pos1[l], pos1[l] + total_displacement[l]) - dist_contact;
        double end = std::max(pos1[l], pos1[l] + total_displacement[l]) + dist_contact;
        zone_dimension[l] = end - base;
        pos2[l] = std::fmod((pos2[l] - base), box_lenght);
        if (pos2[l] < 0) {
            pos2[l] += box_lenght;
        }
        pos2[l] += base;

        // we may need to check for multiple periodicity anyway
        nper[l] = static_cast<int>(std::floor(zone_dimension[l] / box_lenght));
    }
    double res = std::numeric_limits<double>::infinity(); // infinity is too far to care
    for (int i = 0; i <= nper[0]; i++) {
        for (int j = 0; j <= nper[1]; j++) {
            for (int k = 0; k <= nper[2]; k++) {
                std::array<double, 3> pos3{pos2[0] + i * box_lenght,
                                           pos2[1] + j * box_lenght,
                                           pos2[2] + k * box_lenght};
                std::array<double, 3> diff = pos3 - pos1;

                // shortcut
                if (std::abs(diff[0]) > zone_dimension[0]) {
                    continue;
                }
                if (std::abs(diff[1]) > zone_dimension[1]) {
                    continue;
                }
                if (std::abs(diff[2]) > zone_dimension[2]) {
                    continue;
                }
                double proj = diff[0] * displacement_vector[0] +
                              diff[1] * displacement_vector[1] +
                              diff[2] * displacement_vector[2];
                if (proj < 0) {
                    // contact in the past
                    continue;
                }
                bool end_contact = relative_distance_2(pos1 + total_displacement, pos3) <= dist_contact_2;
                if ((!end_contact) && displacement_distance < proj) {
                    // too much distance to travel
                    continue;
                }
                std::array<double, 3> cross{diff[1] * displacement_vector[2] - diff[2] * displacement_vector[1],
                                            diff[2] * displacement_vector[0] - diff[0] * displacement_vector[2],
                                            diff[0] * displacement_vector[1] - diff[1] * displacement_vector[0]};
                double dist_to_axis = std::pow(cross[0], 2) + std::pow(cross[1], 2) + std::pow(cross[2], 2);
                if (dist_to_axis > dist_contact_2) {
                    // too far
                    continue;
                }
                res = std::min(res, proj - std::sqrt(dist_contact_2 - dist_to_axis));
            }
        }
    }
    return SphereContactInfo(res);
}
[[gnu::pure]] HalfSphereListContactInfo distance_to_contact(const Sphere &sphere_1,
                                                            const SphereList &list,
                                                            const std::array<double, 3> &displacement,
                                                            const double distance) noexcept {
    HalfSphereListContactInfo closest_contact; //infinity by default
    for (size_t i_sphere_2 = 0; i_sphere_2 < list.size(); i_sphere_2++) {
        SphereContactInfo potential_contact = distance_to_contact(sphere_1, *list[i_sphere_2], displacement, distance);
        if (potential_contact < closest_contact) {
            closest_contact.other_sphere = std::weak_ptr<Sphere>(list[i_sphere_2]);
            closest_contact.distance = potential_contact.distance;
        }
    }
    return closest_contact;
}
}  // namespace mcac

