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
 - its intersection with an another sphere (volume and surface) (TODO : NOT CORRECTLY COMPUTED)
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
[[gnu::pure]] SphereContactInfo distance_to_contact_old(const Sphere &sphere_1,
                                                        const Sphere &sphere_2,
                                                        const std::array<double, 3> &displacement,
                                                        const double distance) noexcept {
    std::array<std::array<double, 3>, 3> rot_mat = get_rot_mat(displacement);
    return distance_to_contact_old_r(sphere_1, sphere_2, rot_mat, distance);
}
[[gnu::pure]] SphereContactInfo distance_to_contact_old_r(const Sphere &sphere_1,
                                                          const Sphere &sphere_2,
                                                          const std::array<std::array<double, 3>, 3> &rot_mat,
                                                          const double distance) noexcept {
    /*
     * We use a change of axis system
     * the center is placed on the mobil sphere
     * the x axis is chosen to be the movement vector
     *
     * Thus we have two rotation to perform : one around z and one around y
     *
     * Then a collision is easy to detect :
     * The distance to the x axis must be less than the sum of the radius
     *
     * Moreover we are only interested in collision that can happend in the future,
     * wich is also easy to detect (x coordinate must be positive)
     *
     * Finally we consider 27 possible version of the other sphere in order to take into account
     * any periodicity effect.
    */

    double dist_contact = POW_2(sphere_1.get_radius() + sphere_2.get_radius());
    double dist = distance_2(sphere_1, sphere_2);
    double minval = std::numeric_limits<double>::infinity(); // infinity is too far to care
    if (dist <= dist_contact) {
        // already in contact
        minval = 0;
    } else {
        double box_lenght = sphere_1.physicalmodel->box_lenght;
        std::array<double, 3> oldpos_1{sphere_1.get_position()};
        std::array<double, 3> oldpos_2{sphere_2.get_position()};
        std::array<double, 3> pos{};
        for (size_t l = 0; l < 3; ++l) {
            pos[l] = sum(rot_mat[l] * (oldpos_2 - oldpos_1));
        }
        std::array<double, 3> perx{box_lenght * rot_mat[0][0],
                                   box_lenght * rot_mat[1][0],
                                   box_lenght * rot_mat[2][0],};
        std::array<double, 3> pery{box_lenght * rot_mat[0][1],
                                   box_lenght * rot_mat[1][1],
                                   box_lenght * rot_mat[2][1],};
        std::array<double, 3> perz{box_lenght * rot_mat[0][2],
                                   box_lenght * rot_mat[1][2],
                                   box_lenght * rot_mat[2][2],};
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                for (int k = -1; k <= 1; k++) {
                    std::array<double, 3> tmp{pos + i * perx + j * pery + k * perz};

                    // in the future
                    if (tmp[0] < 0.) {
                        continue;
                    }
                    double dist_1 = POW_2(tmp[1]) + POW_2(tmp[2]);
                    bool end_contact = POW_2(tmp[0] - distance) + dist_1 <= dist_contact;
                    if ((!end_contact) && distance < tmp[0]) {
//                         too much distance to travel
                        continue;
                    }

                    // collision is possible
                    if (dist_1 > dist_contact) {
                        continue;
                    }
                    double sol = dist_contact - dist_1;
                    sol = tmp[0] - sqrt(sol);
                    minval = MIN(minval, sol);
                }
            }
        }
    }
    return SphereContactInfo(minval);
}
[[gnu::pure]] std::array<std::array<double, 3>, 3> get_rot_mat(const std::array<double, 3> &displacement) noexcept {
    double anglez = -atan2(displacement[1], displacement[0]);
    std::array<std::array<double, 3>, 3> rotz{};
    rotz[0] = {{cos(anglez), -sin(anglez), 0}};
    rotz[1] = {{sin(anglez), cos(anglez), 0}};
    rotz[2] = {{0, 0, 1}};
    std::array<double, 3> tmp{};
    for (size_t i = 0; i < tmp.size(); ++i) {
        tmp[i] = sum(rotz[i] * displacement);
    }
    double angley = atan2(tmp[2], tmp[0]);
    std::array<std::array<double, 3>, 3> roty{};
    roty[0] = {{cos(angley), 0, sin(angley)}};
    roty[1] = {{0, 1, 0}};
    roty[2] = {{-sin(angley), 0, cos(angley)}};
    std::array<std::array<double, 3>, 3> matrot{};
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            matrot[i][j] = 0;
            for (size_t k = 0; k < 3; ++k) {
                matrot[i][j] += roty[i][k] * rotz[k][j];
            }
        }
    }
    return matrot;
}
[[gnu::pure]] SphereContactInfo distance_to_contact(const Sphere &sphere_1,
                                                    const Sphere &sphere_2,
                                                    const std::array<double, 3> &displacement_vector,
                                                    const double displacement_distance) noexcept {
    //return distance_to_contact_old(sphere_1, sphere_2, displacement_vector, displacement_distance);
    double dist_contact = sphere_1.get_radius() + sphere_2.get_radius();
    double box_lenght = sphere_1.physicalmodel->box_lenght;
    std::array<double, 3> total_displacement = displacement_vector * displacement_distance;
    std::array<double, 3> pos1 = sphere_1.get_position();
    std::array<double, 3> pos2 = sphere_2.get_position();
    std::array<int, 3> nper{0, 0, 0};
    std::array<double, 3> zone_dimension{0., 0., 0.};

    // apply periodicity to sphere_2
    for (size_t l = 0; l < 3; ++l) {
        double base = fmin(pos1[l], pos1[l] + total_displacement[l]) - dist_contact;
        double end = fmax(pos1[l], pos1[l] + total_displacement[l]) + dist_contact;
        zone_dimension[l] = end - base;
        pos2[l] = fmod((pos2[l] - base), box_lenght);
        if (pos2[l] < 0) {
            pos2[l] += box_lenght;
        }
        pos2[l] += base;

        // we may need to check for multiple periodicity anyway
        nper[l] = static_cast<int>(floor(zone_dimension[l] / box_lenght));
    }
    double res = std::numeric_limits<double>::infinity(); // infinity is too far to care
    double dist_contact_2 = POW_2(dist_contact);
    for (int i = 0; i <= nper[0]; i++) {
        for (int j = 0; j <= nper[1]; j++) {
            for (int k = 0; k <= nper[2]; k++) {
                std::array<double, 3> pos3{pos2[0] + i * box_lenght,
                                           pos2[1] + j * box_lenght,
                                           pos2[2] + k * box_lenght};
                std::array<double, 3> diff = pos3 - pos1;

                // shortcut
                if (fabs(diff[0]) > zone_dimension[0]) {
                    continue;
                }
                if (fabs(diff[1]) > zone_dimension[1]) {
                    continue;
                }
                if (fabs(diff[2]) > zone_dimension[2]) {
                    continue;
                }
                if (relative_distance_2(pos1, pos3) <= dist_contact_2) {
                    // already in contact
                    res = 0.;
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
                double dist_to_axis = POW_2(cross[0]) + POW_2(cross[1]) + POW_2(cross[2]);
                if (dist_to_axis > dist_contact_2) {
                    // too far
                    continue;
                }
                res = fmin(res, proj - sqrt(dist_contact_2 - dist_to_axis));
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
    //std::array<std::array<double, 3>, 3> rot_mat{get_rot_mat(displacement)};
    for (size_t i_sphere_2 = 0; i_sphere_2 < list.size(); i_sphere_2++) {
        SphereContactInfo potential_contact = distance_to_contact(sphere_1, list[i_sphere_2], displacement, distance);
        //SphereContactInfo potential_contact = distance_to_contact_old_r(sphere_1, list[i_sphere_2], rot_mat, distance);
        if (potential_contact < closest_contact) {
            closest_contact.other_sphere = i_sphere_2;
            closest_contact.distance = potential_contact.distance;
        }
    }
    return closest_contact;
}
}  // namespace mcac

