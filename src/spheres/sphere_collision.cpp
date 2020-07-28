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



#include "spheres/sphere_collision.hpp"
#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include <cmath>


namespace mcac {
[[gnu::pure]] std::pair<bool, double> sphere_collision(const Sphere &sphere_1, const Sphere &sphere_2,
                                                       const std::array<double, 3> &displacement) noexcept {
    std::array<std::array<double, 3>, 3> rot_mat = get_rot_mat(displacement);
    return sphere_collision_r(sphere_1, sphere_2, rot_mat);
}
[[gnu::pure]] std::pair<bool, double> sphere_collision_r(const Sphere &sphere_1, const Sphere &sphere_2,
                                                         const std::array<std::array<double, 3>, 3> &rot_mat) noexcept {
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

    double dist_contact = std::pow(sphere_1.get_radius() + sphere_2.get_radius(), 2);
    double dist = distance_2(sphere_1, sphere_2);
    double minval = 0.;
    bool collision = true;
    if (dist > dist_contact) {
        double box_lenght = sphere_1.physicalmodel->box_lenght;
        collision = false;
        minval = 10 * box_lenght;
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
                    double dist_1 = std::pow(tmp[1], 2) + std::pow(tmp[2], 2);

                    // collision is possible
                    if (dist_1 > dist_contact) {
                        continue;
                    }
                    collision = true;
                    double sol = dist_contact - dist_1;
                    sol = tmp[0] - std::sqrt(sol);
                    minval = std::min(minval, sol);
                }
            }
        }
    }
    std::pair<bool, double> result = {collision, minval};
    return result;
}
[[gnu::pure]]  std::vector<double> sphere_collisions(const Sphere &sphere_1,
                                                     const SphereList &list,
                                                     const std::array<double, 3> &displacement) noexcept {
    std::array<std::array<double, 3>, 3> rot_mat{get_rot_mat(displacement)};
    std::vector<double> dist_to_collision;
    for (const Sphere *sphere_2 : list) {
        std::pair<bool, double> suspect = sphere_collision_r(sphere_1, *sphere_2, rot_mat);
        if (suspect.first) {
            dist_to_collision.push_back(suspect.second);
        }
    }
    return dist_to_collision;
}
[[gnu::pure]] std::array<std::array<double, 3>, 3> get_rot_mat(const std::array<double, 3> &displacement) noexcept {
    double anglez = -std::atan2(displacement[1], displacement[0]);
    std::array<std::array<double, 3>, 3> rotz{};
    rotz[0] = {{std::cos(anglez), -std::sin(anglez), 0}};
    rotz[1] = {{std::sin(anglez), std::cos(anglez), 0}};
    rotz[2] = {{0, 0, 1}};
    std::array<double, 3> tmp{};
    for (size_t i = 0; i < tmp.size(); ++i) {
        tmp[i] = sum(rotz[i] * displacement);
    }
    double angley = std::atan2(tmp[2], tmp[0]);
    std::array<std::array<double, 3>, 3> roty{};
    roty[0] = {{std::cos(angley), 0, std::sin(angley)}};
    roty[1] = {{0, 1, 0}};
    roty[2] = {{-std::sin(angley), 0, std::cos(angley)}};
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
}  // namespace mcac

