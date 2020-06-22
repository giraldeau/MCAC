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

*/

#include "spheres/sphere.hpp"
#include "tools/tools.hpp"
#include <iostream>


namespace mcac {
/* Getters */

[[gnu::pure]] size_t Sphere::get_index() const noexcept {
    return index_in_storage;
}
[[gnu::pure]] double Sphere::get_volume() const noexcept {
    return *volume;
}
[[gnu::pure]] double Sphere::get_surface() const noexcept {
    return *surface;
}
[[gnu::pure]] double Sphere::get_radius() const noexcept {
    return *r;
}
[[gnu::pure]] std::array<double, 3> Sphere::get_position() const noexcept {
    return {*x, *y, *z};
}
[[gnu::pure]] std::array<double, 3> Sphere::get_relative_position() const noexcept {
    return {*rx, *ry, *rz};
}
/* Modifiers */
void Sphere::set_label(long value) noexcept {
    agg_label = value;
}
void Sphere::decrease_label() noexcept {
    agg_label--;
}
void Sphere::translate(const std::array<double, 3> &trans) noexcept {
    *x += trans[0];
    *y += trans[1];
    *z += trans[2];
}
void Sphere::relative_translate(std::array<double, 3> trans) noexcept {
    *rx += trans[0];
    *ry += trans[1];
    *rz += trans[2];
}
void Sphere::set_position(std::array<double, 3> newposition) noexcept {

//    *x = periodic_position(_newx,physicalmodel->box_lenght);
//    *y = periodic_position(_newy,physicalmodel->box_lenght);
//    *z = periodic_position(_newz,physicalmodel->box_lenght);

    *x = newposition[0];
    *y = newposition[1];
    *z = newposition[2];
}
void Sphere::set_relative_position(std::array<double, 3> newposition) noexcept {
    *rx = newposition[0];
    *ry = newposition[1];
    *rz = newposition[2];
}
void Sphere::init_val() noexcept {
    init_val({0, 0, 0}, 0.);
}
void Sphere::init_val(std::array<double, 3> newposition, double newr) noexcept {
    setpointers();
    set_position(newposition);
    *r = newr;
    *rx = 0.;
    *ry = 0.;
    *rz = 0.;
    update_vol_and_surf();
}
/* #############################################################################################################
 * ############################################# Sphere growing ################################################
 * #############################################################################################################*/
void Sphere::croissance_surface(double dt) noexcept {
    double new_r = physicalmodel->grow(*r, dt);
    double new_r_2 = new_r * new_r;
    double new_r_3 = new_r_2 * new_r;
    *r = new_r;
    *volume = _volume_factor * new_r_3;
    *surface = _surface_factor * new_r_2;
}
void Sphere::print() const noexcept {
    std::cout << "Printing Sphere " << (index_in_storage) << std::endl;
    if (static_cast<bool>(external_storage)) {
        std::cout << "  With external Storage" << std::endl;
    } else {
        std::cout << "  Without external Storage" << std::endl;
    }
    if (agg_label < 0) {
        std::cout << "  This is a virtual Sphere" << std::endl;
    } else {
        std::cout << "  This Sphere is own by the aggregate " << agg_label << std::endl;
    }
    std::cout << "    Position : " << *x << " " << *y << " " << *z << std::endl;
    std::cout << "    Radius   : " << *r << std::endl;
    std::cout << "    Volume   : " << *volume << std::endl;
    std::cout << "    Surface  : " << *surface << std::endl;
    /*
    double* rx;
    double* ry;
    double* rz;
    */
}
/* #############################################################################################################
 * ##################################### compute_volume and surface of a sphere ########################################
 * #############################################################################################################*/


void Sphere::update_vol_and_surf() noexcept {
    if (agg_label > -1) {
        *volume = _volume_factor * POW_3(*r);
        *surface = _surface_factor * POW_2(*r);
    }
}
}  // namespace mcac

