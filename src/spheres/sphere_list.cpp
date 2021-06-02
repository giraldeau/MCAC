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
#include "spheres/sphere.hpp"
#include "spheres/sphere_list.hpp"
#include <iostream>



/* #############################################################################################################
 * #################################                                       #####################################
 * #################################              AGREGATE                 #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/
namespace mcac {
void SphereList::init(const PhysicalModel &physical_model, size_t size) {
    physicalmodel = &physical_model;
    writer = std::make_unique<ThreadedIO>(physical_model.output_dir / "Spheres", physical_model, size);
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>::init(size, *this);
    setpointers();
}
void SphereList::decrease_label() noexcept {
    for (const auto& mysphere : list) {
        mysphere->decrease_label();
    }
}
void SphereList::print() const {
    std::cout << "Printing list of " << size() << " Sphere" << std::endl;
    if (static_cast<bool>(external_storage)) {
        std::cout << "  With external Storage" << std::endl;
    } else {
        std::cout << "  Without external Storage" << std::endl;
    }
    for (const auto& s : list) {
        s->print();
    }
}
/* #############################################################################################################
 * ########################################### grow all spheres ################################################
 * #############################################################################################################*/
void SphereList::croissance_surface(double dt) noexcept {
    for (const auto& mysphere : list) {
        mysphere->croissance_surface(dt);
    }
}
}  // namespace mcac

