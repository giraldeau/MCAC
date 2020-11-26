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

#include "aggregats/aggregat.hpp"
#include "spheres/sphere.hpp"


namespace mcac {
void Sphere::setpointers() {
    x = &(*storage)[SpheresFields::SPHERE_X][index_in_storage];
    y = &(*storage)[SpheresFields::SPHERE_Y][index_in_storage];
    z = &(*storage)[SpheresFields::SPHERE_Z][index_in_storage];
    r = &(*storage)[SpheresFields::SPHERE_R][index_in_storage];
    volume = &(*storage)[SpheresFields::SPHERE_VOLUME][index_in_storage];
    surface = &(*storage)[SpheresFields::SPHERE_SURFACE][index_in_storage];
    rx = &(*storage)[SpheresFields::SPHERE_RX][index_in_storage];
    ry = &(*storage)[SpheresFields::SPHERE_RY][index_in_storage];
    rz = &(*storage)[SpheresFields::SPHERE_RZ][index_in_storage];
}
/** Default constructor in local storage */
Sphere::Sphere() noexcept:
    ElemStorage<SpheresFields::SPHERE_NFIELDS, SphereList>(),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    agg_label(-1),
    electric_charge(0),
    physicalmodel(nullptr) {
    init_val();
}
Sphere::Sphere(const PhysicalModel &physical_model) noexcept:
    ElemStorage<SpheresFields::SPHERE_NFIELDS, SphereList>(),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    agg_label(-1),
    electric_charge(0),
    physicalmodel(&physical_model) {
    init_val();
}
/** Constructor in local storage with initialization */
Sphere::Sphere(const PhysicalModel &physical_model,
               const std::array<double, 3> &newposition,
               double newr) noexcept:
    Sphere(physical_model) {
    init_val(newposition, newr);
}
/** Constructor with external storage */
Sphere::Sphere(SphereList *aggregat, size_t id) noexcept:
    ElemStorage<SpheresFields::SPHERE_NFIELDS, SphereList>(*aggregat, id),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    agg_label(-1),
    electric_charge(0),
    physicalmodel(aggregat->physicalmodel) {
    init_val();
    external_storage->setpointers();
}
Sphere::Sphere(const Aggregate &aggregate) noexcept : Sphere() {
    physicalmodel = aggregate.physicalmodel;
    init_val(aggregate.get_position(), aggregate.get_rmax());
}
/** Destructor */
Sphere::~Sphere() noexcept {
}
/** Copy constructor */
Sphere::Sphere(const Sphere &other, SphereList *aggregat, size_t id) noexcept:
    ElemStorage<SpheresFields::SPHERE_NFIELDS, SphereList>(*aggregat, id),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    agg_label(long(id)),
    electric_charge(0),
    physicalmodel(other.physicalmodel) {
    setpointers();
}
//Sphere::Sphere(const Sphere &other) noexcept:
//    storage_elem<SpheresFields::NFIELD, ListSphere>(other),
//    x(nullptr),
//    y(nullptr),
//    z(nullptr),
//    r(nullptr),
//    rx(nullptr),
//    ry(nullptr),
//    rz(nullptr),
//    volume(nullptr),
//    surface(nullptr),
//    agg_label(other.agg_label),
//    physicalmodel(other.physicalmodel) {
//    setpointers();
//}
///** Move constructor */
//Sphere::Sphere(Sphere &&other) noexcept : /* noexcept needed to enable optimizations in containers */
//    storage_elem<SpheresFields::NFIELD, ListSphere>(std::move(other)),
//    x(nullptr),
//    y(nullptr),
//    z(nullptr),
//    r(nullptr),
//    rx(nullptr),
//    ry(nullptr),
//    rz(nullptr),
//    volume(nullptr),
//    surface(nullptr),
//    physicalmodel(other.physicalmodel),
//    agg_label(other.agg_label) {
//    setpointers();
//    other.x = nullptr;
//    other.y = nullptr;
//    other.z = nullptr;
//    other.r = nullptr;
//    other.rx = nullptr;
//    other.ry = nullptr;
//    other.rz = nullptr;
//    other.volume = nullptr;
//    other.surface = nullptr;
//    other.agg_label = -1;
//}
///** Copy assignment operator */
//Sphere &Sphere::operator=(const Sphere &other) noexcept{
//    Sphere tmp(other);      // re-use copy-constructor
//    *this = std::move(tmp); // re-use move-assignment
//    return *this;
//}
///** Move assignment operator */
//Sphere &Sphere::operator=(Sphere &&other) noexcept{
//    if (physicalmodel && physicalmodel->toBeDestroyed) {
//        delete physicalmodel;
//    }
//    physicalmodel = other.physicalmodel;
//    agg_label = other.agg_label;
//    other.agg_label = -1;
//    other.physicalmodel = nullptr;
//    storage_elem<SpheresFields::NFIELD, ListSphere>::operator=(move(other));
//    setpointers();
//    return *this;
//}
}  // namespace mcac

