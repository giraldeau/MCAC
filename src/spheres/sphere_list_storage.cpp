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

#include "spheres/sphere.hpp"
#include "spheres/sphere_list.hpp"
#include <iostream>


namespace mcac {
void SphereList::setpointers() {
    auto newdeb((*storage)[0].begin());
    auto newfin((*storage)[0].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin)) {
        return;
    }
    for (Sphere *mysphere : list) {
        mysphere->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}
/** Default constructor in local storage */
SphereList::SphereList() noexcept :
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(nullptr),
    last_saved(0),
    physicalmodel(nullptr) {
}
SphereList::SphereList(const PhysicalModel &physical_model, size_t size) noexcept :
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(physical_model, size)),
    last_saved(0),
    physicalmodel(&physical_model) {
    init(physical_model, size);
}
/** Constructor with external storage */
SphereList::SphereList(SphereList *parent, const std::vector<size_t> &index) noexcept:
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>(*parent, index),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(*parent->physicalmodel, size())),
    last_saved(0),
    physicalmodel(parent->physicalmodel) {
    setpointers();
}
/** Copy constructor */
//ListSphere::ListSphere(const ListSphere& other)noexcept:
//    storage_list<SpheresFields::NFIELD,Sphere>(other,*this),
//    ptr_deb(nullptr),
//    ptr_fin(nullptr),
//    writer (new ThreadedIO(*other.physicalmodel, size())),
//    last_saved(other.last_saved),
//    physicalmodel(other.physicalmodel)
//{
//    for (Sphere* s: list)
//    {
//        s->physicalmodel = physicalmodel;
//    }
//    setpointers();
//}
//
SphereList::SphereList(const SphereList &other, SphereList *sphere_list) noexcept:
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>(other, *this, *sphere_list),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(*other.physicalmodel, size())),
    last_saved(other.last_saved),
    physicalmodel(other.physicalmodel) {
    for (Sphere *s: list) {
        s->physicalmodel = physicalmodel;
    }
    setpointers();
}
/** Move constructor */
SphereList::SphereList(SphereList &&other) noexcept:
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>(std::move(other)),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(other.writer),
    last_saved(other.last_saved),
    physicalmodel(other.physicalmodel) {
    other.writer = nullptr;
    setpointers();
}
/** Destructor */
SphereList::~SphereList() noexcept {
    delete writer;
}

/** Copy assignment operator */
//ListSphere& ListSphere::operator= (const ListSphere& other) noexcept
//{
//    ListSphere tmp(other);      // re-use copy-constructor
//    *this = std::move(tmp);     // re-use move-assignment
//    setpointers();
//    return *this;
//}

/** Move assignment operator */
SphereList &SphereList::operator=(SphereList &&other) noexcept {
    physicalmodel = other.physicalmodel;
    delete writer;
    writer = other.writer;
    last_saved = other.last_saved;
    other.physicalmodel = nullptr;
    other.writer = nullptr;
    ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere>::operator=(std::move(other));
    setpointers();
    return *this;
}
}  // namespace mcac

