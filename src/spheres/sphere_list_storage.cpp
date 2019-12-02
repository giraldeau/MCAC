
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



#include "spheres/sphere_list.hpp"
#include <cstdio>
#include <iostream>
#include <utility>


using namespace std;
namespace MCAC {
void SphereList::setpointers() {
    auto newdeb((*Storage)[0].begin());
    auto newfin((*Storage)[0].end());
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
    storage_list<SpheresFields::NFIELD, Sphere>(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(nullptr),
    last_saved(0),
    physicalmodel(nullptr) {
}
SphereList::SphereList(const PhysicalModel &physical_model, size_t size) noexcept :
    storage_list<SpheresFields::NFIELD, Sphere>(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(physical_model, size)),
    last_saved(0),
    physicalmodel(&physical_model) {
    init(physical_model, size);
}
/** Constructor with external storage */
SphereList::SphereList(SphereList &parent, const vector<size_t> &index) noexcept:
    storage_list<SpheresFields::NFIELD, Sphere>(parent, index),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(*parent.physicalmodel, size())),
    last_saved(0),
    physicalmodel(parent.physicalmodel) {
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
SphereList::SphereList(const SphereList &other, SphereList &storage) noexcept:
    storage_list<SpheresFields::NFIELD, Sphere>(other, *this, storage),
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
    storage_list<SpheresFields::NFIELD, Sphere>(move(other)),
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
    if (static_cast<bool>(physicalmodel) && physicalmodel->toBeDestroyed) {
        delete physicalmodel;
    }
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
    if (static_cast<bool>(physicalmodel) && physicalmodel->toBeDestroyed) {
        delete physicalmodel;
    }
    delete writer;
    physicalmodel = other.physicalmodel;
    writer = other.writer;
    last_saved = other.last_saved;
    other.physicalmodel = nullptr;
    other.writer = nullptr;
    storage_list<SpheresFields::NFIELD, Sphere>::operator=(move(other));
    setpointers();
    return *this;
}
}  // namespace MCAC

