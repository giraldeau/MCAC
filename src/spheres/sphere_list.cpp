
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
#include "spheres/sphere.hpp"
#include <iostream>



/* #############################################################################################################
 * #################################                                       #####################################
 * #################################              AGREGATE                 #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/
using namespace std;
namespace MCAC {
void SphereList::init(const PhysicalModel &physical_model, size_t size) {
    delete writer;
    physicalmodel = &physical_model;
    writer = new ThreadedIO(physical_model, size);
    storage_list<SpheresFields::SPHERE_NFIELDS, Sphere>::Init(size, *this);
    setpointers();
}
void SphereList::decrease_label() noexcept {
    for (Sphere *mysphere : list) {
        mysphere->decrease_label();
    }
}
void SphereList::print() const {
    cout << "Printing list of " << size() << " Sphere" << endl;
    if (static_cast<bool>(external_storage)) {
        cout << "  With external Storage" << endl;
    } else {
        cout << "  Without external Storage" << endl;
    }
    for (const Sphere *s : list) {
        s->print();
    }
}
/* #############################################################################################################
 * ########################################### Grow all spheres ################################################
 * #############################################################################################################*/
void SphereList::croissance_surface(double dt) noexcept {
    for (Sphere *mysphere : list) {
        mysphere->croissance_surface(dt);
    }
}
}  // namespace MCAC

