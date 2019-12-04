
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


using namespace std;
namespace MCAC {
/* Getters */

[[gnu::pure]] double Sphere::get_volume() const noexcept {
    return *volume;
}
[[gnu::pure]] double Sphere::get_surface() const noexcept {
    return *surface;
}
[[gnu::pure]] double Sphere::get_radius() const noexcept {
    return *r;
}
[[gnu::pure]] array<double, 3> Sphere::get_position() const noexcept {
    return {*x, *y, *z};
}
[[gnu::pure]] array<double, 3> Sphere::get_relative_position() const noexcept {
    return {*rx, *ry, *rz};
}
/* Modifiers */
void Sphere::set_label(long value) noexcept {
    agg_label = value;
}
void Sphere::decrease_label() noexcept {
    agg_label--;
}
void Sphere::translate(const array<double, 3> &trans) noexcept {
    *x += trans[0];
    *y += trans[1];
    *z += trans[2];
}
void Sphere::relative_translate(array<double, 3> trans) noexcept {
    *rx += trans[0];
    *ry += trans[1];
    *rz += trans[2];
}
void Sphere::set_position(array<double, 3> newposition) noexcept {
    /*
    *x = periodicPosition(_newx,physicalmodel->L);
    *y = periodicPosition(_newy,physicalmodel->L);
    *z = periodicPosition(_newz,physicalmodel->L);
    */
    *x = newposition[0];
    *y = newposition[1];
    *z = newposition[2];
}
void Sphere::init_val() noexcept {
    init_val({0, 0, 0}, 0.);
}
void Sphere::init_val(array<double, 3> newposition, double newr) noexcept {
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
    double new_r = physicalmodel->Grow(*r, dt);
    double new_r_2 = new_r * new_r;
    double new_r_3 = new_r_2 * new_r;
    *r = new_r;
    *volume = _facvol * new_r_3;
    *surface = _facsurf * new_r_2;
}
void Sphere::print() const noexcept {
    cout << "Printing Sphere " << (indexInStorage) << endl;
    if (!static_cast<bool>(external_storage)) {
        cout << "  With external Storage" << endl;
    } else {
        cout << "  Without external Storage" << endl;
    }
    if (agg_label < 0) {
        cout << "  This is a virtual Sphere" << endl;
    } else {
        cout << "  This Sphere is own by the aggregate " << agg_label << endl;
    }
    cout << "    Position : " << *x << " " << *y << " " << *z << endl;
    cout << "    Radius   : " << *r << endl;
    cout << "    Volume   : " << *volume << endl;
    cout << "    Surface  : " << *surface << endl;
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
        *volume = _facvol * POW_3(*r);
        *surface = _facsurf * POW_2(*r);
    }
}
}  // namespace MCAC

