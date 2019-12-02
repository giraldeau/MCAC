
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



#include "spheres/sphere_distance.hpp"
#include <cmath>
#include <iomanip>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
using namespace std;
namespace MCAC {
/* #############################################################################################################
 * ################################# Distance between a sphere and a point #####################################
 * #############################################################################################################*/
[[gnu::pure]]  double distance(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return distance(sphere_1.get_position(), sphere_2.get_position(), sphere_1.physicalmodel->L);
}
[[gnu::pure]]  double distance_2(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return distance_2(sphere_1.get_position(), sphere_2.get_position(), sphere_1.physicalmodel->L);
}
[[gnu::pure]]  double relative_distance(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return relative_distance(sphere_1.get_relative_position(), sphere_2.get_relative_position());
}
[[gnu::pure]]  double relative_distance_2(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    return relative_distance_2(sphere_1.get_relative_position(), sphere_2.get_relative_position());
}
[[gnu::pure]]  double distance(array<double, 3> point_1,
                               array<double, 3> point_2,
                               double box_size) noexcept {
    return sqrt(distance_2(point_1, point_2, box_size));
}
[[gnu::pure]]  double relative_distance(array<double, 3> point_1,
                                        array<double, 3> point_2) noexcept {
    return sqrt(relative_distance_2(point_1, point_2));
}
[[gnu::pure]]  double distance_2(array<double, 3> point_1, array<double, 3> point_2,
                                 double box_size) noexcept {
    double dx(periodicDistance((point_1[0] - point_2[0]), box_size));
    double dy(periodicDistance((point_1[1] - point_2[1]), box_size));
    double dz(periodicDistance((point_1[2] - point_2[2]), box_size));
    return POW_2(dx) + POW_2(dy) + POW_2(dz);
}
[[gnu::pure]]  double relative_distance_2(array<double, 3> point_1,
                                          array<double, 3> point_2) noexcept {
    return POW_2(point_1[0] - point_2[0])
           + POW_2(point_1[1] - point_2[1])
           + POW_2(point_1[2] - point_2[2]);
}
[[gnu::pure]]  bool contact(const Sphere &sphere_1, const Sphere &sphere_2) noexcept {
    //$ Compute signed distance for contact between two spheres
    double distance = distance_2(sphere_1, sphere_2);

    //$ Compute minimum distance for contact
    double dist_contact = POW_2(sphere_1.get_radius() + sphere_2.get_radius());

    // 1e-28 is for rounding error (1e-14 ^ 2)
    return (distance - dist_contact <= 1e-28);
}
}  // namespace MCAC
