#include "spheres/sphere_intersection.hpp"
#include "spheres/sphere_distance.hpp"
#include "cst.hpp"
#include <cmath>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))
using namespace std;
namespace MCAC {
Intersection::Intersection(const Sphere &sphere_1,
                           const Sphere &sphere_2) noexcept : Intersection(sphere_1, sphere_2,
                                                                           distance(sphere_1, sphere_2)) {
};
Intersection::Intersection(const Sphere &sphere_1, const Sphere &sphere_2,
                           double given_dist) noexcept :
    volume_1{0.},
    volume_2{0.},
    surface_1{0.},
    surface_2{0.} {
    dist = given_dist;
    double radius_1 = sphere_1.get_radius();
    double radius_2 = sphere_2.get_radius();

    //$ Check if they are in contact
    if (dist < radius_1 + radius_2) {
        if (dist >= fabs(radius_1 - radius_2)) {
            //$ get_volume of the intersection is returned
            double h_1 = (POW_2(radius_2) - POW_2((radius_1 - dist))) / (2. * dist);
            double h_2 = (POW_2(radius_1) - POW_2((radius_2 - dist))) / (2. * dist);
            volume_1 = _pi * POW_2(h_1) * (3 * radius_1 - h_1) / 3.;
            volume_2 = _pi * POW_2(h_2) * (3 * radius_2 - h_2) / 3.;
            surface_1 = 2 * _pi * radius_1 * h_1;
            surface_2 = 2 * _pi * radius_2 * h_2;
        }
            //$ Check if one is completely absorbed by the other
        else if (radius_1 < radius_2) {
            //$ Volcal = VolJ
            volume_1 = sphere_1.get_volume();
            surface_1 = sphere_1.get_surface();
        } else // if (radius_2 < radius_1)
        {
            //$ Volcal = Voli
            volume_2 = sphere_2.get_volume();
            surface_2 = sphere_2.get_surface();
        }
    }
}
}  // namespace MCAC

