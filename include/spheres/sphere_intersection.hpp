#ifndef INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP_
#define INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP_
#include "spheres/sphere.hpp"
#include <array>


namespace MCAC {
class Intersection {
public:
    double volume_1;
    double volume_2;
    double surface_1;
    double surface_2;
    double dist;
    Intersection(const Sphere &sphere_1, const Sphere &sphere_2, double dist) noexcept;
    Intersection(const Sphere &sphere_1, const Sphere &sphere_2) noexcept;
};  // Intersection
}  // namespace MCAC
#endif //INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP_
