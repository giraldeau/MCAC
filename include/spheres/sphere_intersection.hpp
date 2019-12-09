#ifndef INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP
#define INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP
#include "spheres/sphere.hpp"
#include <array>


namespace mcac {
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
}  // namespace mcac
#endif //INCLUDE_SPHERES_SPHERE_INTERSECTION_HPP
