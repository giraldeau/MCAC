#include "aggregats/aggregat_distance.hpp"
#include "spheres/sphere_collision.hpp"
#include "spheres/sphere_distance.hpp"
#include <limits>


namespace mcac {
double distance(const Aggregate &aggregate_1,
                const Aggregate &aggregate_2,
                std::array<double, 3> direction) noexcept {
    double mindist(std::numeric_limits<double>::infinity());

    //$ For every sphere in the aggregate :
    for (const auto& mysphere: aggregate_1.myspheres) {
        std::vector<double> dists = sphere_collisions(*mysphere, aggregate_2.myspheres, direction);
        if (!dists.empty()) {
            double lmindist = *min_element(dists.begin(), dists.end());
            if (lmindist <= mindist) {
                mindist = lmindist;
            }
        }
    }
    return mindist;
}
bool contact(const Aggregate &aggregate_1,
             const Aggregate &aggregate_2) noexcept {
    /*
        Sphere SphereMe(aggregate_1);
        Sphere SphereOther(aggregate_2);

        if (! SphereMe.contact(SphereOther))
            return false;
    */
    //$ Loop on all the spheres of the other aggregate
    for (const auto& othersphere : aggregate_1.myspheres) {
        //$ For every sphere in the aggregate :
        for (const auto& mysphere: aggregate_2.myspheres) {
            if (contact(*mysphere, *othersphere)) {
                return true;
            }
        }
    }
    return false;
}
}  // namespace mcac
