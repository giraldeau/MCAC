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
#include "aggregats/aggregat_distance.hpp"
#include "spheres/sphere_contact.hpp"
#include "spheres/sphere_distance.hpp"


namespace mcac {
AggregateContactInfo distance_to_contact(const Aggregate &aggregate_1,
                                         const Aggregate &aggregate_2,
                                         const std::array<double, 3>& direction,
                                         const double distance) noexcept {
    AggregateContactInfo closest_contact; //infinity by default
    //$ For every sphere in the aggregate :
    for (size_t i_sphere = 0; i_sphere < aggregate_1.myspheres.size(); i_sphere++) {
        HalfSphereListContactInfo potential_contact = distance_to_contact(aggregate_1.myspheres[i_sphere],
                                                                          aggregate_2.myspheres,
                                                                          direction,
                                                                          distance);
        if (potential_contact < closest_contact) {
            closest_contact.moving_sphere = i_sphere;
            closest_contact.moving_aggregate = aggregate_1.get_label();
            closest_contact.other_aggregate = aggregate_2.get_label();
            closest_contact.distance = potential_contact.distance;
        }
    }
    return closest_contact;
}
bool contact(const Sphere &sphere,
             const Aggregate &aggregate) noexcept {
    Sphere SphereAggregate(aggregate);
    if (!contact(sphere, SphereAggregate)) {
        return false;
    }
    //$ Loop on all the spheres of the other aggregate
    for (const auto &othersphere : aggregate.myspheres) {
        if (!contact(sphere, *othersphere)) {
            return true;
        }
    }
    return false;
}
bool contact(const Aggregate &aggregate_1,
             const Aggregate &aggregate_2) noexcept {
    Sphere SphereMe(aggregate_1);
    Sphere SphereOther(aggregate_2);
    if (!contact(SphereMe, SphereOther)) {
        return false;
    }
    //$ Loop on all the spheres of the other aggregate
    for (const auto &othersphere : aggregate_1.myspheres) {
        if (contact(*othersphere, aggregate_2)) {
            return true;
        }
    }
    return false;
}
}  // namespace mcac
