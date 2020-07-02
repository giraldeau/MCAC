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
#include "spheres/sphere_distance.hpp"
#include "spheres/sphere_intersection.hpp"
#include "tools/tools.hpp"
#include <cmath>


namespace mcac {
Intersection::Intersection(const Sphere &sphere_1,
                           const Sphere &sphere_2) noexcept : Intersection(sphere_1, sphere_2,
                                                                           distance(sphere_1, sphere_2)) {
}
Intersection::Intersection(const Sphere &sphere_1, const Sphere &sphere_2,
                           double given_dist) noexcept :
    volume_1{0.},
    volume_2{0.},
    surface_1{0.},
    surface_2{0.},
    dist{given_dist} {
    if (given_dist <= 0) {
        return;
    }
    double radius_1 = sphere_1.get_radius();
    double radius_2 = sphere_2.get_radius();

    //$ Check if they are in contact
    if (dist < radius_1 + radius_2) {
        if (dist >= std::fdim(radius_1, radius_2)) {
            //$ get_volume of the intersection is returned
            double h_1 = (std::pow(radius_2, 2) - std::pow((radius_1 - dist), 2)) / (2. * dist);
            double h_2 = (std::pow(radius_1, 2) - std::pow((radius_2 - dist), 2)) / (2. * dist);
            volume_1 = _pi * std::pow(h_1, 2) * (3 * radius_1 - h_1) / 3.;
            volume_2 = _pi * std::pow(h_2, 2) * (3 * radius_2 - h_2) / 3.;
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
}  // namespace mcac

