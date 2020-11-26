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
#ifndef INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP
#define INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP 1
#include "aggregat.hpp"
#include "tools/contact_info.hpp"
#include <array>


namespace mcac {
AggregateContactInfo distance_to_contact(const std::shared_ptr<Aggregate> &aggregate_1,
                                         const std::shared_ptr<Aggregate> &aggregate_2,
                                         const std::array<double, 3> &direction,
                                         const double distance) noexcept;
bool contact(const Aggregate &aggregate_1, const Aggregate &aggregate_2) noexcept;
bool contact(const Sphere &sphere, const Aggregate &aggregate) noexcept;
}  // namespace mcac


#endif //INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP
