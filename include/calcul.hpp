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
#ifndef INCLUDE_CALCUL_HPP
#define INCLUDE_CALCUL_HPP 1
#include "aggregats/aggregat_list.hpp"
#include "physical_model/physical_model.hpp"


namespace mcac {
void calcul(PhysicalModel &, AggregatList &);
}  // namespace mcac

#endif //INCLUDE_CALCUL_HPP
