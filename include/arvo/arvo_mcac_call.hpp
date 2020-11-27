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
#ifndef INCLUDE_ARVO_MCAC_CALL_HPP
#define INCLUDE_ARVO_MCAC_CALL_HPP 1
#ifdef WITH_ARVO
#include "spheres/sphere_list.hpp"
#include "arvo/arvo.hpp"
#include <vector>
#include <unordered_map>


namespace mcac {
std::pair<std::vector<double>, std::vector<double>>
arvo_call(const SphereList &spherelist,
          const std::vector<std::unordered_map<size_t, double>> &distances);
}
#endif
#endif //INCLUDE_ARVO_MCAC_CALL_HPP
