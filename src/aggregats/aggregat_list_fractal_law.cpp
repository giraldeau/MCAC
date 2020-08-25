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
#include "aggregats/aggregat_list.hpp"
#include "tools/tools.hpp"


namespace mcac {
std::tuple<bool, double, double, double> AggregatList::get_instantaneous_fractal_law() const {
    std::vector<double> sizes;
    std::vector<double> dg_over_dps;
    sizes.reserve(size());
    dg_over_dps.reserve(size());
    for (const auto& agg : list) {
        sizes.push_back(double(agg->size()));
        dg_over_dps.push_back(*agg->dg_over_dp);
    }
    return linreg(dg_over_dps, sizes);
}
}// namespace mcac

