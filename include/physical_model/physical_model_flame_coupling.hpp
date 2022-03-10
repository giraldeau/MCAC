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
#ifndef INCLUDE_PHYSICAL_MODEL_FLAME_COUPLING_HPP
#define INCLUDE_PHYSICAL_MODEL_FLAME_COUPLING_HPP 1
#include "constants.hpp"
#include <vector>
#include <ctime>


namespace mcac {
/********************************************************************************
* Class for saving the time resolved info of the flame
********************************************************************************/
class FlameCoupling {
private:
public:
    std::vector<double> t_res;                   ///< residence time [s]
    std::vector<double> Temp;                    ///< flame temperature [K]
    std::vector<double> fv;                      ///< soot volume fraction [-]
    std::vector<double> dp_gav;                  ///< primary particle geometric avg. diameter [m]
    std::vector<double> dp_gSTD;                 ///< primary particle geometric std.
    std::vector<double> J_sg;                    ///< surface growth flux [kg/m^3/s]
    std::vector<double> J_nucl;                  ///< nucleation mas flux [kg/m^3/s]
    /// Constructor
    FlameCoupling();
    explicit FlameCoupling(const std::string &flame_file);
};
}  // namespace mcac


#endif //INCLUDE_PHYSICAL_MODEL_FLAME_COUPLING_HPP
