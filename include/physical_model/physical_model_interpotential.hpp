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
#ifndef INCLUDE_PHYSICAL_MODEL_INTERPOTENTIAL_HPP
#define INCLUDE_PHYSICAL_MODEL_INTERPOTENTIAL_HPP 1
#include "constants.hpp"
#include <vector>
#include <ctime>


namespace mcac {
/********************************************************************************
* Class for saving the info related to interaction potentials
********************************************************************************/
class Interpotential {
private:
    int max_charge,min_charge; ///< bounds of the vector charges
    double fixed_Temperature;  ///< Temperature normalized potentials [K]
    size_t N_val_p1, N_val_p2, N_charge;
    std::vector<int> val_charge;        ///< possible charges
    std::vector<double> val_dp1;        ///< diameters particle 1 [m]
    std::vector<double> val_dp2;        ///< diameters particle 2 [m]
    std::vector<std::vector<std::vector<std::vector<double>>>> E_barr;         ///< Energy barrier [dimensionless, factor of k_b*T]
    std::vector<std::vector<std::vector<std::vector<double>>>> E_well;         ///< Energy potential well [dimensionless, factor of k_b*T]
public:
    /// Constructor
    Interpotential();
    explicit Interpotential(const std::string &interpotential_file);
    /* getters */
    double get_fixed_Temperature() const noexcept;
    int get_max_charge() const noexcept;
    int get_min_charge() const noexcept;
    std::pair<double, double> get_Ebar_Ewell(const double dp1, const double dp2, const int charge1, const int charge2);
private:
    void set_min_max_charge();
};
}  // namespace mcac


#endif //INCLUDE_PHYSICAL_MODEL_INTERPOTENTIAL_HPP
