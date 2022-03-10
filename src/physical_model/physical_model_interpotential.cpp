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
#include "constants.hpp"
#include "physical_model/physical_model_interpotential.hpp"
#include "exceptions.hpp"
#include "tools/tools.hpp"
#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <algorithm>    // std::find


namespace fs = std::experimental::filesystem;
namespace mcac {
Interpotential::Interpotential() :
    max_charge(6),
    min_charge(-6),
    fixed_Temperature(),
    N_val_p1(),
    N_val_p2(),
    N_charge(),
    val_charge(),
    val_dp1(),
    val_dp2(),
    E_barr(),
    E_well() {
}
/********************************************************************************
* import the information to simulate inter-potentials
* Filling the class Interpotential with this information
********************************************************************************/
Interpotential::Interpotential(const std::string &interpotential_file) : Interpotential() {
    std::cout << "LOADING THE Interpotential file" << std::endl;
    if (!fs::exists(interpotential_file)) {
        throw IOError(" Interpotential file does not exist: " + interpotential_file);
    }
    std::ifstream ext_FILE;
    ext_FILE.open(interpotential_file, std::ifstream::in);
    double d_temp; // just temporal auxiliary variables
    //WE ASSUME: External DATA is ordered as:
    // 1st line: Fixed temperature
    ext_FILE >> fixed_Temperature;

    // 2nd line: N_val_p1
    ext_FILE >> d_temp;
    N_val_p1 = static_cast<size_t>(d_temp);
    val_dp1.resize(N_val_p1);

    // 3rd line: N_val_p2
    ext_FILE >> d_temp;
    N_val_p2 = static_cast<size_t>(d_temp);
    val_dp2.resize(N_val_p2);

    // 4th line: N_charge
    ext_FILE >> d_temp;
    N_charge = static_cast<size_t>(d_temp);
    val_charge.resize(N_charge);

    //    size_t nlines = N_charge * N_charge * N_val_p1 * N_val_p2;
    E_barr.resize(N_charge);
    E_well.resize(N_charge);
    for (size_t i=0; i<N_charge; i++){
        E_barr[i].resize(N_charge);
        E_well[i].resize(N_charge);
        for (size_t j=0; j<N_charge; j++){
            E_barr[i][j].resize(N_val_p1);
            E_well[i][j].resize(N_val_p1);
            for (size_t k=0; k<N_val_p1; k++){
                E_barr[i][j][k].resize(N_val_p2);
                E_well[i][j][k].resize(N_val_p2);
            }
        }
    }

    // 5th line: val_charge
    for (size_t i=0; i<N_charge; i++) {
        ext_FILE >> d_temp;
        val_charge[i] = static_cast<int>(d_temp);
    }
    // 6th line: val_dp1
    for (size_t i=0; i<N_val_p1; i++) {
        ext_FILE >> val_dp1[i];
    }
    // 7th line: val_dp2
    for (size_t i=0; i<N_val_p2; i++) {
        ext_FILE >> val_dp2[i];
    }
    // >=8th line: [E_barr, E_well]
//    for (int i=0; i<nlines; i++) {
//        ext_FILE >> E_barr[i] >> E_well[i];
//    }
    for (size_t i=0; i<N_val_p1; i++){
        for (size_t j=0; j<N_val_p2; j++){
            for (size_t k=0; k<N_charge; k++){
                for (size_t l=0; l<N_charge; l++){
                    ext_FILE >> E_barr[k][l][i][j] >> E_well[k][l][i][j]; // TODO change the input file
                }
            }
        }
    }
    ext_FILE.close();

    // Set min and max charges
    set_min_max_charge();
}
/********************************************************************************
* Access information
********************************************************************************/
/* getters */
[[gnu::pure]] double Interpotential::get_fixed_Temperature() const noexcept {
    return fixed_Temperature;
}
[[gnu::pure]] int Interpotential::get_max_charge() const noexcept {
    return max_charge;
}
[[gnu::pure]] int Interpotential::get_min_charge() const noexcept {
    return min_charge;
}
/********************************************************************************
* Set min and maximum charge
********************************************************************************/
void Interpotential::set_min_max_charge() {
    min_charge = *std::min_element(val_charge.begin(), val_charge.end());
    max_charge = *std::max_element(val_charge.begin(), val_charge.end());
}
/********************************************************************************
* Update E_bar and E_well
********************************************************************************/
std::pair<double, double> Interpotential::get_Ebar_Ewell(const double dp1, const double dp2, const int charge1, const int charge2) {
    // Search dp1
    auto next_dp1 = std::upper_bound(val_dp1.begin(), val_dp1.end(), dp1);
    auto previous_dp1 = next_dp1 - 1;
    if (next_dp1 == val_dp1.begin() || next_dp1 == val_dp1.end()) {
        throw InterPotentialError("PP diameter out of range");
    }
    // Search dp2
    auto next_dp2 = std::upper_bound(val_dp2.begin(), val_dp2.end(), dp2);
    auto previous_dp2 = next_dp2 - 1;
    if (next_dp2 == val_dp2.begin() || next_dp2 == val_dp2.end()) {
        throw InterPotentialError("PP diameter out of range");
    }
    // Search charge 1
    auto charg1 = std::find(val_charge.begin(), val_charge.end(), charge1);
    if (charg1 == val_charge.end()) {
        throw InterPotentialError("charges out of range");
    }
    // Search charge 2
    auto charg2 = std::find(val_charge.begin(), val_charge.end(), charge2);
    if (charg2 == val_charge.end()) {
        throw InterPotentialError("charges out of range");
    }

    auto loc_previous_dp1 = static_cast<size_t>(std::distance( val_dp1.begin(), previous_dp1 ));
    auto loc_next_dp1 = static_cast<size_t>(std::distance( val_dp1.begin(), next_dp1 ));

    auto loc_previous_dp2 = static_cast<size_t>(std::distance( val_dp2.begin(), previous_dp2 ));
    auto loc_next_dp2 = static_cast<size_t>(std::distance( val_dp2.begin(), next_dp2 ));

    auto loc_charg1 = static_cast<size_t>(std::distance(val_charge.begin(), charg1 ));
    auto loc_charg2 = static_cast<size_t>(std::distance(val_charge.begin(), charg2 ));

    double delta_dp1 = *next_dp1 - *previous_dp1;
    double delta_dp2 = *next_dp2 - *previous_dp2;
    double d_dp1_D_dp1 = (dp1- *previous_dp1)/delta_dp1;
    double d_dp2_D_dp2 = (dp2- *previous_dp2)/delta_dp2;

    // Interpolate E_bar
    double E_bar_prev_prev = E_barr[loc_charg1][loc_charg2][loc_previous_dp1][loc_previous_dp2];
    double E_bar_prev_next = E_barr[loc_charg1][loc_charg2][loc_previous_dp1][loc_next_dp2];
    double E_bar_next_prev = E_barr[loc_charg1][loc_charg2][loc_next_dp1][loc_previous_dp2];
    double E_bar_next_next = E_barr[loc_charg1][loc_charg2][loc_next_dp1][loc_next_dp2];
    double interpolated_E_bar = interpolate_2d(E_bar_prev_prev,E_bar_prev_next,E_bar_next_prev,E_bar_next_next,d_dp1_D_dp1,d_dp2_D_dp2);

    // Interpolate E_well
    double E_well_prev_prev = E_well[loc_charg1][loc_charg2][loc_previous_dp1][loc_previous_dp2];
    double E_well_prev_next = E_well[loc_charg1][loc_charg2][loc_previous_dp1][loc_next_dp2];
    double E_well_next_prev = E_well[loc_charg1][loc_charg2][loc_next_dp1][loc_previous_dp2];
    double E_well_next_next = E_well[loc_charg1][loc_charg2][loc_next_dp1][loc_next_dp2];
    double interpolated_E_well = interpolate_2d(E_well_prev_prev,E_well_prev_next,E_well_next_prev,E_well_next_next,d_dp1_D_dp1,d_dp2_D_dp2);

    return std::make_pair(interpolated_E_bar, interpolated_E_well);
}
}  // namespace mcac
