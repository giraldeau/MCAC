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
#include "physical_model/physical_model_flame_coupling.hpp"
#include "exceptions.hpp"
#include <fstream>
#include <iostream>
#include <experimental/filesystem>


namespace fs = std::experimental::filesystem;
namespace mcac {
FlameCoupling::FlameCoupling() :
    t_res(),
    Temp(),
    fv(),
    dp_gav(),
    dp_gSTD(),
    u_sg(),
    J_nucl() {
}
/********************************************************************************
* import the time resolved information of the flame
* Filling the class flame_coupling with this information
********************************************************************************/
FlameCoupling::FlameCoupling(const std::string &flame_file) : FlameCoupling() {
    std::cout << "LOADING THE FLAME INFORMATION (as a function of time)" << std::endl;
    std::cout << "be careful with the units in the imput file!" << std::endl;
    if (!fs::exists(flame_file)) {
        throw IOError("Flame file does not exist: " + flame_file);
    }
    std::ifstream ext_FILE;
    ext_FILE.open(flame_file, std::ifstream::in);
    //WE ASUME: External DATA CONTAINS 7 COLUMNS AS: [t_res,T,fv,dp_Gav,dp_gstd,u_sg,J_nucl]
    std::cout << "t_res, T, fv, dp_Gav, dp_gstd, u_sg, J_nucl" << std::endl;
    std::array<double, 7> line;
    while (ext_FILE >> line[0] >> line[1] >> line[2] >> line[3] >> line[4] >> line[5] >> line[6]) {
        // becareful with the units in the imput file! -> input-file should be in SI units
        t_res.push_back(line[0]);
        Temp.push_back(line[1]);
        fv.push_back(line[2]);
        dp_gav.push_back(line[3]);
        dp_gSTD.push_back(line[4]);
        u_sg.push_back(line[5]);
        J_nucl.push_back(line[6]);
        std::cout << line[0] << " " << line[1] << " " << line[2] << " " <<
                     line[3] << " " << line[4] << " " << line[5] << " " << line[6] << std::endl;
    }
    ext_FILE.close();
}
}  // namespace mcac

