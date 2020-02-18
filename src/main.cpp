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
#include "calcul.hpp"
#include "constants.hpp"
#include "physical_model/physical_model.hpp"
#include "tools/tools.hpp"
#include <iostream>


int main(int argc, char *argv[]) {
    std::cout << std::endl;
    std::cout << "MCAC  Copyright (C) 2020 CORIA" << std::endl;
    std::cout << "This program comes with ABSOLUTELY NO WARRANTY" << std::endl;
    std::cout << std::endl;
    if (argc <= 1) {
        std::cout << "Missing argument : param file." << std::endl;
        return 1;
    }
    mcac::PhysicalModel physicalmodel(argv[1]);
    mcac::init_random();
    mcac::AggregatList aggregates(&physicalmodel);
    mcac::calcul(&physicalmodel, &aggregates);
    return mcac::ErrorCodes::NO_ERROR;
}
