#include "calcul.hpp"
#include "constants.hpp"
#include "physical_model/physical_model.hpp"
#include "tools/tools.hpp"
#include <iostream>


int main(int argc, char *argv[]) {
    std::cout << "    MCAC  Copyright (C) 2019 CORIA" << std::endl;
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
