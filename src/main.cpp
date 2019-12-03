#include "calcul.hpp"
#include "physical_model.hpp"
#include "tools.hpp"
#include "cst.hpp"
#include <iostream>


int main(int argc, char *argv[]) {
    std::cout << "    MCAC  Copyright (C) 2019 CORIA" << std::endl;
    if (argc <= 1) {
        std::cout << "Missing argument : param file." << std::endl;
        return 1;
    }
    MCAC::PhysicalModel physicalmodel(argv[1]);
    MCAC::InitRandom();
    MCAC::ListAggregat aggregates(&physicalmodel);
    MCAC::Calcul(physicalmodel, aggregates);
    return MCAC::ErrorCodes::VERLET_ERROR;
}
