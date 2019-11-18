#include "calcul.hpp"
#include "physical_model.hpp"
#include <cstdio>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    cout << "    MCAC  Copyright (C) 2019 CORIA" << endl;

    if(argc <=1)
    {
        cout << "Missing argument : param file." << endl;
        return 1;
    }

    string FichierParam = argv[1];

    DLCA::PhysicalModel physicalmodel(FichierParam);
    DLCA::Calcul(physicalmodel);

    return 0;
}
