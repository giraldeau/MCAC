#include "mainwindow.h"

#include <boost/python.hpp>


using namespace boost::python;

BOOST_PYTHON_MODULE(libDLCA)
{
    // Add regular functions to the module.
    def("Calcul", DLCA::Calcul);
    def("LectureParams",DLCA::LectureParams);
    class_<DLCA::PhysicalModel>("PhysicalModel")
//            .add_property("X",)
            ;

    /*
    double Asurfgrowth; // surface growth coefficient
    double dfe, kfe; // dimension fractale et préfacteur fractal
    double xsurfgrowth, coeffB; // surface growth parameter, Bêta
    double lambda, Dpeqmass, rpeqmass, gamma_; // libre parcours moyen d'une sphère
    double P, T, Mu, K, Rho; // pressure, temperature, diffusivity, ?, density
    double Dpm, sigmaDpm;
    double time; // Time
    double X, FV, L; // Temperature, size parameter of the box, Volume ratio, lenght of the box, pressure, density
    double precision;
    double FactorModelBeta;

    int N;// Nombre de sphères initial, bool pour l'activation du module phy, bool pour l'activation de la variation de temps
    int DeltaSauve;
    int root_method;
    int Mode;
    int GridDiv; // Number of Divisions of the box
    bool ActiveModulephysique, ActiveVariationTempo;
    bool use_verlet; // Bool used to chose if the script will run a Verlet list, significantly reducing the cost of Calcul_Distance
    */
}
