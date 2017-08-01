#ifndef PHYSICAL_MODEL
#define PHYSICAL_MODEL

#include <iostream>

using namespace std;

class PhysicalModel
{
public:

    // Physics
    double Asurfgrowth; // surface growth coefficient
    double dfe, kfe; // dimension fractale et préfacteur fractal
    double xsurfgrowth, coeffB; // surface growth parameter, Bêta
    double lambda, Dpeqmass, rpeqmass, gamma_; // libre parcours moyen d'une sphère
    double P, T, Mu, K, Rho; // pressure, temperature, diffusivity, ?, density
    double Dpm, sigmaDpm;
    double temps; // Time
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

    int step;

private:

    // Internal use
    double ModeleBeta(const double rm) const;
    double Dichotomie(const double x0) const;
    double brentq(const double x0) const;
    double Secante(const double x0) const;

public:
    void Init(const double _P, const double _T, const double _dfe, const double _kfe, const double _Dpm, const double _sigmaDpm, const double _xsurfgrowth, const double _coeffB, const double _Rho);
    void SetPrecision(const double _precision);
    void UseDichotomia(void);
    void UseBrent(void);
    void UseSecante(void);
    double ConvertRg2Dm(const double np, const double rg, const double rmoy);
    double ConvertRg2DmFromStart(const double np, const double rg, const double start);

    double Cunningham(const double R) const;
    double Grow(const double R,const double dt) const;
    double diffusivity(const double dm) const;
    double velocity(const double masse) const;
    void print(void) const;
};


double erfccheb(const double z);
double myerfc(const double x);
double inverfc(const double p);
double myerf(const double x);
double inverf(const double p);

#endif // PHYSICAL_MODEL

