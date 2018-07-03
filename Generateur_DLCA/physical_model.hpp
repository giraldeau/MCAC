#ifndef PHYSICAL_MODEL_H
#define PHYSICAL_MODEL_H

#include <cstddef>
#include <ctime>
#include <experimental/filesystem>

namespace DLCA{

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
    double Time; // Time
    double X, FV, L; // Temperature, size parameter of the box, Volume ratio, lenght of the box, pressure, density
    double precision;
    double FactorModelBeta;

    time_t CPUStart, CPULimit;

    size_t GridDiv; // Number of Divisions of the box
    size_t N;// Nombre de sphères initial, bool pour l'activation du module phy, bool pour l'activation de la variation de temps
    size_t AggMin;
    int DeltaSauve;
    int root_method;
    int Mode;
    int Wait, WaitLimit;
    bool ActiveModulephysique, ActiveVariationTempo;
    bool use_verlet; // Bool used to chose if the script will run a Verlet list, significantly reducing the cost of Calcul_Distance

    bool toBeDestroyed;
    std::experimental::filesystem::path CheminSauve;
private:

    // Internal use
    double ModeleBeta(double rm) const;
    double Dichotomie(double x0) const;
    double brentq(double x0) const;
    double Secante(double x0) const;

public:
    PhysicalModel();
    void Init();
    void SetPrecision(double _precision);
    void UseDichotomia();
    void UseBrent();
    void UseSecante();
    double ConvertRg2Dm(size_t np, double rg, double rmoy);
    double ConvertRg2DmFromStart(size_t np, double rg, double start);

    double Cunningham(double R) const;
    double Grow(double R,double dt) const;
    double diffusivity(double dm) const;
    double velocity(double masse) const;
    void print() const;
    bool Finished(const size_t Nagg) const;

    auto xmfWrite() const;

};


double erfccheb(double z);
double inverfc(double p);
double inverf(double p);


template<typename T>
__attribute((const)) inline double periodicDistance(T x, T dim) noexcept
{
    //return periodicPosition(x+0.5*dim,dim)-0.5*dim;
    T dist(x);
    T hdim(0.5*dim);

    while (dist < -hdim)
    {
        dist += dim;
    }
    while (dist >=  hdim)
    {
        dist -= dim;
    }

    return dist;
}

template<typename T>
__attribute((const)) inline T periodicPosition(T x, T dim) noexcept
{
    //return fmod(dim-fmod(dim-x,dim),dim);
    //return fmod(x+dim,dim);

    T pos(x);

    while (pos < 0)
    {
        pos += dim;
    }
    while (pos >= dim)
    {
        pos -= dim;
    }

    return pos;
}

}  // namespace DLCA


#endif // PHYSICAL_MODEL_H

