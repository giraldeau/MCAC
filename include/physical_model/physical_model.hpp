#ifndef PHYSICAL_MODEL_H
#define PHYSICAL_MODEL_H 1
#include <cst.hpp>
#include <io/xmf_includes.hpp>
#include <ctime>
#include <array>
#include <experimental/filesystem>


namespace MCAC {
class PhysicalModel {
public:
    // Physics
    double dfe, kfe;                            // dimension fractale et préfacteur fractal
    double xsurfgrowth, coeffB;                 // surface growth parameter, Bêta
    double Asurfgrowth;                         // surface growth coefficient
    double P, T, lambda, Mu, K, Rho;                    // pressure, temperature, diffusivity, ?, density
    double Dpm, sigmaDpm;
    double Dpeqmass, rpeqmass, gamma_;  // libre parcours moyen d'une sphère
    double Time;                                // Time
    double X, FV, L;                            // size of the box, Volume ratio, lenght of the box, pressure, density
    double precision;
    double FactorModelBeta;
    clock_t CPUStart;
    double CPULimit, PHY_Time_limit;
    size_t NPP_avg_limit;
    size_t GridDiv;     // Number of Divisions of the box
    size_t N;           // Nombre de sphères initial
    size_t AggMin;
    size_t DeltaSauve;
    int root_method;
    int Mode;
    int Wait, WaitLimit;
    bool ActiveModulephysique;
    bool use_verlet;    // Bool used to chose if the script will run a Verlet list

    std::experimental::filesystem::path CheminSauve;
private:
    // Internal use
    double ModeleBeta(double rm) const;
    double Dichotomie(double x0) const;
    double brentq(double x0) const;
    double Secante(double x0) const;
public:
    explicit PhysicalModel(const std::string &FichierParam);
    void SetPrecision(double _precision);
    void UseDichotomia();
    void UseBrent();
    void UseSecante();
    double ConvertRg2Dm(size_t np, double rg, double rmoy);
    double ConvertRg2DmFromStart(size_t np, double rg, double start);
    double Cunningham(double R) const;
    double Grow(double R, double dt) const;
    double friction_coeff(size_t npp) const;
    double friction_coeff2(double rgg) const;
    double diffusivity(double) const;
    double velocity(double masse) const;
    double relax_time(double masse, double) const;
    void print() const;
    bool Finished(size_t Nagg, double NPP_avg) const;
    XMF_OUTPUT xmf_write() const;
};

double erfccheb(double z);
double inverfc(double p);
double inverf(double p);
template<typename T>
__attribute((const)) inline double periodicDistance(T x, T dim) noexcept {
    //return periodicPosition(x+0.5*dim,dim)-0.5*dim;
    T dist(x);
    T hdim(0.5 * dim);
    while (dist < -hdim) {
        dist += dim;
    }
    while (dist >= hdim) {
        dist -= dim;
    }
    return dist;
}
__attribute((const)) inline std::array<double, 3> periodicDistance(std::array<double, 3> x, double dim) noexcept {
    return {periodicDistance(x[0], dim),
            periodicDistance(x[1], dim),
            periodicDistance(x[2], dim)};
}
__attribute((const)) inline double periodicPosition(double x, double dim) noexcept {
    //return fmod(dim-fmod(dim-x,dim),dim);
    //return fmod(x+dim,dim);

    double pos(x);
    while (pos < 0) {
        pos += dim;
    }
    while (pos >= dim) {
        pos -= dim;
    }
    return pos;
}
__attribute((const)) inline std::array<double, 3> periodicPosition(std::array<double, 3> x, double dim) noexcept {
    return {periodicPosition(x[0], dim), periodicPosition(x[1], dim), periodicPosition(x[2], dim)};
}
bool locale_with_dots();
double latof(const char *_char);
bool dirExists(const char *path);
std::experimental::filesystem::path extractPath(const std::string &file);
std::tuple<bool, double, double, double> linreg(const std::vector<double> &x, const std::vector<double> &y);
}  // namespace MCAC


#endif // PHYSICAL_MODEL_H

