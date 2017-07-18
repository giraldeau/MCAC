#ifndef PHYSICAL_MODEL
#define PHYSICAL_MODEL

class PhysicalModel
{
private:

    // Physical constants
    double Asurfgrowth; // surface growth coefficient
    double dfe, kfe; // dimension fractale et préfacteur fractal
    double xsurfgrowth, coeffB; // surface growth parameter, Bêta
    double lambda, Dpeqmass, rpeqmass, gamma_; // libre parcours moyen d'une sphère
    double P, T, Mu, K, Rho; // pressure, temperature, diffusivity, ?, density
    double Dpm, sigmaDpm;

    // Options
    double precision, root_dicho, root_sec, root_brent;

    // Internal use
    double FactorModelBeta;

    double ModeleBeta(const double rm, const double np, const double rg) const;
    double Dichotomie(const double np, const double rg,const double rpmoy,const double x0) const;
    double brentq(const double np, const double rg,const double rpmoy,const double x0) const;
    double Secante(const double np,const double rg,const double rpmoy, const double xsave) const;

public:
    void Init(const double _P, const double _T, const double _dfe, const double _kfe, const double _Dpm, const double _sigmaDpm, const double _xsurfgrowth, const double _coeffB, const double _Rho);
    void SetPrecision(const double _precision);
    void UseDichotomia(void);
    void UseBrent(void);
    void UseSecante(void);
    double ConvertRg2Dm(const double np, const double rg, const double rpmoy);
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
