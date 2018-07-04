#include "physical_model.hpp"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))


using namespace std;

namespace DLCA{

const double PI = atan(1.0)*4;

PhysicalModel::PhysicalModel() :
    Asurfgrowth(0),
    dfe(0),
    kfe(0),
    xsurfgrowth(0),
    coeffB(0),
    lambda(0),
    Dpeqmass(0),
    rpeqmass(0),
    gamma_(0),
    P(0),
    T(0),
    Mu(0),
    K(0),
    Rho(0),
    Dpm(0),
    sigmaDpm(0),
    Time(0),
    X(0),
    FV(0),
    L(DBL_MAX),
    precision(0),
    FactorModelBeta(0),
    CPUStart(0),
    CPULimit(0),
    GridDiv(0),
    N(0),
    AggMin(0),
    DeltaSauve(0),
    root_method(0),
    Mode(0),
    Wait(0),
    WaitLimit(0),
    ActiveModulephysique(false),
    ActiveVariationTempo(false),
    use_verlet(false),
    toBeDestroyed(true),
    CheminSauve()
{}




void PhysicalModel::Init()
{

    L = X*Dpm*1E-9;
    Time=0;

    K = 1.38066E-23;
    lambda = 66.5E-9*(101300/P)*(T/293.15)*(1+110/293.15)/(1+110/T);
    Dpeqmass = Dpm*exp(1.5*log(sigmaDpm)*log(sigmaDpm)); //Diamètre équivalent massique moyen des monomères
                                                         //donné par l'équation de Hatch-Choate
    rpeqmass = (Dpeqmass*1E-9)/2.0; //Rayon équivalent massique moyen des monomères
    gamma_ = 1.378*(0.5+0.5*erf(((lambda/rpeqmass)+4.454)/10.628));
    Mu = 18.203E-6*(293.15+110)/(T+110)*pow(T/293.15,1.5);

    Asurfgrowth = coeffB*1E-3;

    use_verlet = true; // Bool used to chose if the script will run a Verlet list, significantly reducing the cost of Calcul_Distance
    GridDiv = 10;      // Number of Divisions of the box

    AggMin = MAX(AggMin, 1);
    time(&CPUStart);
    SetPrecision(1e-5);
    UseSecante();

    print();

    toBeDestroyed = false;
}

bool PhysicalModel::Finished(const size_t Nagg) const
{
    if (Nagg <= AggMin)
    {
        cout << "We reach the AggMin condition" << endl << endl;
        return true;
    }
    if (WaitLimit > 0 && Wait >= WaitLimit)
    {
        cout << "We reach the WaitLimit condition" << endl << endl;
        return true;
    }

    if (CPULimit > 0)
    {
        time_t currentCPU;
        time(&currentCPU);
        if (currentCPU - CPUStart >= CPULimit)
        {
            cout << "We reach the CPULimit condition" << endl << endl;
            return true;
        }
    }

    return false;
}

void PhysicalModel::SetPrecision(const double _precision)
{
    precision = _precision;
}

void PhysicalModel::UseDichotomia()
{
    root_method = 0;
}

void PhysicalModel::UseBrent()
{
    root_method = 1;
}
void PhysicalModel::UseSecante()
{
    root_method = 2;
}

void PhysicalModel::print() const
{

    string RootMethod;
    switch(root_method)
    {
    case 0 :
        RootMethod = "Dichotomia";
        break;
    case 1:
        RootMethod = "Brent";
        break;
    default:
        cout << "Root method unknown, using secante" << endl;
       /* FALLTHRU */
    case 2 :
        RootMethod = "Secante";
    }


    cout << "Physical parameters:" << endl
         << " Initial Nagg : " << N << endl
         << " Box size     : " << L << endl
         << " X            : " << X << endl
         << " Pressure     : " << P << endl
         << " Temperature  : " << T << endl
         << " diffusivity  : " << Mu << endl
         << " K            : " << K << endl
         << " FV           : " << FV << endl
         << " density      : " << Rho << endl
         << " Dpm          : " << Dpm << endl
         << " sigmaDpm     : " << sigmaDpm << endl
         << " Asurfgrowth  : " << Asurfgrowth << endl
         << " xsurfgrowth  : " << xsurfgrowth << endl
         << " coeffB       : " << coeffB << endl
         << " dfe          : " << dfe << endl
         << " kfe          : " << kfe << endl
         << " lambda       : " << lambda << endl
         << " Dpeqmass     : " << Dpeqmass << endl
         << " rpeqmass     : " << rpeqmass << endl
         << " gamma_       : " << gamma_ << endl
         << endl
         << "Options for Pysical model: " << endl
         << " precision   : " << precision << endl
         << " root method : " << RootMethod << endl
         << " Mode : " << Mode << endl
         << endl
         << "Ending calcul when:" << endl
         << " - There is " << AggMin << " aggregats left or less" << endl;
    if (WaitLimit > 0)
    {
        cout << " - It has been " << WaitLimit << " iterations without collision" <<endl;
    }
    if (CPULimit > 0)
    {
        cout << " - The simulations is running for more than " << CPULimit << " seconds" << endl;
    }
    cout << endl;
}


//###############################################################################################################################

 __attribute__((pure)) double PhysicalModel::Cunningham(const double R) const //Facteur correctif de Cunningham
{
    double A = 1.142;
    double B = 0.558;
    double C = 0.999;
    return 1+A*lambda/R+B*lambda/R*exp(-C*R/lambda);
}


//######################################## Fonctions pour le calcul du diamètre de mobilité #####################################

/*
 Fonction permettant de retrouver le rayon de mobilité en régime transitoire
 On obtient le bon rayon de mobilité lorsque la fonction retourne 0
*/

double PhysicalModel::ModeleBeta(const double rm) const
{
    if (rm<0.) {throw 42;}
    return Cunningham(rm) - FactorModelBeta*rm;
}

double PhysicalModel::Dichotomie(const double x0) const
{
    double 	rmin, rmax, rmed, frmed, frmin, frmax;
    rmin = x0/100;   //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    rmax = 2*x0; //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm

    frmin = ModeleBeta(rmin);
    frmax = ModeleBeta(rmax);

    rmed = (rmin+rmax)/2 ;
    frmed = ModeleBeta(rmed);

    while (fabs(frmed)>precision)
    {
        if (frmin*frmax>=0)
        {
            cout << "Intervalle incorrect : " << frmin << " " << frmax << endl;
            return -1;
        }
        if (frmin*frmed < 0)
        {
            rmax = rmed;
            frmax = frmed;
        }
        else
        {
            rmin = rmed;
            frmin = frmed;
        }
        rmed = (rmin+rmax)/2 ;
        frmed = ModeleBeta(rmed);
    }
    return rmed;
}

double PhysicalModel::brentq(const double x0) const
{
    double xa,xb,xtol,rtol;
    int iter=500;

    xa = x0/100;   //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    xb = 2*x0; //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm
    xtol = rtol = precision;

    double xpre = xa, xcur = xb;
    double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0.;
    double stry, dpre, dblk;
    int i;

    fpre = ModeleBeta(xpre);
    if (fabs(fpre) < precision) {
        return xpre;
    }
    fcur = ModeleBeta(xcur);
    if (fabs(fcur) < precision) {
        return xcur;
    }

    if (fpre*fcur > 0)
    {
        cout << "Intervalle incorrect : " << fpre << " " << fcur << endl;
        return -1;
    }

    for (i = 0; i < iter; i++) {
        if (fpre*fcur < 0) {
            xblk = xpre;
            fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (fabs(fblk) < fabs(fcur)) {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        double delta = (xtol + rtol*fabs(xcur))/2;
        double sbis = (xblk - xcur)/2;
        //if (fcur == 0 || fabs(sbis) < delta) {
        if (fabs(fcur) < precision) {
            return xcur;
        }

        if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
            if (abs(xpre - xblk) <= 1e-16) {
                /* interpolate */
                stry = -fcur*(xcur - xpre)/(fcur - fpre);
            }
            else {
                /* extrapolate */
                dpre = (fpre - fcur)/(xpre - xcur);
                dblk = (fblk - fcur)/(xblk - xcur);
                stry = -fcur*(fblk*dblk - fpre*dpre)
                    /(dblk*dpre*(fblk - fpre));
            }
            if (2*fabs(stry) < MIN(fabs(spre), 3*fabs(sbis) - delta)) {
                /* good short step */
                spre = scur;
                scur = stry;
            } else {
                /* bisect */
                spre = sbis;
                scur = sbis;
            }
        }
        else {
            /* bisect */
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur; fpre = fcur;
        if (fabs(scur) > delta) {
            xcur += scur;
        }
        else {
            xcur += (sbis > 0 ? delta : -delta);
        }

        fcur = ModeleBeta(xcur);
    }
    return xcur;
}

double PhysicalModel::Secante(const double x0) const
{
    double x[2];
    double fx[2];
    bool i=1;

    x[0]=x0+4.e-13;
    fx[0] = ModeleBeta(x[0]);

    x[1]=x0;

    try
    {
        while ((fabs(fx[!i])>precision))
        {
            fx[i]=ModeleBeta(x[i]);
            x[!i]=x[i]-fx[i]*(x[i]-x[!i])/(fx[i]-fx[!i]);
            i = !i;
        }
    }
    catch (int)
    {
        return brentq(x0);
    }
    return x[i];
}

__attribute((pure)) double PhysicalModel::ConvertRg2Dm(const size_t np, const double rg,const double rmoy)
{
    double start = 2*rmoy*pow(np,gamma_/dfe);
    return ConvertRg2DmFromStart(np,rg,start);
}

double PhysicalModel::ConvertRg2DmFromStart(const size_t np, const double rg, const double start)
{
    FactorModelBeta = pow(kfe,-1/dfe)*pow(np,(1-gamma_)/dfe)*Cunningham(rpeqmass)/rg;
/*
    npeqmass = kfe*pow(rg/rpeqmass,dfe);
    double nptmp = pow(npeqmass*pow(np,gamma_-1),1./gamma_);
    FactorModelBeta = 1./(rpeqmass*pow(nptmp,gamma_/dfe)/Cunningham(rpeqmass));
*/
    switch (root_method)
    {
    case 2:
        return Secante(start*0.5)*2;
    case 1:
        return brentq(start*0.5)*2;
    default:
        return Dichotomie(start*0.5)*2;
    }
}

 __attribute__((pure)) double PhysicalModel::Grow(const double R,const double dt) const
{
    return R + Asurfgrowth*pow(R, xsurfgrowth-2)*dt;
}

__attribute((pure)) double PhysicalModel::diffusivity(const double dm) const
{
    double cc = Cunningham(dm/2);
    return K*T/3/PI/Mu/dm*cc;
}
 __attribute__((pure)) double PhysicalModel::velocity(const double masse) const
{
    return sqrt(8*K*T/PI/masse);
}

__attribute((const)) double inverfc(const double p)
{
        double x,t,pp;
        if (p >= 2.0){ return -100.;}
        if (p <= 0.0){ return 100.;}
        pp = (p < 1.0)? p : 2. - p;
        t = sqrt(-2.*log(pp/2.));
        x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
        for (int j=0;j<2;j++) {
                double err = erfc(x) - pp;
                x += err/(1.12837916709551257*exp(-(x*x))-x*err);
        }
        return (p < 1.0? x : -x);
}
__attribute((const)) double inverf(const double p) {return inverfc(1.-p);}

}  // namespace DLCA

