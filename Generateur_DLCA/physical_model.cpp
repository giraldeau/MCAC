#include "physical_model.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

const double PI = atan(1.0)*4;

void PhysicalModel::Init(const double _P, const double _T, const double _dfe, const double _kfe, const double _Dpm, const double _sigmaDpm, const double _xsurfgrowth, const double _coeffB, const double _Rho)
{
    dfe=_dfe;
    kfe=_kfe;
    P = _P;
    T=_T;
    Dpm=_Dpm;
    sigmaDpm=_sigmaDpm;
    xsurfgrowth = _xsurfgrowth;
    coeffB = _coeffB;
    Rho = _Rho;

    K = 1.38066E-23;
    lambda = 66.5E-9*(101300/P)*(T/293.15)*(1+110/293.15)/(1+110/T);
    Dpeqmass = Dpm*exp(1.5*log(sigmaDpm)*log(sigmaDpm)); //Diamètre équivalent massique moyen des monomères
                                                         //donné par l'équation de Hatch-Choate
    rpeqmass = (Dpeqmass*1E-9)/2.0; //Rayon équivalent massique moyen des monomères
    gamma_ = 1.378*(0.5+0.5*myerf(((lambda/rpeqmass)+4.454)/10.628));
    Mu = 18.203E-6*(293.15+110)/(T+110)*pow(T/293.15,1.5);

    Asurfgrowth = coeffB*1E-3;

    use_verlet = true; // Bool used to chose if the script will run a Verlet list, significantly reducing the cost of Calcul_Distance
    GridDiv = 10;      // Number of Divisions of the box


    SetPrecision(1e-4);
    UseSecante();

    print();
}

void PhysicalModel::SetPrecision(const double _precision)
{
    precision = _precision;
}

void PhysicalModel::UseDichotomia(void)
{
    root_method = 0;
}

void PhysicalModel::UseBrent(void)
{
    root_method = 1;
}
void PhysicalModel::UseSecante(void)
{
    root_method = 2;
}

void PhysicalModel::print(void) const
{

    string RootMethod;
    if (root_method==0)
        RootMethod = "Dichotomia";
    else if (root_method==1)
        RootMethod = "Brent";
    else if (root_method==2)
        RootMethod = "Secante";


    cout << "Physical parameters:" << endl
         << " Pressure    : " << P << endl
         << " Temperature : " << T << endl
         << " diffusivity : " << Mu << endl
         << " K           : " << K << endl
         << " density     : " << Rho << endl
         << " Dpm         : " << Dpm << endl
         << " sigmaDpm    : " << sigmaDpm << endl
         << " Asurfgrowth : " << Asurfgrowth << endl
         << " xsurfgrowth : " << xsurfgrowth << endl
         << " coeffB      : " << coeffB << endl
         << " dfe         : " << dfe << endl
         << " kfe         : " << kfe << endl
         << " lambda      : " << lambda << endl
         << " Dpeqmass    : " << Dpeqmass << endl
         << " rpeqmass    : " << rpeqmass << endl
         << " gamma_      : " << gamma_ << endl
         << endl
         << "Options for Pysical model: " << endl
         << " precision   : " << precision << endl
         << " root method : " << RootMethod << endl;
}


//###############################################################################################################################

double PhysicalModel::Cunningham(const double R) const //Facteur correctif de Cunningham
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

double PhysicalModel::ModeleBeta(const double rm,const double np, const double rg) const
{
    if (rm<0.) {throw 42;}
    return Cunningham(rm) - FactorModelBeta*rm;
}

double PhysicalModel::Dichotomie(const double np, const double rg, const double rpmoy, const double x0) const
{
    double 	rmin, rmax, rmed, frmed, frmin, frmax;
    int ite=1;
    rmin = x0/100;   //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    rmax = 2*x0; //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm

    frmin = ModeleBeta(rmin, np, rg);
    frmax = ModeleBeta(rmax, np, rg);

    if (frmin*frmax>=0) {printf("Intervalle incorrect : %e %e \n",frmin,frmax); return -1;} //Intervalle incorrect

    rmed = (rmin+rmax)/2 ;
    frmed = ModeleBeta(rmed, np, rg);

    while (fabs(frmed)>precision)
    {

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
        frmed = ModeleBeta(rmed, np, rg);
        ite++;
    }
    return rmed;
}

double PhysicalModel::brentq(const double np, const double rg, const double rpmoy, const double x0) const
{
    double xa,xb,xtol,rtol;
    int iter=500;

    xa = x0/100;   //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    xb = 2*x0; //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm
    xtol = rtol = precision;

    double xpre = xa, xcur = xb;
    double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0., sbis;
    /* the tolerance is 2*delta */
    double delta;
    double stry, dpre, dblk;
    int i;

    fpre = ModeleBeta(xpre, np, rg);
    fcur = ModeleBeta(xcur, np, rg);
    if (fpre*fcur > 0) {printf("Intervalle incorrect : %e %e \n",fpre,fcur); return -1;} //Intervalle incorrect

    if (fpre == 0) {
        return xpre;
    }
    if (fcur == 0) {
        return xcur;
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

        delta = (xtol + rtol*fabs(xcur))/2;
        sbis = (xblk - xcur)/2;
        //if (fcur == 0 || fabs(sbis) < delta) {
        if (fabs(fcur) < precision) {
            return xcur;
        }

        if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
            if (xpre == xblk) {
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

        fcur = ModeleBeta(xcur, np, rg);
    }
    return xcur;
}

double PhysicalModel::Secante(const double np,const double rg,const double rpmoy, const double xsave) const
{
    double x[2];
    double fx[2];
    bool i=1;

    x[0]=xsave+4.e-13;
    fx[0] = ModeleBeta(x[0], np, rg);

    x[1]=xsave;

    try
    {
        while ((fabs(fx[!i])>precision))
        {
            fx[i]=ModeleBeta(x[i],np,rg);
            x[!i]=x[i]-fx[i]*(x[i]-x[!i])/(fx[i]-fx[!i]);
            i = !i;
        }
    }
    catch (int e)
    {
        printf("FAILSAFE\n");
        return brentq(np,rg,rpmoy,xsave);
    }
    return x[i];
}

double PhysicalModel::ConvertRg2Dm(const double np, const double rg,const double rpmoy)
{
    FactorModelBeta = pow(kfe,-1/dfe)*pow(np,(1-gamma_)/dfe)*Cunningham(rpeqmass)/rg;
    double x0 = rpmoy*pow(np,gamma_/dfe);
/*
    npeqmass = kfe*pow(rg/rpeqmass,dfe);
    double nptmp = pow(npeqmass*pow(np,gamma_-1),1./gamma_);
    FactorModelBeta = 1./(rpeqmass*pow(nptmp,gamma_/dfe)/Cunningham(rpeqmass));
*/
    if (root_method==0)
        return Dichotomie(np,rg,rpmoy,x0)*2;
    else if (root_method==1)
        return brentq(np,rg,rpmoy,x0)*2;
    else if (root_method==2)
        return Secante(np,rg,rpmoy,x0)*2;

    return  Dichotomie(np,rg,rpmoy,x0)*2;
}

double PhysicalModel::Grow(const double R,const double dt) const
{
    return R + Asurfgrowth*pow(R, xsurfgrowth-2)*dt;
}

double PhysicalModel::diffusivity(const double dm) const
{
    double cc = Cunningham(dm/2);
    return K*T/3/PI/Mu/dm*cc;
}
double PhysicalModel::velocity(const double masse) const
{
    return sqrt(8*K*T/PI/masse);
}

//###################################################### Fonction Erf ###########################################################
const int ncof=28;

const double cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
        1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
        3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
        -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
        6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
        9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
        -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};

double erfccheb(const double z)
{
        int j;
        double t,ty,tmp,d=0.,dd=0.;
        //if (z < 0.) throw("erfccheb requires nonnegative argument");
        t = 2./(2.+z);
        ty = 4.*t - 2.;
        for (j=ncof-1;j>0;j--) {
                tmp = d;
                d = ty*d - dd + cof[j];
                dd = tmp;
        }
        return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}
double myerfc(const double x)
{
                if (x >= 0.) return erfccheb(x);
                else return 2.0 - erfccheb(-x);
}

double inverfc(const double p)
{
        double x,err,t,pp;
        if (p >= 2.0) return -100.;
        if (p <= 0.0) return 100.;
        pp = (p < 1.0)? p : 2. - p;
        t = sqrt(-2.*log(pp/2.));
        x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
        for (int j=0;j<2;j++) {
                err = myerfc(x) - pp;
                x += err/(1.12837916709551257*exp(-(x*x))-x*err);
        }
        return (p < 1.0? x : -x);
}
double myerf(const double x) { return 1-myerfc(x); }
double inverf(const double p) {return inverfc(1.-p);}
