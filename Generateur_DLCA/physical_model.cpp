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

using namespace std;

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
    gamma_ = 1.378*(0.5+0.5*erf(((lambda/rpeqmass)+4.454)/10.628));
    Mu = 18.203E-6*(293.15+110)/(T+110)*pow(T/293.15,1.5);

    Asurfgrowth = coeffB*1E-3;

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
    root_dicho = true;
    root_brent = false;
    root_sec = false;
}

void PhysicalModel::UseBrent(void)
{
    root_dicho = false;
    root_brent = true;
    root_sec = false;
}
void PhysicalModel::UseSecante(void)
{
    root_dicho = false;
    root_brent = false;
    root_sec = true;
}

void PhysicalModel::print(void) const
{
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
         << " root_dicho  : " << root_dicho << endl
         << " root_sec    : " << root_sec << endl
         << " root_brent  : " << root_brent << endl;
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
    if (root_dicho) return Dichotomie(np,rg,rpmoy,x0)*2;
    if (root_sec) return Secante(np,rg,rpmoy,x0)*2;

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
