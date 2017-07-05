#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <Sphere.h>
#include <string.h>
#include <qdir.h>
#include <QString>
#include <stdio.h>
#include <string>
#include <cmath>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

QString FichierParam = "___";
QString pathParam = "___";
QString FichierSuiviTempo = "___";
QString pathSuiviTempo = "___";

int N, ActiveModulephysique, ActiveVariationTempo;// Nombre de sphères initial, bool pour l'activation du module phy, bool pour l'activation de la variation de temps
char CheminSauve[500];
char commentaires[500];

using namespace std;

const double PI = atan(1.0)*4; // 3.14159
double T, X, FV, L, Mu, K, P, Rho; // Ture, Paramètre de taille de la boîte, Fraction volumique de suie, longueur de la boîte, ?, ?, ?, Masse volumique
double Asurfgrowth; // surface growth coefficient
double dfe, kfe; // dimension fractale et préfacteur fractal
double xsurfgrowth, coeffB; // surface growth parameter, Bêta
int puissancep;
double lambda, Dpeqmass, rpeqmass, gamma_; // libre parcours moyen d'une sphère
double Dpm, sigmaDpm; //
double temps;
double* Vectdir; // direction aléatoire
double* TriCum;
double* Translate;
double* DistTab;
double* NombreAlea;
double* TpT;
double** PosiGravite; // Position du centre de gravité
double** Aggregate; // Tableau des aggrégats
double** tab;
int DeltaSauve;
int NSauve;
int Mode;
int NumAgg;
int compteur;
int nb_line;
int secondes;
int NAgg; // Nombre d'aggrégats (1 à l'initialisation)
int iValTab=0;
int* Monoi; //
int* MonoSel; //  Tableaux d'indices de sphères appartenant à un aggrégat
int* MonoRep; //
double** IdPossible;
int* IdDistTab;
int* IndexPourTri;
bool* Select;
char com[500];
bool with_dots;
int* TamponValeurs;
int** AggLabels;
bool root_dicho = false;
bool root_sec = true;
//double precision = 1E-1; //Précision recherchée
double precision = 1E-4; //Précision recherchée
//double precision = 1E-12; //Précision recherchée
double FactorModelBeta;
MainWindow* GUI;

Sphere* spheres;
Sphere s1,s2;

void InitRandom()
{
    time_t t;
    time(&t);
    //srand(t);
    srand(0);
}

double Random()
{
    double v = rand();
    v = v/RAND_MAX;
    return v;
}

double Maxi2D(int colonne, int nmax)
{
    //Maximum of a column in the Aggregate table
    int i;
    double m = Aggregate[1][colonne];

    for (i = 2;i <= nmax; i++)
        if (Aggregate[i][colonne] > m)    m = Aggregate[i][colonne];

    return m;
}

double MinEtIndex(double* tableau, int size, int& position)
{
    //Minimum of a table
    int i;
    double m;
    position = 1;
    m = tableau[position];

    for (i = 2;i <= size; i++)
    {
        if (tableau[i] < m)
        {
            position = i;
            m = tableau[position];
        }
    }

    return m;
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

double erfccheb(double z)
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
double erfc(double x)
{
                if (x >= 0.) return erfccheb(x);
                else return 2.0 - erfccheb(-x);
}

double inverfc(double p)
{
        double x,err,t,pp;
        if (p >= 2.0) return -100.;
        if (p <= 0.0) return 100.;
        pp = (p < 1.0)? p : 2. - p;
        t = sqrt(-2.*log(pp/2.));
        x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
        for (int j=0;j<2;j++) {
                err = erfc(x) - pp;
                x += err/(1.12837916709551257*exp(-(x*x))-x*err);
        }
        return (p < 1.0? x : -x);
}
double erf(double x) { return 1-erfc(x); }
double inverf(double p) {return inverfc(1.-p);}
//###############################################################################################################################

double Cunningham(double R) //Facteur correctif de Cunningham
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

double ModeleBeta(double rm, double np, double rg)
{
    if (rm<0.) {throw 42;}

    //return rm*Cunningham(rpeqmass) - rg*Cunningham(rm)*pow(kfe,1/dfe)*pow(np,(gamma_-1)/dfe) ;

    return Cunningham(rm) - FactorModelBeta*rm;
    //return rm - FactorModelBeta*Cunningham(rm);
}

double Dichotomie(double np, double rg,double rpmoy,double x0)
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

double brentq(double np, double rg,double rpmoy,double x0)
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
            if (2*fabs(stry) < fmin(fabs(spre), 3*fabs(sbis) - delta)) {
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

double Secante(double np,double rg,double rpmoy, double xsave)
{
    double 	x0,fx0,x1, x2, fx1;
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

double ConvertRg2Dm(double np, double rg,double rpmoy)
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
//################################################## Recherche de sphères #############################################################################

int SelectLabelEgal(int id, int* resu)
{
    int i, n;

    n = 0;
/*
    for (i = 1; i <= N; i++)
        if (spheres[i].Label == id)
        {
            n++;
            resu[n] = i;
        }
*/
    for(i=1;i<=AggLabels[id][0];i++)
    {
        resu[i]=AggLabels[id][i];
    }
    n=AggLabels[id][0];

    return n;
}

int SelectLabelSuperieur(int id, int* resu)
{
    int i,j;
    /*
    int n = 0;
    for (i = 1;i <= N; i++)
        if (spheres[i].Label > id)
        {
            n++;
            resu[n] = i;
        }

    return n
    */
    int m = 0;

    for(i=id;i<=NAgg;i++)
    {

        for(j=1;j<=AggLabels[i][0];j++)
        {

            resu[m+j]=AggLabels[i][j];
        }
        m+=AggLabels[i][0];
    }
    return m;
}

//############# Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ##############

double RayonGiration(int id, double &rmax, double &Tv, int &Nc, double &cov, double &volAgregat, double &surfAgregat, double** PosiGravite)
{
    double dist, rpmoy, dbordabord, li, r, Arg, Brg, terme;
    int i, j, k, nmonoi;
    double * tabVol;
    double* tabSurf;

    cov = 0.0;
    Nc = 0;
    volAgregat = surfAgregat = terme = 0.0;
    rmax = 0.0;
    Arg = Brg = 0.0;
    Tv = 0.0;
    //$ Identification of the spheres in Agg id
    nmonoi = SelectLabelEgal(id, Monoi);

    tabVol = new double[nmonoi+1];
    tabSurf = new double[nmonoi+1];

    for (k = 1; k <= 3; k++)    PosiGravite[id][k] = 0.0; //Initialisation

    for (i = 1; i <= nmonoi; i++)
    {
        tabVol[i] =0.;
        tabSurf[i] =0.;
    }


    for (i = 1; i <= nmonoi; i++) //Pour les i sphérules constituant l'agrégat n°id
    {
        //$ Calculation of the volume and surface of monomere i of Agg id
        tabVol[i] += spheres[Monoi[i]].Volume(); //Calculation of the volume of i
        tabSurf[i] += spheres[Monoi[i]].Surface();    //Calculation of the surface of i

        for (j = i+1; j <= nmonoi; j++) //for the j spheres composing Aggregate n°id
        {

            double voli, volj, surfi, surfj;
            voli = volj = surfi = surfj = 0.;

            //$ Calculation of the intersection between the spheres i and j
            dist = spheres[Monoi[i]].Intersection(spheres[Monoi[j]],voli,volj,surfi,surfj);

            rpmoy = (spheres[Monoi[i]].Radius()+spheres[Monoi[j]].Radius())/2.0; //Mean Radius between i and j monomeres
            dbordabord = ((dist-2.0*rpmoy)/(2.0*rpmoy))*1E6; //distance between the two particles
            //$ Check if i is covering j
            //$ [dbordabord <= 1]
            if (dbordabord <= 1.)
            {
                //$ Calculation of the Number of contacts
                cov = cov - dbordabord; //Coefficient of total Covering of Agg id
                Nc = Nc + 1; //Nombre de coordination (nombre de points de contact entre deux sphérules)
            }

            //$ The volume and surface covered by j is substracted from those of i
            tabVol[i] = tabVol[i] - voli;    //Calcul du volume de la sphérule i moins le volume de
                                             //la calotte due à la sphérule j
            tabSurf[i] = tabSurf[i] - surfi; //Calcul de la surface de la sphérule i moins la surface de
                                             //la calotte due à la sphérule j

            //$ The volume and surface covered by i is substracted from those of j
            tabVol[j] = tabVol[j] - volj;    //Calcul du volume de la sphérule j moins le volume de
                                             //la calotte due à la sphérule i
            tabSurf[j] = tabSurf[j] - surfj; //Calcul de la surface de la sphérule j moins la surface de
                                             //la calotte due à la sphérule i
        }
        //$ Calculation of the total volume and surface of the aggregate
        volAgregat = volAgregat + tabVol[i];    //Total Volume of Agg id
        surfAgregat = surfAgregat + tabSurf[i]; //Total Surface of Agg id

        terme = terme + tabVol[i]/spheres[Monoi[i]].Volume();

        //$ Calculation of the position of the center of mass
        double* pos = spheres[Monoi[i]].Position();
        for (k = 1; k <= 3; k++)
            PosiGravite[id][k] = PosiGravite[id][k] + pos[k]*tabVol[i]; //Somme des Vi*xi
    }
    Tv = 1 - terme /nmonoi;

    Nc = Nc/2;//and determine the coefficient of mean covering of Agg Id
    //$ Check if there is more than one monomere in the aggregate
    //$ [nmonoi == 1]
    if (nmonoi == 1 || Nc == 0)
    {
        //$ Cov = 0
        cov = 0;
    }
    else
    {
        //$ Determination of the coefficient of mean covering, using the one determined in the precedent loop
        cov = cov/((double)Nc)/2.0;
    }
    //$ Filling of PosiGravite

    for (k = 1; k <= 3; k++)
    {
        PosiGravite[id][k] = PosiGravite[id][k]/volAgregat;
    } //Centre of mass of Agg Id
    //$ Determination of the maximal radius of Agg Id and the Radius of gyration

    for (i = 1; i <= nmonoi; i++)
    {
        //$ Determination of the distance between each monomere and the center of mass of Agg Id
        li = spheres[Monoi[i]].Distance(PosiGravite[id]); //Distance entre les centres de masse de la sphérule i et de l'agrégat n°id

        r = li + spheres[Monoi[i]].Radius();

        //$ Calculation of rmax
        rmax=fmax(rmax,r);

        //$ Calculation of Rg
        Arg = Arg + tabVol[i]*pow(li, 2);
        Brg = Brg + tabVol[i]*pow(spheres[Monoi[i]].Radius(), 2);
    }

    delete[] tabVol;
    delete[] tabSurf;

    double rg = fabs((Arg+3.0/5.0*Brg)/volAgregat);
    volAgregat=fabs(volAgregat);

    return sqrt(rg); //Rayon de giration de l'agrégat n°id
}
//###############################################################################################################################

//######### Mise à jour des paramètres physiques d'un agrégat (rayon de giration, masse, nombre de sphérules primaires) #########
void ParametresAgg(int Agg)
{
    int i, np, Nc;
    double masse, dm, cc, diff, vit, lpm, Tv, cov, rmax, rg, rpmoy, rpmoy2, rpmoy3;
    double volAgregat, surfAgregat;

    rpmoy = rpmoy2 = rpmoy3 = 0.0;
    rg = RayonGiration(Agg, rmax, Tv, Nc, cov, volAgregat, surfAgregat, PosiGravite);

    masse = Rho*volAgregat; //Masse réelle de l'agrégat n°Agg

    np = SelectLabelEgal(Agg, MonoSel); //Liste des sphérules constituant l'agrégat n°Agg
    //npeqmass = kfe*pow(rg/rpeqmass,dfe);

    for (i = 1; i <= np; i++)
    {
        rpmoy = rpmoy + spheres[MonoSel[i]].Radius(); //Somme des rayons des sphérules de l'agrégat n°Agg
        rpmoy2 = rpmoy2 + pow(spheres[MonoSel[i]].Radius(), 2);
        rpmoy3 = rpmoy3 + pow(spheres[MonoSel[i]].Radius(), 3);
    }


    rpmoy = rpmoy/((double)np);   //Calcul du rayon moyen de l'agrégat n°Agg
    rpmoy2 = rpmoy2/((double)np); //Calcul du rayon moyen d'ordre 2 de l'agrégat n°Agg
    rpmoy3 = rpmoy3/((double)np); //Calcul du rayon moyen d'ordre 3 de l'agrégat n°Agg
    //printf("Np   %d   Rpmoy  %e   lambda  %e   Gamma  %e   dfe  %e\n",np,rpmoy,lambda,gamma_,dfe);
    dm = ConvertRg2Dm(np,rg,rpmoy);
    cc = Cunningham(dm/2);
    diff = K*T/3/PI/Mu/dm*cc;
    vit = sqrt(8*K*T/PI/masse);
    lpm = 8*diff/PI/vit;

    Aggregate[Agg][0] = rg; //Rayon de giration
    Aggregate[Agg][1] = np; //Nombre de spherules par aggregat
    Aggregate[Agg][2] = Nc; //Nombre de coordination
    Aggregate[Agg][3] = dm; //Diametre de mobilité

    if (ActiveModulephysique == 1)
    {
        Aggregate[Agg][4] = lpm;     //Libre parcours moyen
        Aggregate[Agg][5] = lpm/vit; //Durée du deplacement
    }
    else
    {
        Aggregate[Agg][4] = Dpm*1E-9;
        Aggregate[Agg][5] = 1E-6;
    }

    Aggregate[Agg][6] = rmax;          //Rayon de la sphere qui contient l'agregat
    Aggregate[Agg][7] = volAgregat;    //Estimation du volume de l'agrégat
    Aggregate[Agg][8] = surfAgregat;   //Estimation de la surface libre de l'agrégat
    Aggregate[Agg][9] = np*4.0*PI*rpmoy3/3.0; //Volume de l'agrégat sans recouvrement      (Avant c'était Tv : Taux de recouvrement volumique)
    Aggregate[Agg][10] = cov;          //Paramètre de recouvrement
    Aggregate[Agg][11] = np*4.0*PI*rpmoy2;  //Surface libre de l'agrégat sans recouvrement       (Avant c'était surfAgregat/volAgregat : Estimation du rapport surface/volume de l'agrégat)
}
//###############################################################################################################################

//############################################# Conditions aux limites périodiques ##############################################
void ReplacePosi(int id)
{
    int i,j,k,nr;

    nr = SelectLabelEgal(id,MonoRep);

    double* trans = new double[4];
    bool move=false;

    for (i = 1; i <= 3; i++)
    {
        if (PosiGravite[id][i] > L)
        {
            trans[i] = - L;
            move = true;
        }
        else if (PosiGravite[id][i] < 0)
        {
            trans[i] = L;
            move = true;
        }
        else
        {
            trans[i] = 0;
        }
    }
    if (move)
    {
        for (i = 1; i <= 3; i++)
            PosiGravite[id][i] = PosiGravite[id][i] + trans[i];

        for (j = 1; j <= nr; j++)
        {
            spheres[MonoRep[j]].Translate(trans);
        }
    }
    delete[] trans;
}
//###############################################################################################################################

void SupprimeLigne(int ligne)
{
    int i, j;

    for (i = ligne + 1; i<= N; i++)
    {
        for (j = 0; j <= 11; j++)
            Aggregate[i-1][j] = Aggregate[i][j];
        for (j = 1; j<= 3; j++)
            PosiGravite[i-1][j] = PosiGravite[i][j];

    }

    for (i=ligne+1;i<=NAgg;i++)
    {
        delete[] AggLabels[i-1];
        AggLabels[i-1]=new int[AggLabels[i][0]+1];
        for (j=0;j<=AggLabels[i][0];j++)
        {
            AggLabels[i-1][j]=AggLabels[i][j];
        }
    }



    NAgg--;
}

void InsertionSort(int n, double arr[], int index[])
{
    int i, j, id;
    double a;
    //$ TOTO

    for (i = 1; i <= n; i++)
        index[i] = i;
    for (j = 2; j <= n; j++)
    {
        a = arr[j];
        id = index[j];
        i = j-1;

        while (i > 0 && arr[i] > a)
        {
            arr[i+1] = arr[i];
            index[i+1] = index[i];
            i--;
        }

        arr[i+1] = a;
        index[i+1] = id;
    }
}



void quickSort(double arr[], int index[], int left, int right) {

      int i = left, j = right;
      int itmp;
      double pivot = arr[(left + right) / 2];
      double dtmp;

      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  dtmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = dtmp;
                  itmp = index[i];
                  index[i] = index[j];
                  index[j] = itmp;
                  i++;
                  j--;
            }
      };

      /* recursion */
      if (left < j)
            quickSort(arr, index, left, j);
      if (i < right)
            quickSort(arr, index, i, right);
}

void quickSort(int n, double arr[], int index[])
{
      quickSort(arr, index, 1, n);
}

void MonTri(int n, double arr[], int index[])
{
    int i, j, id;
    double a;

    for (i = 1; i <= n; i++)
        index[i] = i;

    //InsertionSort(n, arr, index);
    quickSort(n, arr, index);

}

int Probabilite(bool trier,double &deltatemps)
{
  /*!
   *  \image html cFProbabilitebd.png
   */


    double valAlea,max;
    int i,n,nret;

    //$ Get the maximum timestep
    max = Maxi2D(5,NAgg);

    if (trier)
    {
    //$ Sort the timesteps
        for (i=1; i <= NAgg; i++)
            TpT[i] = max/Aggregate[i][5];

        MonTri(NAgg, TpT, IndexPourTri); //$

    //$ Accumulate the timesteps
        TriCum[1] = TpT[1];
        for (i=2; i <= NAgg; i++)
            TriCum[i] = TriCum[i-1]+TpT[i];
    }
    else
    {
    //$ Keep previously sorted TriCum
    }

    //$ Pick a random sphere
    valAlea=Random()*TriCum[NAgg];
    n = 0;
    for (i=1; i<= NAgg; i++)
    {
        if (TriCum[i] < valAlea)
            n++;
    }
    n++;

    if (n > NAgg)     n = NAgg;

    nret = IndexPourTri[n];
    deltatemps = max/TriCum[NAgg];

    return nret;
}

//############################################# Calcul de la distance inter-agrégats ############################################
double Distance_Aggregate(int s, int nmonoi, double lpm)
{
    int agg, j, k, knum, np, l;
    double dist, dc, ret;
    Sphere spheredecale;

    ret = 1.0;
    agg = IdPossible[s][1];

    np = SelectLabelEgal(agg, MonoSel); //Liste des sphérules constituant l'agrégat n°Agg

    double* trans = new double[4];
    trans[1] = IdPossible[s][2]*L;
    trans[2] = IdPossible[s][3]*L;
    trans[3] = IdPossible[s][4]*L;

    //$ Loop on all the spheres of the other aggregate
    for (j = 1; j <= np; j++)
    {
        //$ spheredecale is used to replace the sphere into the corresponding box
        spheredecale.Update(spheres[MonoSel[j]]);
        spheredecale.Translate(trans);

        //$ For every sphere in the aggregate :
        for (knum = 1; knum <= nmonoi; knum++)
        {
            k=Monoi[knum];
            //$ Check if j and k are contact and if not, what is the distance between them
            dist=spheres[k].Collision(spheredecale, Vectdir, lpm, dc);
            if (dist < ret)
            {
                ret = dist;
            }

        }
    }
    delete[] trans;
    return ret;
}
//###############################################################################################################################

//########################################## Determination of the contacts between agrgates ##########################################
void CalculDistance(int id, double &distmin, int &aggcontact)
{
    double lpm,dist;
    double dc;
    int nmonoi;
    int i,s;
    int npossible;
    int dx,dy,dz;
    nmonoi =0;
    lpm = Aggregate[id][4]; // mean free path of the agregate labeled id
    npossible = 0;
    aggcontact = 0;
    distmin = 1.0;

    //+++++++++++++++++++++++++++++++++ Determination of the potential contacts (= in this part, we're considering aggreghates as spheres  with a diameter of Aggregate[6],+++++++++++++++++
    //+++++++++++++++++++++++++++++++++ the distance between the center of gravity of the aggregate and the furthest sphere in the agregate +++++++++++++++++++++++++++++++++++++++++++++++++++

    //$ Find potential contacts between agregates

    s1.Update(PosiGravite[id], Aggregate[id][6]); // Represents the sphere containing the agregate we're testing

    //$ [For all other agregate]
    for (i = 1;i <= NAgg; i++)
    {
        if (i != id)
        {
            double inix,iniy,iniz,inir;
            inix = PosiGravite[i][1];
            iniy = PosiGravite[i][2];
            iniz = PosiGravite[i][3];
            inir = Aggregate[i][6]; //represents the different agregates

            //$ [3 imbricated loops on dx,dy,dz to look into the 27 boxes]
            for (dx = -1;dx <= 1; dx++)
            {
                for (dy = -1; dy <= 1; dy++)
                {
                    for (dz = -1; dz <= 1; dz++)
                    {
                        s2.Update(inix+L*dx,iniy+L*dy,iniz+L*dz,inir);

                        // checks if the two spheres will be in contact while
                         //... the first one is moving
                        //$ Intersection check between agregates
                        dist = s1.Collision(s2, Vectdir, lpm, dc);

                        //$ [Potential Collision]
                        if (dist <= lpm)
                        {
                            //$ Aggregate is stocked into IdPossible
                            npossible++; // Number of aggregates that could be hit
                            IdPossible[npossible][1] = i;  //Label of an aggregate that could be in contact with the one moving
                            IdPossible[npossible][2] = dx; //X coordinate of the "box" where this agregate was
                            IdPossible[npossible][3] = dy; //Y coordinate of the "box" where this agregate was
                            IdPossible[npossible][4] = dz; //Z coordinate of the "box" where this agregate was
                        }
                    }
                }
            }
        }
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //$ Number of agregates possibly in contact
    if (npossible > 0)
    {
        nmonoi = SelectLabelEgal(id, Monoi); //Liste of the monomeres in the aggregate id
        //$ loop on the agregates potentially in contact
        for (s = 1; s <= npossible; s++) //For every aggregate that could be in contact
        {
            dist = Distance_Aggregate(s, nmonoi, lpm);
            if (dist != -1 && dist < distmin)
            {
                distmin = dist;
                aggcontact = s; //Prise en compte de l'image par translation d'un agrégat cible
            }
        }
    }
}
//###############################################################################################################################

//################################################# Réunion de deux agrégats ####################################################
int Reunit(int AggI, int AggJ, int &err)
{
    int i, nselect, numreject, numstudy,nmonoi;

    err = 0;

    nmonoi = SelectLabelEgal(AggI, Monoi); //Liste of the monomeres in the aggregate id

    for (i = 1; i <= nmonoi; i++)
            spheres[Monoi[i]].Translate(Translate);

    if (AggI < AggJ)
    {
        numstudy = AggI;
        numreject = AggJ;
    }
    else
    {
        numstudy = AggJ;
        numreject = AggI;
    }

    TamponValeurs= new int[AggLabels[AggI][0]+AggLabels[AggJ][0]+1];
    TamponValeurs[0]=AggLabels[AggI][0]+AggLabels[AggJ][0];
    for(i=1;i<=AggLabels[AggI][0];i++)
    {
        TamponValeurs[i]= AggLabels[AggI][i];
    }
    for(i=AggLabels[AggI][0]+1;i<=AggLabels[AggJ][0]+AggLabels[AggI][0];i++)
    {
        TamponValeurs[i]=AggLabels[AggJ][i-AggLabels[AggI][0]];
    }

    delete[] AggLabels[numstudy];
    AggLabels[numstudy] = new int [TamponValeurs[0]+1];
    for (i=0;i<=TamponValeurs[0];i++)
    {
        AggLabels[numstudy][i]=TamponValeurs[i];
    }



    nselect = SelectLabelEgal(numreject, MonoSel);

    for (i = 1; i <= nselect; i++)
    {
        spheres[MonoSel[i]].SetLabel(numstudy);
    }

    SupprimeLigne(numreject);
    delete[] TamponValeurs;

    nselect = SelectLabelSuperieur(numreject, MonoSel);

    for (i = 1; i <= nselect; i++)
        spheres[MonoSel[i]].DecreaseLabel();

    return numstudy;
}
//###############################################################################################################################

//################################### Test de superposition des sphérules dans un agrégat #######################################
int CalculSuperposition(int id)
{
    double d, dsuperpos, mem;
    int i, j, nmonoi;

    nmonoi = SelectLabelEgal(id, Monoi);
    mem = 0;

    for (i = 1;i <= nmonoi; i++)
    {
        for (j = 1;j <= i-1; j++)
        {
            d = spheres[Monoi[j]].Distance(spheres[Monoi[i]]);
            double rj = spheres[Monoi[j]].Radius();
            double ri = spheres[Monoi[i]].Radius();
            dsuperpos = (d-(rj+ri))/(rj+ri);
            if (dsuperpos < mem)
            {
                mem = dsuperpos;
            }
        }
    }
    mem = mem*1E6;

    if ((int)mem < 0)
    {
        printf("Superposition pour Numfichier %d,   agrégat n°%d    mem=%d\n", NSauve, id, (int)mem);
        return 1;
    }
    else return 0;
}
//###############################################################################################################################

//####################################### Croissance de surface des particules primaires ########################################
void CroissanceSurface(double dt)
{
    int i;

    for (i = 1; i <= N; i++)
    {
        double newradius = spheres[i].Radius() + Asurfgrowth*pow(spheres[i].Radius(), xsurfgrowth-2)*dt;
        spheres[i].Update(spheres[i].Position(),newradius );
    }
}
//###############################################################################################################################

//################################ Lecture des données physiques d'entrée variables dans le temps ###############################
int LectureSuiviTempo()
{
    FILE* f = NULL;
    f = fopen(qPrintable(FichierSuiviTempo), "rt");
    int i, j, k;
    char skip_line[500], data[500];

    nb_line=0;

    if (f != NULL)
    {
        //On lit la première ligne pour la sauter
        fgets(skip_line, 500, f);

        while(fgets(data, 500, f) != NULL)
        {
            nb_line++; //On compte le nombre de lignes ...
        }

        printf("\nLe fichier contient %d lignes.\n", nb_line); //... et on l'affiche à l'écran
        printf("\n");

        //On se replace au début du fichier ...
        fseek(f, 0, SEEK_SET);

        //... puis on relit la première ligne pour la sauter
        fgets(skip_line, 500, f);

        //On déclare un tableau de type double de deux colonnes qui contiendra les valeurs lues
        tab = new double* [nb_line];
        for(k = 0; k < nb_line; k++)
        {
            tab[k] = new double[2];
        }

        //On lit chaque ligne suivante du fichier
        for (i = 0; i < nb_line; i++)
        {
            //On lit la ligne
            if (fgets(data, 500, f) != NULL)
            {
                //On crée une chaine qui contiendra le token
                char *token;

                //On découpe selon les tabulations
                token = strtok(data, "\t");

                while(token != NULL)
                {
                    for (j = 0; j < 2; j++)
                    {
                        //On enregistre la valeur convertie en double dans le tableau
                        tab[i][j] = atof(token);

                        //On affiche le tableau
                        //printf("ligne %d : tableau[%d][%d] = %.6f\n",i, i, j, tab[i][j]);

                        token = strtok(NULL, "\t");
                    }
                }
            }
        }
        /*
        for (i = 0; i < nb_line; i++)
        {
            for (j = 0; j < 2; j++)
            {
                cout << tab[i][j] << "  ";
            }
            cout << endl;
        }
        */
        fclose(f);
        printf("\nLecture terminee\n");
    }
    else
    {
        cout << "Impossible d'ouvrir le fichier de donnees de suivi temporel" << endl;
        exit(1);
    }

    return nb_line;
}
//###############################################################################################################################

//########################################### Recherche d'une valeur dans un tableau ############################################
int rechercheValTab() //Programme qui cherche les données physiques dans le tableau de suivi temporel et qui retourne 1 si on atteint la limite de ce fichier
{
    int test = 0;

    while (tab[iValTab][0] < temps && iValTab < (nb_line-1))
    {
        iValTab++;
//        cout << "temps residence = " << temps << "s"
//             << "       temps fichier = " << tab[iValTab][0] << "s"
//             << "    T = " << tab[iValTab][1] << "K" << endl;
    }

    T = tab[iValTab-1][1];

    if (iValTab == nb_line-1)     test = 1;

    return test;
}
//###############################################################################################################################

void Init()
{
    int i, j, k, test, testmem;
    double x, masse, surface, Dp, Cc, Diff, Vit, lpm, rg, dist;

    compteur = 0;
    testmem = 0;
    L = X*Dpm*1E-9;
    Mu = 18.203E-6*(293.15+110)/(T+110)*pow(T/293.15,1.5);
    K = 1.38066E-23;
    lambda = 66.5E-9*(101300/P)*(T/293.15)*(1+110/293.15)/(1+110/T);
    Dpeqmass = Dpm*exp(1.5*log(sigmaDpm)*log(sigmaDpm)); //Diamètre équivalent massique moyen des monomères
                                                           //donné par l'équation de Hatch-Choate
    rpeqmass = (Dpeqmass*1E-9)/2.0; //Rayon équivalent massique moyen des monomères
    gamma_ = 1.378*(0.5+0.5*erf(((lambda/rpeqmass)+4.454)/10.628));

    spheres = new Sphere[N+1];
    Translate = new double[4];
    Vectdir = new double[4];
    TriCum = new double[N+1];
    Select = new bool[N+1];
    Monoi = new int[N+1];
    MonoSel = new int[N+1];
    MonoRep = new int[N+1];
    IdPossible = new double* [N+1];
    DistTab = new double[N+1];
    IdDistTab = new int[N+1];
    PosiGravite = new double* [N+1];
    Aggregate = new double* [N+1];
    TpT = new double[N+1];
    IndexPourTri = new int[N+1];

    AggLabels= new int*[N+1];

    for (i=1;i<=N;i++)
    {
        AggLabels[i]= new int[2];
        AggLabels[i][0]=1;
        AggLabels[i][1]=i;

    }
    AggLabels[0]=new int[2];
    AggLabels[0][0]=0;
    AggLabels[0][1]=0;

    for (i = 1; i <= N; i++)
    {
        PosiGravite[i] = new double[4];
        Aggregate[i] = new double[12];
        IdPossible[i] = new double[5];
    }

    for (i = 1; i <= N; i++)
    {

        for (j = 1; j<= 3 ; j++)
            PosiGravite[i][j] = Random()*L;

        x = Random(); //Tirage aléatoire

        if (Mode == 1)
            Dp = Dpm+sqrt(2.0)*sigmaDpm*inverf(2*x-1); //Loi normale
        else
            Dp = exp(log(Dpm)+sqrt(2.0)*log(sigmaDpm)*inverf(2*x-1)); //Loi log-normale

        Dp = Dp/1E9;

        if (Dp <= 0)  Dp = Dpm*1E-9;

        spheres[i].Update(PosiGravite[i], Dp/2);
        spheres[i].SetLabel(i);

        //++++++++++++ Test de superposition des sphérules lors de leur génération aléatoire ++++++++++++
        test=0;
        for (k = 1; k <= i-1; k++)
        {
            dist = spheres[k].Distance(spheres[i]); // Calcule la distance centre à centre entre le monomère k et tous les autres
            if (dist <= spheres[k].Radius()+spheres[i].Radius())
                test++;
        }

        if (test > 0)
            i--;

        testmem = testmem + test; //Comptabilise le nombre d'échecs à positionner une sphère sans superposition

        if (testmem > N)
        {
            printf("Impossible de générer tous les monomères sans superposition.\n");
            printf("La fraction volumique doit être diminuée.\n");
            exit(0);
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        masse = Rho*PI*pow(Dp,3)/6;
        surface = PI*pow(Dp,2);
        Cc = Cunningham(Dp/2);
        Diff = K*T/3/PI/Mu/Dp*Cc;
        Vit = sqrt(8*K*T/PI/masse);
        lpm = 8*Diff/PI/Vit;
        rg = sqrt(3.0/5.0)*Dp/2;

        if (ActiveModulephysique==1)
        {
            Aggregate[i][0] = rg;                   //Rayon de gyration
            Aggregate[i][1] = 1;                    //Nommbre de sphérules par agrégat Np
            Aggregate[i][2] = 0;                    //Nombre de coordination Nc
            Aggregate[i][3] = Dp;                   //Diamètre de mobilité
            Aggregate[i][4] = lpm;                  //Libre parcours moyen
            Aggregate[i][5] = lpm/Vit;              //Durée du déplacement
            Aggregate[i][6] = Dp/2;                 //Rayon de la sphère d'enveloppe de l'agrégat réunifié
            Aggregate[i][7] = masse/Rho;            //Volume estimé de l'agrégat réunifié
            Aggregate[i][8] = surface;              //Surface estimée de l'agrégat réunifié
            Aggregate[i][9] = PI*pow(Dp, 3)/6;      //Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
            Aggregate[i][10] = 0;                   //Coefficient de pénétration (paramètre de recouvrement Cov)
            Aggregate[i][11] = PI*pow(Dp, 2);       //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)
        }
        else
        {
            Aggregate[i][0] = rg;                   //Rayon de gyration
            Aggregate[i][1] = 1;                    //Nommbre de sphérules par agrégat Np
            Aggregate[i][2] = 0;                    //Nombre de coordination Nc
            Aggregate[i][3] = Dp;                   //Diamètre de mobilité
            Aggregate[i][4] = Dpm*1E-9;             //Libre parcours moyen
            Aggregate[i][5] = 1E-6;                 //Durée du déplacement
            Aggregate[i][6] = Dp/2;                 //Rayon de la sphère d'enveloppe de l'agrégat réunifié
            Aggregate[i][7] = masse/Rho;            //Volume de l'agrégat réunifié
            Aggregate[i][8] = surface;              //Surface de l'agrégat réunifié
            Aggregate[i][9] = PI*pow(Dp, 3);        //Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
            Aggregate[i][10] = 0;                   //Coefficient de pénétration (paramètre de recouvrement Cov)
            Aggregate[i][11] = PI*pow(Dp, 2);       //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)
        }
    }

    NAgg = N;
}

void Fermeture()
{

    for (int i = 0; i <= N; i++)
    {
        delete[] PosiGravite[i];
        delete[] Aggregate[i];
        delete[] IdPossible[i];

    }
    delete[] spheres;
    delete[] Translate;
    delete[] Vectdir;
    delete[] TriCum;
    delete[] Select;
    delete[] Monoi;
    delete[] MonoSel;
    delete[] MonoRep;
    delete[] IdPossible;
    delete[] DistTab;
    delete[] IdDistTab;
    delete[] PosiGravite;
    delete[] Aggregate;
    delete[] TpT;
    delete[] IndexPourTri;
    delete[] tab;
}

void SauveASCII(int value, int id)
{
    int i,j;
    char NomComplet[500];
    FILE *f;

    locale::global(locale("C"));

    return;

    sprintf(NomComplet, "%s/Sphere%05d.txt", CheminSauve, value);
    f = fopen(NomComplet, "w");
    fprintf(f, "%d  N_[]\n", N);
    fprintf(f, "%10.3f  FV_[ppm]\n", FV*1E6);
    fprintf(f, "%10.3f  X_[]\n", X);
    fprintf(f, "%10.3f  Dpm_[nm]\n", Dpm);
    fprintf(f, "%10.3f  sigmaDpm_[nm]\n", sigmaDpm);
    fprintf(f, "%d  NAgg_[]\n", NAgg);
    fprintf(f, "%10.6f  Temps_[µs]\n", temps*1E6);
    fprintf(f, "Label\t      Rp(nm)\t      X(nm)\t     Y(nm)\t     Z(nm)\n");

    for (i=1;i<=N;i++)
        fprintf(f,"%s\n", spheres[i].str(1e9).c_str());

    fclose(f);

    sprintf(NomComplet,"%s/Agg%05d.txt", CheminSauve, value);
    f = fopen(NomComplet, "w");
    fprintf(f,"%d  N_[]\n",N);
    fprintf(f,"%1.3f  FV_[ppm]\n", FV*1E6);
    fprintf(f,"%1.3f  X_[]\n", X);
    fprintf(f,"%1.3f  Dpm_[nm]\n", Dpm);
    fprintf(f,"%1.3f  sigmaDpm_[nm]\n", sigmaDpm);
    fprintf(f,"%d  NAgg_[]\n", NAgg);
    fprintf(f,"%1.6f  Temps_[µs]\n", temps*1E6);
    fprintf(f,"Rg_[nm]\tNp_[]\tNc_[]\tDm_[nm]\tlpm_[nm]\tdeltat_[µs]\tRgeo_[nm]\tXG(nm)\tYG(nm)\tZG(nm)\tV_[1E-25m3]\tS_[1E-16m2]\tVOlWO_[1E-25m3]\tcov[]\tSurfWO_[1E-16m2]\n");

    for (i = 1; i <= NAgg; i++)
    {
        fprintf(f, "%10.3f\t", Aggregate[i][0]*1E9);
        fprintf(f, "%d\t", (int)Aggregate[i][1]);
        fprintf(f, "%d\t", (int)Aggregate[i][2]); //Nombre de coordination
        fprintf(f, "%10.3f\t", Aggregate[i][3]*1E9);
        fprintf(f, "%10.3f\t", Aggregate[i][4]*1E9);
        fprintf(f, "%10.3f\t", Aggregate[i][5]*1E6);
        fprintf(f, "%10.3f\t", Aggregate[i][6]*1E9);

        for (j = 1; j <= 3; j++)
            fprintf(f, "%10.3f\t", PosiGravite[i][j]*1E9);

        fprintf(f, "%10.3f\t", Aggregate[i][7]*1E25);  //Estimation du volume d'un agrégat
        fprintf(f, "%10.3f\t", Aggregate[i][8]*1E16);  //Estimation de la surface d'un agrégat
        fprintf(f, "%10.3f\t", Aggregate[i][9]*1E25);       //Volume d'un agrégat sans recouvrement (Avant c'était Taux de recouvrement volumique)
        fprintf(f, "%10.3f\t", Aggregate[i][10]);      //Coefficient de pénétration (paramètre de recouvrement Cov)
        fprintf(f, "%10.3f\n", Aggregate[i][11]*1E16); //Surface d'un agrégat sans recouvrement (Avant c'était Rapport surface/volume)
    }
    fclose(f);

    sprintf(NomComplet, "%s/AFractalPlot.txt", CheminSauve);
    if (Aggregate[id][1] > 10)
    {
        f=fopen(NomComplet, "a");
        fprintf(f, "%e\t%e\n", Aggregate[id][0]/(Dpm*1E-9), Aggregate[id][1]);
        fclose(f);
    }
}


void test_locale()
{
    double testfloat = 1.5;
    string teststr1 = "1.5";
    string teststr2 = "1,5";
    double test1=atof(teststr1.c_str());
    double test2=atof(teststr2.c_str());

    if (fabs(test1-testfloat)<1e-3)
        with_dots = true;
    else if (fabs(test2-testfloat)<1e-3)
            with_dots = false;
    else
    {
        printf("What locale are you using ?\n");
        exit(1);
    }
}

double latof(const char* _char)
{
    string mystring = _char;
    if (!with_dots)
    {
        int f = mystring.find(".");
        if (f>0)
            mystring.replace(f, 1, ",");
    }
    return atof(mystring.c_str());
}


void LectureParams()
{
    FILE* f;
    char sauve[500];

    //char t[500] ;
    f = fopen(qPrintable(FichierParam), "rt");
    QDir lDir;
    if (f == NULL)
    {
        N = 2500;
        T = 1500;
        Dpm = 30;
        sigmaDpm = 0.0;
        FV = 3E-3;
        Rho = 1.8E3;
        P = 101300;
        Mode = 1;
        DeltaSauve = 0;
        sprintf(sauve, "%s", "DATA");
    }
    else
    {

        test_locale();

        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        N=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        T=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        Dpm=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        sigmaDpm=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        FV=latof(commentaires)*1E-6;
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        P=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        Rho=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        Mode=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        DeltaSauve=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",sauve,com);
        fclose(f);
    }
    if (Mode == 1)
        X = pow(N*PI/6.0/FV*(1.0+3.0*sigmaDpm*sigmaDpm/Dpm/Dpm),1.0/3.0); //Loi normale
    else
        X = pow(N*PI/6.0/FV*exp(9.0/2.0*log(sigmaDpm)*log(sigmaDpm)),1.0/3.0); //Loi log-normale

    strcat(CheminSauve,qPrintable(pathParam));
    strcat(CheminSauve,"//");
    strcat(CheminSauve,sauve);
    lDir.mkdir(CheminSauve);

    sprintf(commentaires, "N=%d \nT=%1.3f \nDpm=%1.3f \nsigmaDpm=%1.3f \nFV=%1.3e\nX=%1.3f \nP=%1.3f\nMode=%d\nRho=%1.3f \nDeltaSauve=%d\nCheminSauve=%s\n", N, T, Dpm, sigmaDpm, FV, X, P, Mode, Rho, DeltaSauve, CheminSauve);
}

void Calcul() //Coeur du programme
{
    double deltatemps, distmin, lpm;
    double thetarandom, phirandom;
    int aggcontact, newnumagg, finfichiersuivitempo, finmem = 0;
    //int tmp,superpo;
    int i, j, co, err;
    bool contact;
    time_t t, t0;

    if (ActiveModulephysique)
    {
        sprintf(commentaires, "Le module physique est activé.\n");
        if (GUI == NULL)
            printf("%s",commentaires);
        else
            GUI->print(commentaires);
    }
    else
    {
        sprintf(commentaires, "Le module physique n'est pas activé.\n");
        if (GUI == NULL)
            printf("%s",commentaires);
        else
            GUI->print(commentaires);
    }

    if (ActiveVariationTempo)
    {
        LectureSuiviTempo();
        sprintf(commentaires, "Le fichier de données de suivi temporel est lu.\n");
        if (GUI == NULL)
            printf("%s",commentaires);
        else
            GUI->print(commentaires);
    }
    else
    {
        sprintf(commentaires, "Le fichier de données sélectionné est le fichier 'params.txt'.\n");
        if (GUI == NULL)
            printf("%s",commentaires);
        else
            GUI->print(commentaires);
    }

    sprintf(commentaires, "\nDimension fractale : %1.2f\nPréfacteur fractal : %1.2f\nB = %1.2f\nx = %1.2f\n", dfe, kfe, coeffB, xsurfgrowth);
    if (GUI == NULL)
        printf("%s",commentaires);
    else
        GUI->print(commentaires);

    Asurfgrowth = coeffB*1E-3;
    sprintf(commentaires, "Coefficient de croissance de surface : %e\n", Asurfgrowth);
    if (GUI == NULL)
        printf("%s",commentaires);
    else
        GUI->print(commentaires);

    InitRandom();


    Init();

    co=0;
    //superpo=0;
    NSauve=0;
    temps=0;
    deltatemps=0;
    time(&t0);
    contact=true;
    int end = fmax(5,N/200);
    int it_without_contact=0;
    int lim_it_without_contact = 100;

    end = 20;

    printf("\n");
    printf("Ending calcul when there is less than %d aggregats or %d iterations without contact\n", end,lim_it_without_contact);

    //$ Loop on the N monomeres
    while (NAgg > end && it_without_contact<100) //Pour N=1000 le calcul s'arrête quand il reste 5 agrégats
    {                       //Pour N=2000 le calcul s'arrête quand il reste 10 agrégats
        qApp->processEvents(); //Permet de rafraichir la fenêtre Qt
        time(&t);
        secondes = t-t0;

        // -- Generating a random direction --

        thetarandom = Random()*PI*2;
        phirandom = acos(1-2*Random());
        Vectdir[1] = sin(phirandom)*cos(thetarandom);
        Vectdir[2] = sin(phirandom)*sin(thetarandom);
        Vectdir[3] = cos(phirandom);

        // -- --



        if (ActiveModulephysique)
        {
            //$ Choice of an aggregate according to his MFP
            NumAgg = Probabilite(contact, deltatemps);// Choice of an agrgegate, the probability of said agrgegate to be chosen proportionnal to his lpm

            temps = temps + deltatemps; // Time incrementation with a value given by Probabilite

            if (ActiveVariationTempo)
            {
                //$ ?
                finfichiersuivitempo = rechercheValTab();
                if (finmem != finfichiersuivitempo)
                {
                    sprintf(commentaires, "Attention le suivi temporel est plus long que celui du fichier lu.\n");
                    if (GUI == NULL)
                        printf("%s",commentaires);
                    else
                        GUI->print(commentaires);
                }
                finmem = finfichiersuivitempo;
            }

            //$ Surface Growth
            CroissanceSurface(deltatemps);

            //$ Aggregates parameter update

            for (i = 1; i<= NAgg; i++)
            {
                ParametresAgg(i);
            }
            //$ looking for potential contacts
            CalculDistance(NumAgg, distmin, aggcontact);
            lpm = Aggregate[NumAgg][4];
            contact = (aggcontact != 0);
        }
        else
        {
            //$ Random Choice of an aggregate
            NumAgg = int(Random()*(double)NAgg)+1;
            deltatemps = 0.0;
            temps = temps + 1E-9;

            //$ looking for potential contacts

            CalculDistance(NumAgg, distmin, aggcontact);

            lpm = Dpm*1E-9; //On fixe le lpm au Dpm
            contact = (aggcontact != 0);
        }

        if (contact)
        {
            //$ Translation of the aggregate
            for (i = 1; i <= 3; i++)
            {
                Translate[i] = Vectdir[i]*distmin;
            }
            //$ The aggregate in contact is moved to the right box

            int nmonoi = SelectLabelEgal(IdPossible[aggcontact][1], Monoi);

            double* trans = new double[4];
            trans[1] = IdPossible[aggcontact][2]*L;
            trans[2] = IdPossible[aggcontact][3]*L;
            trans[3] = IdPossible[aggcontact][4]*L;

            //the aggregate that's been tested in contact with the one moving is replaced in the "box" of the space where the contact happened
            //It gets in box on one side, then has to go out on the other one.
            for (j = 1; j <= nmonoi; j++)
                    spheres[Monoi[j]].Translate(trans);
            delete[] trans;

            //$ Aggregates in contact are reunited
            newnumagg = Reunit(NumAgg, IdPossible[aggcontact][1], err);

            ParametresAgg(newnumagg);

            temps = temps-deltatemps*(1-distmin/lpm);
            /*
            //+++++ Test de superposition des monomères dans un agrégat lors d'un contact +++++
            tmp = CalculSuperposition(newnumagg);
            if (tmp == 1)     {superpo++;}
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            */
             printf("NAgg=%d  temps=%5.1f E-6 s     CPU=%d sec    after %d it\n", NAgg, temps*1E6, secondes,it_without_contact);
            it_without_contact = 0;
             if (!(GUI == NULL)            )
                GUI->progress(N-NAgg+1);

        }
        else
        {
            //$ Translation of the aggregate on his full lpm distance
            for (i = 1; i <= 3; i++)
                Translate[i] = Vectdir[i]*lpm;

            int nmonoi = SelectLabelEgal(NumAgg, Monoi);

            for (i = 1; i <= nmonoi; i++)
            {
                spheres[Monoi[i]].Translate(Translate);

            }
            for (j = 1; j <= 3; j++)
            {
                PosiGravite[NumAgg][j] = PosiGravite[NumAgg][j] + Translate[j];
            }

            newnumagg = NumAgg;
            it_without_contact++;
        }
        ReplacePosi(newnumagg);
        if (DeltaSauve>0)
        {
            co++;
            if (co > DeltaSauve)
            {
                co = 0;
                if (contact == false)
                    SauveASCII(NSauve++, newnumagg);
            }
        }

        if (contact && DeltaSauve >= 0)
            SauveASCII(NSauve++, newnumagg);
    }

    printf("Nombre total d'aggregats : %d\n",NAgg);
/*
    cout << "L=" << L*1E9
         <<"     lambda=" << lambda*1E9
         <<"     Dpeqmass=" << Dpeqmass
         <<"     rpeqmass=" << rpeqmass*1E9
         <<"     gamma=" << gamma <<endl;

*/
    for (i = 1; i <= NAgg; i++)
    {
        printf("%d\t", i);
        for (j = 1; j <= 3; j++)
            printf("%e\t", PosiGravite[i][j]*1E9);
        printf("\t%e\t%e\n",Aggregate[i][7]*1E25,Aggregate[i][9]*1E25);
    }

    Fermeture();
    if (GUI == NULL)
        printf("%s",CheminSauve);
    else
        GUI->print(CheminSauve);

    sprintf(commentaires,"\nFin du calcul  ...\n");
    if (GUI == NULL)
        printf("%s",commentaires);
    else
        GUI->print(commentaires);
}







MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButton, SIGNAL(clicked()), this, SLOT(BoutonQuitter()));
    connect(ui->BoutonRep, SIGNAL(clicked()), this, SLOT(BoutonRechercheParam()));
    connect(ui->BoutonExecDLCA, SIGNAL(clicked()), this, SLOT(ExecuterDLCA()));
    connect(ui->ActModulePhysique, SIGNAL(clicked()), this, SLOT(ModulePhysique()));
    connect(ui->SuiviTempo, SIGNAL(clicked()), this, SLOT(BoutonRechercheSuiviTempo()));
    GUI = this;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::BoutonQuitter()
{
    this->close();
    exit(0);
}

void MainWindow::BoutonRechercheParam()
{
    ui->progressBar->setVisible(false);
    FichierParam = QFileDialog::getOpenFileName(this,"Sélectionner le fichier de données DLCA",
                                                        "C:/Users/dlca/Desktop/DLCA_sous_Qt/",
                                                        "Fichier de paramètres (*.txt)");
    QFileInfo tmp2 = FichierParam;
    pathParam = tmp2.absolutePath(); //Cette variable ne retient que le chemin du fichier param

    //Affiche le chemin dans la zone de texte
    ui->AfficheurRep->append(FichierParam);
    LectureParams();
    ui->AfficheurRep->append(commentaires);
    ui->ActModulePhysique->setEnabled(true);
    ui->BoutonExecDLCA->setEnabled(true);
}


void MainWindow::BoutonRechercheSuiviTempo()
{
    FichierSuiviTempo = QFileDialog::getOpenFileName(this,"Sélectionner le fichier de suivi temporel",
                                                        "C:/Users/dlca/Desktop/DLCA_sous_Qt/",
                                                        "Fichier de paramètres (*.txt)");
    QFileInfo tmp3 = FichierSuiviTempo;
    pathSuiviTempo = tmp3.absolutePath(); //Cette variable ne retient que le chemin du fichier Suivi Tempo

    //Affiche le chemin dans la zone de texte
    ui->AfficheurRep->append(FichierSuiviTempo);
}

void MainWindow::ModulePhysique()
{
    ui->SuiviTempo->setEnabled(true);

    ui->DoubleB->setEnabled(true);
    ui->Doublex->setEnabled(true);
    ui->DimensionFractale->setEnabled(true);
    ui->PrefacteurFractal->setEnabled(true);
    ui->doubleSpinBoxCoeffCroissance->setEnabled(true);
    ui->doubleSpinBoxPuissancex->setEnabled(true);
    ui->doubleSpinBoxDimFractale->setEnabled(true);
    ui->doubleSpinBoxPrefFractal->setEnabled(true);
}

void MainWindow::ExecuterDLCA()
{

    ui->progressBar->setMaximum(N);
    ui->progressBar->setVisible(true);
    ui->progressBar->setValue(0);

    coeffB = ui->doubleSpinBoxCoeffCroissance->value();
    xsurfgrowth = ui->doubleSpinBoxPuissancex->value();
    dfe = ui->doubleSpinBoxDimFractale->value();
    kfe = ui->doubleSpinBoxPrefFractal->value();

    if (ui->ActModulePhysique->isChecked())     ActiveModulephysique = 1;
    else    ActiveModulephysique = 0;

    if (ui->SuiviTempo->isChecked())    ActiveVariationTempo = 1;
    else    ActiveVariationTempo = 0;

    Calcul();

    ui->progressBar->setVisible(false);
}

void MainWindow::print(char* str){
    ui->AfficheurRep->append(str);
}
void MainWindow::progress(int value){
    ui->progressBar->setValue(value);
}





int No_GUI(int argc, char *argv[]){
    FichierParam = argv[1];
    QFileInfo tmp2 = FichierParam;
    pathParam = tmp2.absolutePath(); //Cette variable ne retient que le chemin du fichier param
    LectureParams();

    coeffB = 1.0;
    xsurfgrowth = 2.0;
    dfe = 1.80;
    kfe = 1.40;
    ActiveModulephysique = 1;
    ActiveVariationTempo = 0;

    Calcul();

    return 0;
}

