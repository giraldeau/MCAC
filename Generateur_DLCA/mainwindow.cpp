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
double Asurfgrowth; // ?
double dfe, kfe; // dimension fractale et préfacteur fractal
double xsurfgrowth, coeffB; // N, Bêta
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


Sphere* spheres;
Sphere s1,s2;

void InitRandom()
{
    time_t t;
    time(&t);
    srand(t);
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

double Distance(double* V1, double* V2) //Calcule la distance centre à centre entre deux sphères
{
    double dx=V2[1]-V1[1];
    double dy=V2[2]-V1[2];
    double dz=V2[3]-V1[3];

    return sqrt(dx*dx+dy*dy+dz*dz);
}

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
    //double df = 1.72;
    //double kf = 1.47; //pow(5./3.,df/2.0);

    return rg/rm*Cunningham(rm) - pow(kfe,-1/dfe)*pow(np,(1-gamma_)/dfe)*Cunningham(rpeqmass);
}

double Dichotomie (double np, double rg)
{
    double 	rmin, rmax, rmed, frmed, frmin, frmax, precision;

    rmin = 0.0;   //pow(np/1.5,1/1.8)*rp/40; //borne inférieure de rm
    rmax = 5E-6; //pow(np/1.5,1/1.8)*rp*40; //bornes de recherche de rm
    precision = 0.01E-9; //Précision recherchée

    frmin = ModeleBeta(rmin, np, rg);
    frmax = ModeleBeta(rmax, np, rg);

    if (frmin*frmax>=0) {printf("Intervalle incorrect.\n"); return -1;} //Intervalle incorrect
    while (rmax-rmin>precision)
    {
        rmed = (rmin+rmax)/2 ;
        frmed = ModeleBeta(rmed, np, rg);
        if (frmed==0)   break; //zero trouvé
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
    }

    return rmed;
}

double ConvertRg2Dm(double np, double rg)
{
    return  Dichotomie (np, rg)*2; //Retourne le diamètre de mobilité
}
//################################################## Recherche de sphères #############################################################################

int SelectLabelEgal(int id, int* resu)
{
    int i, n;
    n = 0;

    for (i = 1; i <= N; i++)
        if (spheres[i].Label == id)
        {
            n++;
            resu[n] = i;
        }

    return n;
}

int SelectLabelSuperieur(int id, int* resu)
{
    int i, n;
    n = 0;

    for (i = 1;i <= N; i++)
        if (spheres[i].Label > id)
        {
            n++;
            resu[n] = i;
        }

    return n;
}

//############# Calcul du volume, de la surface, du centre de masse et du rayon de giration réévalués d'un agrégat ##############

double RayonGiration(int id, double &rmax, double &Tv, int &Nc, double &cov, double &volAgregat, double &surfAgregat, double** PosiGravite)
{
    double dist, rpmoy, dbordabord, li, r, Arg, Brg, terme;
    int i, j, k, nmonoi;
    double* tabVol;
    double* tabSurf;

    cov = 0.0;
    Nc = 0;
    volAgregat = surfAgregat = terme = 0.0;
    rmax = 0.0;
    Arg = Brg = 0.0;
    Tv = 0.0;
    nmonoi = SelectLabelEgal(id, Monoi);

    tabVol = new double [nmonoi+1];
    tabSurf = new double [nmonoi+1];

    for (k = 1; k <= 3; k++)    PosiGravite[id][k] = 0.0; //Initialisation

    for (i = 1; i <= nmonoi; i++) //Pour les i sphérules constituant l'agrégat n°id
    {
        tabVol[i] = 4.0*PI*pow(spheres[Monoi[i]].r, 3.0)/3.0; //Calcul du volume de chaque sphérule i
        tabSurf[i] = 4.0*PI*pow(spheres[Monoi[i]].r, 2.0);    //Calcul de la surface de chaque sphérule i

        for (j = 1; j <= nmonoi; j++) //Pour les j sphérules constituant l'agrégat n°id
        {
            if (i != j) //On évite l'interpénétration d'une sphérule avec elle-même
            {
                dist = spheres[Monoi[i]].Distance(spheres[Monoi[j]]);  //Calcul de la distance centre à centre entre deux monomères
                rpmoy = (spheres[Monoi[i]].r+spheres[Monoi[j]].r)/2.0; //Calcul du rayon moyen entre deux monomères
                dbordabord = ((dist-(spheres[Monoi[i]].r+spheres[Monoi[j]].r))/(spheres[Monoi[i]].r+spheres[Monoi[j]].r))*1E6; //Calcul de la distance inter-particulaire

                if ((int)dbordabord <= 1)
                {
                    cov = cov + ((2.0*rpmoy-dist)/(2.0*rpmoy)); //Coefficient de recouvrement total de l'agrégat n°id
                    Nc = Nc + 1; //Nombre de coordination (nombre de points de contact entre deux sphérules)
                }

                tabVol[i] = tabVol[i] - spheres[Monoi[i]].VolumeCalotteij(spheres[Monoi[j]]);    //Calcul du volume de la sphérule i moins le volume de
                                                                                                 //la calotte due à la sphérule j
                tabSurf[i] = tabSurf[i] - spheres[Monoi[i]].SurfaceCalotteij(spheres[Monoi[j]]); //Calcul de la surface de la sphérule i moins la surface de
                                                                                                 //la calotte due à la sphérule j
            }
        }
        volAgregat = volAgregat + tabVol[i];    //Volume complet de l'agrégat n°id
        surfAgregat = surfAgregat + tabSurf[i]; //Surface libre complète de l'agrégat n°id

        terme = terme + tabVol[i]/pow(spheres[Monoi[i]].r, 3.0);

        Tv = 1 - (3.0/(4.0*nmonoi*PI))*terme;

        for (k = 1; k <= 3; k++)
            PosiGravite[id][k] = PosiGravite[id][k] + spheres[Monoi[i]].pos[k]*tabVol[i]; //Somme des Vi*xi
    }

    Nc = Nc/2;

    if (nmonoi == 1)      cov = 0;
    else        cov = cov/((double)Nc)/2.0; //Calcul du vrai paramètre de recouvrement moyen de l'agrégat n°id

    for (k = 1; k <= 3; k++)    PosiGravite[id][k] = PosiGravite[id][k]/volAgregat; //Centre de masse de l'agrégat n°id

    for (i = 1; i <= nmonoi; i++)
    {
        li = Distance(spheres[Monoi[i]].pos, PosiGravite[id]); //Distance entre les centres de masse de la sphérule i et de l'agrégat n°id

        r = li + spheres[Monoi[i]].r;
        if (r > rmax)
            rmax = r;

        Arg = Arg + tabVol[i]*pow(li, 2.0);
        Brg = Brg + tabVol[i]*pow(spheres[Monoi[i]].r, 2.0);
    }

    delete[] tabVol;
    delete[] tabSurf;

    return sqrt((Arg+3.0/5.0*Brg)/volAgregat); //Rayon de giration de l'agrégat n°id
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

    for (i = 1; i <= N; i++)
    {
        if (spheres[i].Label == Agg)
        {
            rpmoy = rpmoy + spheres[i].r; //Somme des rayons des sphérules de l'agrégat n°Agg
            rpmoy2 = rpmoy2 + pow(spheres[i].r, 2.0);
            rpmoy3 = rpmoy3 + pow(spheres[i].r, 3.0);
        }
    }
    rpmoy = rpmoy/((double)np);   //Calcul du rayon moyen de l'agrégat n°Agg
    rpmoy2 = rpmoy2/((double)np); //Calcul du rayon moyen d'ordre 2 de l'agrégat n°Agg
    rpmoy3 = rpmoy3/((double)np); //Calcul du rayon moyen d'ordre 3 de l'agrégat n°Agg

    dm = ConvertRg2Dm(np,rg);
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
    Aggregate[Agg][9] = 4.0*PI*rpmoy3/3.0; //Volume de l'agrégat sans recouvrement      (Avant c'était Tv : Taux de recouvrement volumique)
    Aggregate[Agg][10] = cov;          //Paramètre de recouvrement
    Aggregate[Agg][11] = 4.0*PI*rpmoy2;  //Surface libre de l'agrégat sans recouvrement       (Avant c'était surfAgregat/volAgregat : Estimation du rapport surface/volume de l'agrégat)
}
//###############################################################################################################################

//############################################# Conditions aux limites périodiques ##############################################
void ReplacePosi(int id)
{
    int i,j,k,nr;

    nr = SelectLabelEgal(id,MonoRep);

    for (i = 1; i <= 3; i++)
    {
        if (PosiGravite[id][i] > L)
        {
            PosiGravite[id][i] = PosiGravite[id][i] - L;
            for (j = 1; j <= nr; j++)
            {
                k = MonoRep[j];
                spheres[k].pos[i] = spheres[k].pos[i] - L;
            }
        }

        if (PosiGravite[id][i] < 0)
        {
            PosiGravite[id][i] = PosiGravite[id][i] + L;
            for (j = 1; j <= nr; j++)
            {
                k = MonoRep[j];
                spheres[k].pos[i] = spheres[k].pos[i] + L;
            }
        }
    }
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

    NAgg--;
}

void MonTri(int n, double arr[], int index[])
{
    int i, j, id;
    double a;

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

int Probabilite(bool trier,double &deltatemps)
{
    double valAlea,max;
    int i,n,nret;
    max = Maxi2D(5,NAgg);

    if (trier)
    {
        for (i=1; i <= NAgg; i++)
            TpT[i] = max/Aggregate[i][5];

        MonTri(NAgg, TpT, IndexPourTri);
        TriCum[1] = TpT[1];

        for (i=2; i <= NAgg; i++)
            TriCum[i] = TriCum[i-1]+TpT[i];
    }

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
    int agg, j, k, knum, memsuperpo;
    double dist, dc, ret;
    Sphere spheredecale;

    ret = 1.0;
    memsuperpo = 0;
    agg = IdPossible[s][1];

    for (j = 1; j <= N; j++)
    {
        if (spheres[j].Label == agg) //Sphérules de l'agrégat cible
        {
            spheredecale.Update(spheres[j].pos[1]+IdPossible[s][2]*L,spheres[j].pos[2]+IdPossible[s][3]*L,spheres[j].pos[3]+IdPossible[s][4]*L,spheres[j].r);
            for (knum = 1; knum <= nmonoi; knum++) //Sphérules de l'agrégat courant
            {
                k=Monoi[knum];
                dist=spheres[k].Intersection(spheredecale, Vectdir, lpm, dc);
                if (dist == -1)
                    memsuperpo = 1;
                else if (dist < ret)
                    ret = dist;
            }
        }
    }

    if (memsuperpo == 1)
        ret = 1;  //S'il y a au moins deux sphérules en superposition, on ignore l'agglomération des deux agrégats

    return ret;
}
//###############################################################################################################################

//########################################## Détermination des contacts inter-agrégats ##########################################
void CalculDistance(int id, double &distmin, int &aggcontact)
{
    double lpm,dist;
    double dc;
    int nmonoi;
    int i,s;
    int npossible;
    int dx,dy,dz;

    lpm = Aggregate[id][4]; // libre parcours moyen de l'agrégat n° id
    npossible = 0;
    aggcontact = 0;
    distmin = 1.0;

    //+++++++++++++++++++++++++++++++++ Détermination des contacts inter-agrégats possibles +++++++++++++++++++++++++++++++++++++++++++++++++++
    s1.Update(PosiGravite[id], Aggregate[id][6]);
    for (dx = -1;dx <= 1; dx++)
    {
        for (dy = -1; dy <= 1; dy++)
        {
            for (dz = -1; dz <= 1; dz++)
            {
                for (i = 1;i <= NAgg; i++)
                {
                    if (i != id)
                    {
                        s2.Update(PosiGravite[i][1]+dx*L, PosiGravite[i][2]+dy*L, PosiGravite[i][3]+dz*L, Aggregate[i][6]); //Envisage les autres objets ...
                        dist = s1.Intersection(s2, Vectdir, lpm, dc);                                                       //... et leurs images par translation

                        if (dist <= lpm)
                        {
                            npossible++;
                            IdPossible[npossible][1] = i;  //Label d'un agrégat potentiellement contactable
                            IdPossible[npossible][2] = dx; //Prise en compte d'un décalage spatial de l'agrégat suivant x
                            IdPossible[npossible][3] = dy; //Prise en compte d'un décalage spatial de l'agrégat suivant y
                            IdPossible[npossible][4] = dz; //Prise en compte d'un décalage spatial de l'agrégat suivant z
                        }
                    }
                }
            }
        }
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (npossible > 0)
    {
        nmonoi = SelectLabelEgal(id, Monoi); //Liste des monomères qui constituent l'agrégat id

        for (s = 1; s <= npossible; s++) //Pour tous les contacts inter-agrégats possibles
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
    int i, nselect, numreject, numstudy;

    err = 0;

    for (i = 1; i <= N; i++)
        if (spheres[i].Label == AggI)
            spheres[i].Translate(Translate);

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

    nselect = SelectLabelEgal(numreject, MonoSel);

    for (i = 1; i <= nselect; i++)
    {
        spheres[MonoSel[i]].Label = numstudy;
    }

    SupprimeLigne(numreject);

    nselect = SelectLabelSuperieur(numreject, MonoSel);

    for (i = 1; i <= nselect; i++)
        spheres[MonoSel[i]].Label--;

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
            d = Distance(spheres[Monoi[j]].pos,spheres[Monoi[i]].pos);
            dsuperpos = (d-(spheres[Monoi[j]].r+spheres[Monoi[i]].r))/(spheres[Monoi[j]].r+spheres[Monoi[i]].r);
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
        spheres[i].Update(spheres[i].pos, spheres[i].r+Asurfgrowth*pow(spheres[i].r, xsurfgrowth-2)*dt);
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
            dist = Distance(spheres[k].pos,spheres[i].pos); // Calcule la distance centre à centre entre le monomère k et tous les autres
            if (dist <= spheres[k].r+spheres[i].r)
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

        masse = Rho*PI*pow(Dp,3.0)/6;
        surface = PI*pow(Dp,2.0);
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
            Aggregate[i][9] = PI*pow(Dp, 3.0)/6; //Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
            Aggregate[i][10] = 0;                   //Coefficient de pénétration (paramètre de recouvrement Cov)
            Aggregate[i][11] = PI*pow(Dp, 2.0); //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)
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
            Aggregate[i][9] = PI*pow(Dp, 3.0);//Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
            Aggregate[i][10] = 0;                   //Coefficient de pénétration (paramètre de recouvrement Cov)
            Aggregate[i][11] = PI*pow(Dp, 2.0); //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)
        }
    }

    NAgg = N;
}

void Fermeture()
{
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
        fprintf(f,"%5d\t %10.5f\t %10.5f\t %10.5f\t %10.5f\n", spheres[i].Label, spheres[i].r*1E9, spheres[i].pos[1]*1E9, spheres[i].pos[2]*1E9, spheres[i].pos[3]*1E9);

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

void MainWindow::Calcul() //Coeur du programme
{
    double deltatemps, distmin, lpm;
    double thetarandom, phirandom;
    int aggcontact, newnumagg, finfichiersuivitempo, finmem = 0;
    //int tmp,superpo;
    int i, j, co, err;
    bool contact;
    time_t t, t0;


    if (ui->ActModulePhysique->isChecked())     ActiveModulephysique = 1;
    else    ActiveModulephysique = 0;

    if (ui->SuiviTempo->isChecked())    ActiveVariationTempo = 1;
    else    ActiveVariationTempo = 0;

    if (ActiveModulephysique)
    {
        sprintf(commentaires, "Le module physique est activé.\n");
        ui->AfficheurRep->append(commentaires);
    }
    else
    {
        sprintf(commentaires, "Le module physique n'est pas activé.\n");
        ui->AfficheurRep->append(commentaires);
    }

    if (ActiveVariationTempo)
    {
        LectureSuiviTempo();
        sprintf(commentaires, "Le fichier de données de suivi temporel est lu.\n");
        ui->AfficheurRep->append(commentaires);
    }
    else
    {
        sprintf(commentaires, "Le fichier de données sélectionné est le fichier 'params.txt'.\n");
        ui->AfficheurRep->append(commentaires);
    }

    Init();
    co=0;
    //superpo=0;
    NSauve=0;
    temps=0;
    time(&t0);
    contact=true;

    printf("\n");

    while (NAgg > N*5/1000) //Pour N=1000 le calcul s'arrête quand il reste 5 agrégats
    {                       //Pour N=2000 le calcul s'arrête quand il reste 10 agrégats
        qApp->processEvents(); //Permet de rafraichir la fenêtre Qt
        time(&t);
        secondes = t-t0;

        thetarandom = Random()*PI*2;
        phirandom = acos(1-2*Random());
        Vectdir[1] = sin(phirandom)*cos(thetarandom);
        Vectdir[2] = sin(phirandom)*sin(thetarandom);
        Vectdir[3] = cos(phirandom);

        if (ActiveModulephysique)
        {
            NumAgg = Probabilite(contact, deltatemps);
            temps = temps + deltatemps;

            if (ActiveVariationTempo)
            {
                finfichiersuivitempo = rechercheValTab();
                if (finmem != finfichiersuivitempo)
                {
                    sprintf(commentaires, "Attention le suivi temporel est plus long que celui du fichier lu.\n");
                    ui->AfficheurRep->append(commentaires);
                }
                finmem = finfichiersuivitempo;
            }

            CroissanceSurface(deltatemps);

            for (i = 1; i<= NAgg; i++)
                ParametresAgg(i);

            CalculDistance(NumAgg, distmin, aggcontact);
            lpm = Aggregate[NumAgg][4];
            contact = (aggcontact != 0);
        }
        else
        {
            NumAgg = int(Random()*(double)NAgg)+1;
            deltatemps = 0.0;
            temps = temps + 1E-9;
            CalculDistance(NumAgg, distmin, aggcontact);
            lpm = Dpm*1E-9; //On fixe le lpm au Dpm
            contact = (aggcontact != 0);
        }

        if (contact)
        {
            for (i = 1; i <= 3; i++)
                Translate[i] = Vectdir[i]*distmin;

            for (j = 1; j <= N; j++)
                if (spheres[j].Label == IdPossible[aggcontact][1])
                {
                    spheres[j].Update(spheres[j].pos[1]+IdPossible[aggcontact][2]*L,spheres[j].pos[2]+IdPossible[aggcontact][3]*L,spheres[j].pos[3]+IdPossible[aggcontact][4]*L,spheres[j].r);
                }

            newnumagg = Reunit(NumAgg, IdPossible[aggcontact][1], err);
            ParametresAgg(newnumagg);
            temps = temps-deltatemps*(1-distmin/lpm);
            /*
            //+++++ Test de superposition des monomères dans un agrégat lors d'un contact +++++
            tmp = CalculSuperposition(newnumagg);
            if (tmp == 1)     {superpo++;}
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            */
             printf("NAgg=%d  temps=%5.1f E-6 s     CPU=%d sec\t\n", NAgg, temps*1E6, secondes);

             ui->progressBar->setValue(N-NAgg+1);
        }
        else
        {
            for (i = 1; i <= 3; i++)
                Translate[i] = Vectdir[i]*lpm;
            for (i = 1; i <= N; i++)
                if (spheres[i].Label == NumAgg)
                    spheres[i].Translate(Translate);
            for (j = 1; j <= 3; j++)
                PosiGravite[NumAgg][j] = PosiGravite[NumAgg][j] + Translate[j];

            newnumagg = NumAgg;
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

    //printf("Nombre total de superpositions : %d\n",superpo);
    /*
    cout << "L=" << L*1E9
         <<"     lambda=" << lambda*1E9
         <<"     Dpeqmass=" << Dpeqmass
         <<"     rpeqmass=" << rpeqmass*1E9
         <<"     gamma=" << gamma <<endl;
    */
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
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        N=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        T=atof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        Dpm=atof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        sigmaDpm=atof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        FV=atof(commentaires)*1E-6;
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        P=atof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        Rho=atof(commentaires);
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
    coeffB = ui->doubleSpinBoxCoeffCroissance->value();
    xsurfgrowth = ui->doubleSpinBoxPuissancex->value();
    dfe = ui->doubleSpinBoxDimFractale->value();
    kfe = ui->doubleSpinBoxPrefFractal->value();
    sprintf(commentaires, "\nDimension fractale : %1.2f\nPréfacteur fractal : %1.2f\nB = %1.2f\nx = %1.2f\n", dfe, kfe, coeffB, xsurfgrowth);
    ui->AfficheurRep->append(commentaires);

    Asurfgrowth = coeffB*1E-3;
    sprintf(commentaires, "Coefficient de croissance de surface : %e\n", Asurfgrowth);
    ui->AfficheurRep->append(commentaires);

    InitRandom();

    ui->progressBar->setMaximum(N);
    ui->progressBar->setVisible(true);
    ui->progressBar->setValue(0);
    Calcul();
    Fermeture();
    ui->progressBar->setVisible(false);
    ui->AfficheurRep->append(CheminSauve);
    sprintf(commentaires,"\nFin du calcul  ...\n");
    ui->AfficheurRep->append(commentaires);
}
