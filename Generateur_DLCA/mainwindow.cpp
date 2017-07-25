#include "mainwindow.h"
#include <Sphere.h>
#include <aggregat.h>
#include <physical_model.h>
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string.h>
#include <qdir.h>
#include <QString>
#include <stdio.h>
#include <string>
#include <list>
#include <cmath>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

MainWindow* GUI;
Verlet verlet; // Verlet's List
ListSphere spheres;
PhysicalModel physicalmodel;

const double PI = atan(1.0)*4; // 3.14159
double* Vectdir; // Direction of the translation of an aggregate
int NumAgg; // Number of the aggregate in translation
int secondes; // CPU Time variable
int NSauve; // last file number


// Nagg
int NAgg; // Nombre d'aggrégats (1 à l'initialisation)
Aggregate *Aggregates; // Position of the center of gravity
double** OldAggregates; // Array of the different aggregates and their caracteristics
int** AggLabels; // 2 dimensionnal array. It stocks in its first dimensions the Ids of the different aggregates, The second dimensions stocks at the index 0, the number of spheres in the
                // Aggregate, and then , the labels of saif spheres
double* TriCum; // Cumulated probabilities
vector<int> IndexPourTri;
double RayonAggMax=0.; // Value of the radius of the bigger aggregate in the box


// Collision Agg
int** IdPossible;  // Array of the Ids of the aggregates that could be in contact with the one moving


// suivit tempo
double** tab; // generic array
int iValTab=0;
int nb_line;

QString FichierParam = "___";
QString pathParam = "___";
QString FichierSuiviTempo = "___";

char CheminSauve[500];
char commentaires[500];


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
    double m = OldAggregates[1][colonne];

    for (i = 2;i <= nmax; i++)
        if (OldAggregates[i][colonne] > m)    m = OldAggregates[i][colonne];

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

//############################################# Conditions aux limites périodiques ##############################################

void ReplacePosi(int id)
{// This function will relocate an aggregate when it gets out of the box limiting the space
    array<double, 4> trans;
    bool move=false;

    const array<double, 4> pos = Aggregates[id].Position();

    //$ for every dimension
    for (int i = 1; i <= 3; i++)
    {
        //$ Check if it is getting out
        if (pos[i] > physicalmodel.L)
        {
            trans[i] = - physicalmodel.L;
            move = true;
        }
        else if (pos[i] < 0)
        {
            trans[i] = physicalmodel.L;
            move = true;
        }
        else
        {
            trans[i] = 0;
        }
    }

    //$ If it is getting out
    if (move)
    {
        //$ Update the position of aggregate
        Aggregates[id].Translate(trans);
    }
}
//###############################################################################################################################



//############# Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ##############

double RayonGiration(int id, double &rmax, double &Tv, int &Nc, double &cov, double &volAgregat, double &surfAgregat)
{// This function determines the Gyration Radius of the Aggregate Id.
    double dist, rpmoy, dbordabord, li, r, Arg, Brg, terme;
    int i, j, k, nmonoi;
    double * tabVol;
    double* tabSurf;

    cov = 0.0;// Coeficient of mean covering
    Nc = 0; // Number of contacts
    volAgregat = surfAgregat = terme = 0.0; // Volume and surface of Agg Id
    rmax = 0.0; // Maximum radius of the aggregate, this corresponds to the distance between the center of mass of the aggregate and the edge of the furthest ball from said center.
                // It is used to assimilate the aggregate to a sphere when checking for intersections
    Arg = Brg = 0.0; // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient, they are used  used in the final formula of the Radius of Gyration
    Tv = 0.0;

    //$ Identification of the spheres in Agg id

    //$ Determination of the monomeres in Agg Id
    ListSphere mySpheres = ListSphere(spheres, AggLabels[id]);
    nmonoi = mySpheres.size();

    tabVol = new double[nmonoi+1];
    tabSurf = new double[nmonoi+1];

    //$ Initialisation of the arrays of volume, surface of each sphere, and the center of mass
    array<double, 4> newpos={0.,0.,0.,0.};

    for (i = 1; i <= nmonoi; i++)
    {
        tabVol[i] =0.;
        tabSurf[i] =0.;
    }

    //$ For the Spheres i in Agg Id
    for (i = 1; i <= nmonoi; i++)
    {
        //$ Calculation of the volume and surface of monomere i of Agg id
        tabVol[i] += mySpheres[i].Volume(); //Calculation of the volume of i
        tabSurf[i] += mySpheres[i].Surface();    //Calculation of the surface of i

        for (j = i+1; j <= nmonoi; j++) //for the j spheres composing Aggregate n°id
        {

            double voli, volj, surfi, surfj;
            voli = volj = surfi = surfj = 0.;

            //$ Calculation of the intersection between the spheres i and j
            dist = mySpheres[i].Intersection(mySpheres[j],voli,volj,surfi,surfj);

            rpmoy = (mySpheres[i].Radius()+mySpheres[j].Radius())/2.0; //Mean Radius between i and j monomeres
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

        terme = terme + tabVol[i]/mySpheres[i].Volume();

        //$ Calculation of the position of the center of mass
        const array<double, 4> pos = mySpheres[i].Position();
        for (k = 1; k <= 3; k++)
        {
            newpos[k] += pos[k]*tabVol[i];
        }
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
        cov = cov/(double(Nc))/2.0;
    }
    //$ Filling of PosiGravite

    for (k = 1; k <= 3; k++)
    {
        newpos[k] = newpos[k]/volAgregat;
    } //Centre of mass of Agg Id


    Aggregates[id].Position(newpos);

    ReplacePosi(id);

    //$ Determination of the maximal radius of Agg Id and the Radius of gyration

    for (i = 1; i <= nmonoi; i++)
    {
        //$ Determination of the distance between each monomere and the center of mass of Agg Id
        li = mySpheres[i].Distance(Aggregates[id].Position()); //Distance entre les centres de masse de la sphérule i et de l'agrégat n°id

        r = li + mySpheres[i].Radius();

        //$ Calculation of rmax
        rmax=MAX(rmax,r);

        //$ Calculation of Rg
        Arg = Arg + tabVol[i]*pow(li, 2);
        Brg = Brg + tabVol[i]*pow(mySpheres[i].Radius(), 2);
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
{// This function will update the parameter of Agg
    int i, Nc;
    double masse, dm, diff, vit, lpm, Tv, cov, rmax, rg, rpmoy, rpmoy2, rpmoy3;
    double volAgregat, surfAgregat;

    rpmoy = rpmoy2 = rpmoy3 = 0.0;
    //$ Determination of the Radius of gyration of Agg using RayonGiration()
    rg = RayonGiration(Agg, rmax, Tv, Nc, cov, volAgregat, surfAgregat);

    masse = physicalmodel.Rho*volAgregat;//Determination of the real mass of Agg
    //$Determination of the spheres in Agg
    ListSphere MonoSel(spheres, AggLabels[Agg]);
    int np = MonoSel.size();

    //$ Determination of the mean radius of degrees 1,2,3 of Agg
    for (i = 1; i <= np; i++)
    {
        rpmoy = rpmoy + MonoSel[i].Radius(); //Sum of the radius of each sphere in Agg
        rpmoy2 = rpmoy2 + pow(MonoSel[i].Radius(), 2);
        rpmoy3 = rpmoy3 + pow(MonoSel[i].Radius(), 3);
    }


    rpmoy = rpmoy/(double(np));
    rpmoy2 = rpmoy2/(double(np));
    rpmoy3 = rpmoy3/(double(np));


    //$ Determination of Dm using ConvertRg2Dm
    dm = physicalmodel.ConvertRg2Dm(np,rg,rpmoy);

    //$Calculation fo it's mean free path, speed, and diff

    diff = physicalmodel.diffusivity(dm);
    vit =  physicalmodel.velocity(masse);
    lpm = 8*diff/PI/vit;


    //$ Update of Aggregate[] with the parameters
    OldAggregates[Agg][0] = rg; //Gyration Radius
    OldAggregates[Agg][1] = np; //Number of spheres in Agg
    OldAggregates[Agg][2] = Nc; //Cumber of contacts
    OldAggregates[Agg][3] = dm; //Mobility Diameter

    if (physicalmodel.ActiveModulephysique == 1)
    {
        OldAggregates[Agg][4] = lpm;     //Mean Free Path
        OldAggregates[Agg][5] = lpm/vit; //Displacement duration
    }
    else
    {
        OldAggregates[Agg][4] = physicalmodel.Dpm*1E-9;
        OldAggregates[Agg][5] = 1E-6;
    }
    if(rmax>RayonAggMax)
    {RayonAggMax=rmax;}
    OldAggregates[Agg][6] = rmax;          //Radius of the sphere containing Agg
    OldAggregates[Agg][7] = volAgregat;    //Etimation of the Aggregate's volume
    OldAggregates[Agg][8] = surfAgregat;   //Estimation of the sufrace of the aggregate
    OldAggregates[Agg][9] = np*4.0*PI*rpmoy3/3.0; //Volume of the aggregate without considering the spheres covering each other (Avant c'était Tv : Taux de recouvrement volumique)
    OldAggregates[Agg][10] = cov;          //Covering Parameter
    OldAggregates[Agg][11] = np*4.0*PI*rpmoy2;  //Free surface of the aggregate (without covering)(Avant c'était surfAgregat/volAgregat : Estimation du rapport surface/volume de l'agrégat)
}


void SupprimeLigne(int ligne)
{// This functions deletes a line in the arrays Aggregates Agglabels, it is called in Reunit(), when 2 aggregates are in contact and merge into one aggregate
    int i, j;
    //printf("SupLigne  : ");

    for (i = ligne + 1; i<= NAgg; i++)
    {

        for (j = 0; j <= 11; j++)
            OldAggregates[i-1][j] = OldAggregates[i][j];
        const array<double, 4> newpos = Aggregates[i].Position();

        Aggregates[i-1].Position(newpos);

    }
    verlet.Remove(NAgg,Aggregates[NAgg].VerletIndex());

    for (i=ligne+1;i<=NAgg;i++)
    {
        delete[] AggLabels[i-1]; // Considering the 2nd dimension of Agglabels doesn't always have the same size, we have to delete
        AggLabels[i-1]=new int[AggLabels[i][0]+1];//and reallocate each of these sub-arrays
        for (j=0;j<=AggLabels[i][0];j++)
        {
            AggLabels[i-1][j]=AggLabels[i][j];
        }
    }
    NAgg--;


    //$ Update of the labels of every sphere that is in an aggregate indexed highger than the one absorbed
    ListSphere SpheresToReLabel(spheres, AggLabels,ligne,NAgg);
    int nselect = SpheresToReLabel.size();

    for (i = 1; i <= nselect; i++)
        SpheresToReLabel[i].DecreaseLabel();

    for (i=ligne;i<=NAgg;i++)
        Aggregates[i].UpdatesSpheres(spheres, AggLabels[i]);
}

void InsertionSort(int n, double arr[], int index[])
{

    for (int i = 1; i <= n; i++)
        index[i] = i;
    for (int j = 2; j <= n; j++)
    {
        double a = arr[j];
        int id = index[j];
        int i = j-1;

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



void quickSort(double arr[], vector<int>& index, int left, int right) {

      int i = left, j = right;
      double pivot = arr[(left + right) / 2];

      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  double dtmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = dtmp;
                  int itmp = index[i];
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

void quickSort(int n, double arr[], vector<int>& index)
{
      quickSort(arr, index, 1, n);
}

void MonTri(int n, double arr[], vector<int>& index)
{
    for (int i = 1; i <= n; i++)
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
        double TpT[physicalmodel.N+1];

    //$ Sort the timesteps
        for (i=1; i <= NAgg; i++)
            TpT[i] = max/OldAggregates[i][5];

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
double Distance_Aggregate(int id, int s, double lpm)
{
    int agg, j, knum;
    double dist, ret;

    ret = 1.0;
    agg = IdPossible[s][1];

    //Liste des sphérules constituant l'agrégat n°Agg

    ListSphere mySpheres(spheres, AggLabels[id]);
    int nmonoi = mySpheres.size();
    ListSphere MonoSel(spheres, AggLabels[agg]);
    int np = MonoSel.size();

    double* trans = new double[4];
    trans[1] = IdPossible[s][2]*physicalmodel.L;
    trans[2] = IdPossible[s][3]*physicalmodel.L;
    trans[3] = IdPossible[s][4]*physicalmodel.L;

    //$ Loop on all the spheres of the other aggregate
    for (j = 1; j <= np; j++)
    {
        //$ spheredecale is used to replace the sphere into the corresponding box
        Sphere spheredecale(MonoSel[j]);
        spheredecale.Translate(trans);

        //$ For every sphere in the aggregate :
        for (knum = 1; knum <= nmonoi; knum++)
        {
            //$ Check if j and k are contact and if not, what is the distance between them
            dist=mySpheres[knum].Collision(spheredecale, Vectdir, lpm);
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
{   _List_iterator<int> p;
    double lpm,dist;
    double dc,tampon;
    int nmonoi;
    int bornei1,bornei2,bornej1,bornej2,bornek1,bornek2;
    int i,s,j,k;
    int npossible;
    int dx,dy,dz;
    nmonoi =0;
    lpm = OldAggregates[id][4]; // mean free path of the agregate labeled id
    npossible = 0;
    aggcontact = 0;
    distmin = 1.0;

    //+++++++++++++++++++++++++++++++++ Determination of the potential contacts (= in this part, we're considering aggreghates as spheres  with a diameter of Aggregate[6],+++++++++++++++++
    //+++++++++++++++++++++++++++++++++ the distance between the center of gravity of the aggregate and the furthest sphere in the agregate +++++++++++++++++++++++++++++++++++++++++++++++++++

    //$ Find potential contacts between agregates

    Sphere s1(physicalmodel, Aggregates[id].Position(), OldAggregates[id][6]); // Represents the sphere containing the agregate we're testing
    //$ [3 imbricated loops on dx,dy,dz to look into the 27 boxes]

    if (! physicalmodel.use_verlet)
    {
        //$ [For all other agregate]
        for (i = 1;i <= NAgg; i++)
        {
            if (i != id)
            {
                double inix,iniy,iniz,inir;
                const array<double, 4> pos = Aggregates[i].Position();
                inix = pos[1];
                iniy = pos[2];
                iniz = pos[3];
                inir = OldAggregates[i][6]; //represents the different agregates

                //$ [3 imbricated loops on dx,dy,dz to look into the 27 boxes]
                for (dx = -1;dx <= 1; dx++)
                {
                    for (dy = -1; dy <= 1; dy++)
                    {
                        for (dz = -1; dz <= 1; dz++)
                        {
                            Sphere s2(physicalmodel,inix+physicalmodel.L*dx,iniy+physicalmodel.L*dy,iniz+physicalmodel.L*dz,inir);

                            // checks if the two spheres will be in contact while
                             //... the first one is moving
                            //$ Intersection check between agregates
                            dist = s1.Collision(s2, Vectdir, lpm);

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
    }
    else
    {


        double mindist = (OldAggregates[id][6]+RayonAggMax);
        const array<double, 4> pos = Aggregates[id].Position();

        // Détermination des bornes
        if(Vectdir[1]>=0)
        {
            tampon=pos[1]-mindist;
            bornei1=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+1;
            tampon=pos[1]+mindist+lpm*Vectdir[1];
            bornei2=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+2;
        }
        else
        {
            tampon=pos[1]-mindist+lpm*Vectdir[1];
            bornei1=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+1;
            tampon=pos[1]+mindist;
            bornei2=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+2;
        }

        if(Vectdir[2]>=0)
        {
            tampon=pos[2]-mindist;
            bornej1=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+1;
            tampon=pos[2]+mindist+lpm*Vectdir[2];
            bornej2=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+2;
        }
        else
        {
            tampon=pos[2]-mindist+lpm*Vectdir[2];
            bornej1=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+1;
            tampon=pos[2]+mindist;
            bornej2=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+2;
        }

        if(Vectdir[3]>=0)
        {
            tampon=pos[3]-mindist;
            bornek1=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+1;
            tampon=pos[3]+mindist+lpm*Vectdir[3];
            bornek2=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+2;
        }
        else
        {
            tampon=pos[3]-mindist+lpm*Vectdir[3];
            bornek1=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+1;
            tampon=pos[3]+mindist;
            bornek2=floor(tampon*physicalmodel.GridDiv/physicalmodel.L)+2;
        }

        bornei1 = fmax(bornei1,-physicalmodel.GridDiv+1) ; bornei2 = fmin(bornei2,2*physicalmodel.GridDiv);
        bornej1 = fmax(bornej1,-physicalmodel.GridDiv+1) ; bornej2 = fmin(bornej2,2*physicalmodel.GridDiv);
        bornek1 = fmax(bornek1,-physicalmodel.GridDiv+1) ; bornek2 = fmin(bornek2,2*physicalmodel.GridDiv);

        // ///////
        for (i=bornei1+physicalmodel.GridDiv;i<=bornei2+physicalmodel.GridDiv;i++)
        {
            for (j=bornej1+physicalmodel.GridDiv;j<=bornej2+physicalmodel.GridDiv;j++)
            {
                for (k=bornek1+physicalmodel.GridDiv;k<=bornek2+physicalmodel.GridDiv;k++)
                {

                    dx=floor((i-1)/physicalmodel.GridDiv)-1;
                    dy=floor((j-1)/physicalmodel.GridDiv)-1;
                    dz=floor((k-1)/physicalmodel.GridDiv)-1;

                    list<int>* cell = verlet.GetCell(i-dx*physicalmodel.GridDiv,
                                                    j-dy*physicalmodel.GridDiv,
                                                    k-dz*physicalmodel.GridDiv);

                    for(p=cell->begin();p!=cell->end();p++)
                    {
                        if (*p != id)
                        {
                            double inix,iniy,iniz,inir;
                            const array<double, 4> pos = Aggregates[*p].Position();
                            inix = pos[1];
                            iniy = pos[2];
                            iniz = pos[3];
                            inir = OldAggregates[*p][6]; //represents the different agregates

                            Sphere s2(physicalmodel,inix+physicalmodel.L*dx,iniy+physicalmodel.L*dy,iniz+physicalmodel.L*dz,inir);

                            dist = s1.Collision(s2, Vectdir, lpm);
                            if (dist <= lpm)
                            {
                                npossible++; // Number of aggregates that could be hit
                                IdPossible[npossible][1] = *p;  //Label of an aggregate that could be in contact with the one moving
                                IdPossible[npossible][2] = dx; //X coordinate of the "box" where this agregate was
                                IdPossible[npossible][3] = dy; //Y coordinate of the "box" where this agregate was
                                IdPossible[npossible][4] = dz; //Z coordinate of the "box" where this agregate was
                            }
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
        //$ loop on the agregates potentially in contact
        for (s = 1; s <= npossible; s++) //For every aggregate that could be in contact
        {
            dist = Distance_Aggregate(id, s, lpm);
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
{// This function will merge the aggregates AggI and AggJ
    int i, numreject, numstudy;
    int* TamponValeurs; // Buffer variable

    err = 0;

    //$ Check wich one has the smallest index
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


    //$ Creation of the new sub array in Agglabels
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

    //$ Reallocation of the subarray that will contain the labels of the spheres in the reunited aggregate
    delete[] AggLabels[numstudy];
    AggLabels[numstudy] = new int [TamponValeurs[0]+1];
    for (i=0;i<=TamponValeurs[0];i++)
    {
        AggLabels[numstudy][i]=TamponValeurs[i];
    }

    ListSphere SpheresToDelete(spheres, AggLabels[numreject]);
    int nselect = SpheresToDelete.size();
    //$ Update of the labels of the spheres that were in the deleted aggregate
    for (i = 1; i <= nselect; i++)
    {
        SpheresToDelete[i].SetLabel(numstudy);
    }

    //$ Deletionn of the aggregate that was absorbed, using SupprimeLigne()

    SupprimeLigne(numreject);
    delete[] TamponValeurs;

    Aggregates[numstudy].UpdatesSpheres(spheres, AggLabels[numstudy]);

    //$ Index of the Reunited aggregate is returned
    return numstudy;
}
//###############################################################################################################################

//################################### Test de superposition des sphérules dans un agrégat #######################################
int CalculSuperposition(int id)
{
    double d, dsuperpos, mem;
    int i, j, nmonoi;

    ListSphere mySpheres = ListSphere(spheres, AggLabels[id]);
    nmonoi = mySpheres.size();
    mem = 0;

    for (i = 1;i <= nmonoi; i++)
    {
        for (j = 1;j <= i-1; j++)
        {
            d = mySpheres[j].Distance(mySpheres[i]);
            double rj = mySpheres[j].Radius();
            double ri = mySpheres[i].Radius();
            dsuperpos = (d-(rj+ri))/(rj+ri);
            if (dsuperpos < mem)
            {
                mem = dsuperpos;
            }
        }
    }
    mem = mem*1E6;

    if (mem < 0)
    {
        printf("Superposition pour Numfichier %d,   agrégat n°%d    mem=%d\n", NSauve, id, int(mem));
        return 1;
    }
    else return 0;
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

    while (tab[iValTab][0] < physicalmodel.temps && iValTab < (nb_line-1))
    {
        iValTab++;
//        cout << "temps residence = " << temps << "s"
//             << "       temps fichier = " << tab[iValTab][0] << "s"
//             << "    T = " << tab[iValTab][1] << "K" << endl;
    }

    physicalmodel.T = tab[iValTab-1][1];

    if (iValTab == nb_line-1)     test = 1;

    return test;
}
//###############################################################################################################################

void Init()
{
    int i, j, k, test, testmem;
    double x, masse, surface, Dp, Diff, Vit, lpm, rg, dist;

    testmem = 0;
    NAgg = physicalmodel.N;
    spheres.Init(physicalmodel, NAgg);
    Vectdir = new double[4];
    TriCum = new double[NAgg+1];
    IdPossible = new int* [NAgg+1];
    Aggregates = new Aggregate[NAgg+1];
    OldAggregates = new double* [NAgg+1];
    IndexPourTri.assign(NAgg+1,0);


    verlet.Init(physicalmodel.GridDiv);

    // Agglabels
    AggLabels= new int*[NAgg+1];// Array containing the labels of the spheres in each aggregate

    for (i=1;i<=NAgg;i++)
    {                               //_____
        AggLabels[i]= new int[2];   //     |
        AggLabels[i][0]=1;          //     |--- Initialisation of Agglabels, at the start there are N aggregates of 1 sphere.
        AggLabels[i][1]=i;          //_____|
                                    //
    }
    // Agglabels[0] isn't used
    AggLabels[0]=new int[2];
    AggLabels[0][0]=0;
    AggLabels[0][1]=0;

    for (i = 1; i <= NAgg; i++)
    {
        OldAggregates[i] = new double[12];
        IdPossible[i] = new int[5];
    }

    for (i = 1; i <= NAgg; i++)
    {          

        array<double, 4> newpos;
        //random position
        for (j = 1; j<= 3 ; j++)
            newpos[j] = Random()*physicalmodel.L;

        //random size
        x = Random();

        if (physicalmodel.Mode == 1)
            Dp = physicalmodel.Dpm+sqrt(2.0)*physicalmodel.sigmaDpm*inverf(2*x-1); //Loi normale
        else
            Dp = exp(log(physicalmodel.Dpm)+sqrt(2.0)*log(physicalmodel.sigmaDpm)*inverf(2*x-1)); //Loi log-normale

        Dp = Dp/1E9;
        if (Dp <= 0)  Dp = physicalmodel.Dpm*1E-9;

        //++++++++++++ Test de superposition des sphérules lors de leur génération aléatoire ++++++++++++
        test=0;
        for (k = 1; k <= i-1; k++)
        {
            dist = spheres[k].Distance(newpos); // Calcule la distance centre à centre entre le monomère k et tous les autres
            if (dist <= spheres[k].Radius()+Dp/2)
                test++;
        }
        testmem = testmem + test; //Comptabilise le nombre d'échecs à positionner une sphère sans superposition

        if (test > 0)
            i--;
        else
        {
            Aggregates[i].Init(physicalmodel,verlet,newpos,i,spheres,Dp/2);
        }

        if (testmem > NAgg)
        {
            printf("Impossible de générer tous les monomères sans superposition.\n");
            printf("La fraction volumique doit être diminuée.\n");
            exit(0);
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        masse = physicalmodel.Rho*PI*pow(Dp,3)/6;
        surface = PI*pow(Dp,2);
        Diff = physicalmodel.diffusivity(Dp);
        Vit =  physicalmodel.velocity(masse);
        lpm = 8*Diff/PI/Vit;
        rg = sqrt(3.0/5.0)*Dp/2;

        if(Dp/2>RayonAggMax)
        {
            RayonAggMax=Dp/2;
        }
        if (physicalmodel.ActiveModulephysique==1)
        {
            OldAggregates[i][0] = rg;                   //Rayon de gyration
            OldAggregates[i][1] = 1;                    //Nommbre de sphérules par agrégat Np
            OldAggregates[i][2] = 0;                    //Nombre de coordination Nc
            OldAggregates[i][3] = Dp;                   //Diamètre de mobilité
            OldAggregates[i][4] = lpm;                  //Libre parcours moyen
            OldAggregates[i][5] = lpm/Vit;              //Durée du déplacement
            OldAggregates[i][6] = Dp/2;                 //Rayon de la sphère d'enveloppe de l'agrégat réunifié
            OldAggregates[i][7] = masse/physicalmodel.Rho;            //Volume estimé de l'agrégat réunifié
            OldAggregates[i][8] = surface;              //Surface estimée de l'agrégat réunifié
            OldAggregates[i][9] = PI*pow(Dp, 3)/6;      //Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
            OldAggregates[i][10] = 0;                   //Coefficient de pénétration (paramètre de recouvrement Cov)
            OldAggregates[i][11] = PI*pow(Dp, 2);       //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)
        }
        else
        {
            OldAggregates[i][0] = rg;                   //Rayon de gyration
            OldAggregates[i][1] = 1;                    //Nommbre de sphérules par agrégat Np
            OldAggregates[i][2] = 0;                    //Nombre de coordination Nc
            OldAggregates[i][3] = Dp;                   //Diamètre de mobilité
            OldAggregates[i][4] = physicalmodel.Dpm*1E-9;             //Libre parcours moyen
            OldAggregates[i][5] = 1E-6;                 //Durée du déplacement
            OldAggregates[i][6] = Dp/2;                 //Rayon de la sphère d'enveloppe de l'agrégat réunifié
            OldAggregates[i][7] = masse/physicalmodel.Rho;            //Volume de l'agrégat réunifié
            OldAggregates[i][8] = surface;              //Surface de l'agrégat réunifié
            OldAggregates[i][9] = PI*pow(Dp, 3);        //Volume de l'agrégat réunifié sans recouvrement (Avant c'était 0; : Taux de recouvrement volumique)
            OldAggregates[i][10] = 0;                   //Coefficient de pénétration (paramètre de recouvrement Cov)
            OldAggregates[i][11] = PI*pow(Dp, 2);       //Surface de l'agrégat réunifié sans recouvrement (Avant c'était surface/(masse/Rho); : Surface/volume de l'agrégat réunifié)
        }
    }
}

void Fermeture()
{

    for (int i = 1; i <= physicalmodel.N; i++)
    {
        delete[] OldAggregates[i];
        delete[] IdPossible[i];

    }
    delete[] Vectdir;
    delete[] TriCum;
    delete[] IdPossible;
    delete[] Aggregates;
    delete[] OldAggregates;
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
    fprintf(f, "%d  N_[]\n", physicalmodel.N);
    fprintf(f, "%10.3f  FV_[ppm]\n", physicalmodel.FV*1E6);
    fprintf(f, "%10.3f  X_[]\n", physicalmodel.X);
    fprintf(f, "%10.3f  Dpm_[nm]\n", physicalmodel.Dpm);
    fprintf(f, "%10.3f  sigmaDpm_[nm]\n", physicalmodel.sigmaDpm);
    fprintf(f, "%d  NAgg_[]\n", NAgg);
    fprintf(f, "%10.6f  Temps_[µs]\n", physicalmodel.temps*1E6);
    fprintf(f, "Label\t      Rp(nm)\t      X(nm)\t     Y(nm)\t     Z(nm)\n");

    for (i=1;i<=physicalmodel.N;i++)
        fprintf(f,"%s\n", spheres[i].str(1e9).c_str());

    fclose(f);

    sprintf(NomComplet,"%s/Agg%05d.txt", CheminSauve, value);
    f = fopen(NomComplet, "w");
    fprintf(f,"%d  N_[]\n",physicalmodel.N);
    fprintf(f,"%1.3f  FV_[ppm]\n", physicalmodel.FV*1E6);
    fprintf(f,"%1.3f  X_[]\n", physicalmodel.X);
    fprintf(f,"%1.3f  Dpm_[nm]\n", physicalmodel.Dpm);
    fprintf(f,"%1.3f  sigmaDpm_[nm]\n", physicalmodel.sigmaDpm);
    fprintf(f,"%d  NAgg_[]\n", NAgg);
    fprintf(f,"%1.6f  Temps_[µs]\n", physicalmodel.temps*1E6);
    fprintf(f,"Rg_[nm]\tNp_[]\tNc_[]\tDm_[nm]\tlpm_[nm]\tdeltat_[µs]\tRgeo_[nm]\tXG(nm)\tYG(nm)\tZG(nm)\tV_[1E-25m3]\tS_[1E-16m2]\tVOlWO_[1E-25m3]\tcov[]\tSurfWO_[1E-16m2]\n");

    for (i = 1; i <= NAgg; i++)
    {
        fprintf(f, "%10.3f\t", OldAggregates[i][0]*1E9);
        fprintf(f, "%d\t", int(OldAggregates[i][1]));
        fprintf(f, "%d\t", int(OldAggregates[i][2])); //Nombre de coordination
        fprintf(f, "%10.3f\t", OldAggregates[i][3]*1E9);
        fprintf(f, "%10.3f\t", OldAggregates[i][4]*1E9);
        fprintf(f, "%10.3f\t", OldAggregates[i][5]*1E6);
        fprintf(f, "%10.3f\t", OldAggregates[i][6]*1E9);

        const array<double, 4> pos = Aggregates[i].Position();
        for (j = 1; j <= 3; j++)
            fprintf(f, "%10.3f\t", pos[j]*1E9);

        fprintf(f, "%10.3f\t", OldAggregates[i][7]*1E25);  //Estimation du volume d'un agrégat
        fprintf(f, "%10.3f\t", OldAggregates[i][8]*1E16);  //Estimation de la surface d'un agrégat
        fprintf(f, "%10.3f\t", OldAggregates[i][9]*1E25);       //Volume d'un agrégat sans recouvrement (Avant c'était Taux de recouvrement volumique)
        fprintf(f, "%10.3f\t", OldAggregates[i][10]);      //Coefficient de pénétration (paramètre de recouvrement Cov)
        fprintf(f, "%10.3f\n", OldAggregates[i][11]*1E16); //Surface d'un agrégat sans recouvrement (Avant c'était Rapport surface/volume)
    }
    fclose(f);

    sprintf(NomComplet, "%s/AFractalPlot.txt", CheminSauve);
    if (OldAggregates[id][1] > 10)
    {
        f=fopen(NomComplet, "a");
        fprintf(f, "%e\t%e\n", OldAggregates[id][0]/(physicalmodel.Dpm*1E-9), OldAggregates[id][1]);
        fclose(f);
    }
}


bool locale_with_dots()
{
    static bool tested = false;
    static bool with_dots;

    if (tested)
        return with_dots;
    else
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
        return with_dots;
    }
}

double latof(const char* _char)
{
    string mystring = _char;
    if (!locale_with_dots())
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
        physicalmodel.N = 2500;
        physicalmodel.T = 1500;
        physicalmodel.Dpm = 30;
        physicalmodel.sigmaDpm = 0.0;
        physicalmodel.FV = 3E-3;
        physicalmodel.Rho = 1.8E3;
        physicalmodel.P = 101300;
        physicalmodel.Mode = 1;
        physicalmodel.DeltaSauve = 0;
        sprintf(sauve, "%s", "DATA");
    }
    else
    {
        char com[500]; // Char array used in the ASCII Save

        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.N=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.T=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.Dpm=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.sigmaDpm=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.FV=latof(commentaires)*1E-6;
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.P=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.Rho=latof(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.Mode=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",commentaires,com);
        physicalmodel.DeltaSauve=atoi(commentaires);
        fgets(com,500,f);
        sscanf(com,"%s  %s",sauve,com);
        fclose(f);
    }
    if (physicalmodel.Mode == 1)
        physicalmodel.X = pow(physicalmodel.N*PI/6.0/physicalmodel.FV*(1.0+3.0*physicalmodel.sigmaDpm*physicalmodel.sigmaDpm/physicalmodel.Dpm/physicalmodel.Dpm),1.0/3.0); //Loi normale
    else
        physicalmodel.X = pow(physicalmodel.N*PI/6.0/physicalmodel.FV*exp(9.0/2.0*log(physicalmodel.sigmaDpm)*log(physicalmodel.sigmaDpm)),1.0/3.0); //Loi log-normale

    strcat(CheminSauve,qPrintable(pathParam));
    strcat(CheminSauve,"//");
    strcat(CheminSauve,sauve);
    lDir.mkdir(CheminSauve);

    sprintf(commentaires, "N=%d \nT=%1.3f \nDpm=%1.3f \nsigmaDpm=%1.3f \nFV=%1.3e\nX=%1.3f \nP=%1.3f\nMode=%d\nRho=%1.3f \nDeltaSauve=%d\nCheminSauve=%s\n", physicalmodel.N, physicalmodel.T, physicalmodel.Dpm, physicalmodel.sigmaDpm, physicalmodel.FV, physicalmodel.X, physicalmodel.P, physicalmodel.Mode, physicalmodel.Rho, physicalmodel.DeltaSauve, CheminSauve);
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


    physicalmodel.L = physicalmodel.X*physicalmodel.Dpm*1E-9;
    physicalmodel.Init(physicalmodel.P,physicalmodel.T,physicalmodel.dfe,
                       physicalmodel.kfe,physicalmodel.Dpm,physicalmodel.sigmaDpm,
                       physicalmodel.xsurfgrowth,physicalmodel.coeffB,physicalmodel.Rho);


    if (physicalmodel.ActiveModulephysique)
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

    if (physicalmodel.ActiveVariationTempo)
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

    sprintf(commentaires, "\nDimension fractale : %1.2f\nPréfacteur fractal : %1.2f\nB = %1.2f\nx = %1.2f\n",
                          physicalmodel.dfe, physicalmodel.kfe, physicalmodel.coeffB, physicalmodel.xsurfgrowth);
    if (GUI == NULL)
        printf("%s",commentaires);
    else
        GUI->print(commentaires);

    double Asurfgrowth = physicalmodel.coeffB*1E-3;
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
    physicalmodel.temps=0;
    deltatemps=0;
    time(&t0);
    contact=true;
    int end = MAX(5,physicalmodel.N/200);
    int it_without_contact=0;
    int lim_it_without_contact = 200;

    end = 20;

    printf("\n");
    printf("Ending calcul when there is less than %d aggregats or %d iterations without contact\n", end,lim_it_without_contact);

    //$ Loop on the N monomeres
    while (NAgg > end && it_without_contact<lim_it_without_contact) //Pour N=1000 le calcul s'arrête quand il reste 5 agrégats
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

        if (physicalmodel.ActiveModulephysique)
        {
            //$ Choice of an aggregate according to his MFP
            NumAgg = Probabilite(contact, deltatemps);// Choice of an agrgegate, the probability of said agrgegate to be chosen proportionnal to his lpm

            physicalmodel.temps = physicalmodel.temps + deltatemps; // Time incrementation with a value given by Probabilite

            if (physicalmodel.ActiveVariationTempo)
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

            //$ Surface Growth of all spheres
            spheres.CroissanceSurface(deltatemps);

            //$ Aggregates parameter update

            for (i = 1; i<= NAgg; i++)
            {
                ParametresAgg(i);
            }
            lpm = OldAggregates[NumAgg][4];
        }
        else
        {
            //$ Random Choice of an aggregate
            NumAgg = int(Random()*double(NAgg))+1;
            deltatemps = 0.0;
            physicalmodel.temps = physicalmodel.temps + 1E-9;
            lpm = physicalmodel.Dpm*1E-9; //On fixe le lpm au Dpm
        }

        double distmove = lpm;

        //$ looking for potential contacts
        CalculDistance(NumAgg, distmin, aggcontact);

        contact = (aggcontact != 0);
        if (contact)
            distmove = distmin;

        //$ Translation of the aggregate
        for (i = 1; i <= 3; i++)
        {
            Vectdir[i] = Vectdir[i]*distmove;
        }

        Aggregates[NumAgg].Translate(Vectdir);

        if (contact)
        {

            //the aggregate that's been tested in contact with the one moving is replaced in the "box" of the space where the contact happened
            //It gets in box on one side, then has to go out on the other one.

            double* trans = new double[4];
            trans[1] = IdPossible[aggcontact][2]*physicalmodel.L;
            trans[2] = IdPossible[aggcontact][3]*physicalmodel.L;
            trans[3] = IdPossible[aggcontact][4]*physicalmodel.L;

            Aggregates[IdPossible[aggcontact][1]].Translate(trans);

            delete[] trans;

            //$ Aggregates in contact are reunited
            newnumagg = Reunit(NumAgg, IdPossible[aggcontact][1], err);

            ParametresAgg(newnumagg);

            physicalmodel.temps = physicalmodel.temps-deltatemps*(1-distmin/lpm);
            /*
            //+++++ Test de superposition des monomères dans un agrégat lors d'un contact +++++
            tmp = CalculSuperposition(newnumagg);
            if (tmp == 1)     {superpo++;}
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            */
             printf("NAgg=%d  temps=%5.1f E-6 s     CPU=%d sec    after %d it\n", NAgg, physicalmodel.temps*1E6, secondes,it_without_contact);
            it_without_contact = 0;
             if (!(GUI == NULL)            )
                GUI->progress(physicalmodel.N-NAgg+1);

        }
        else
        {
            newnumagg = NumAgg;
            it_without_contact++;
        }
        ReplacePosi(newnumagg);


        if (physicalmodel.DeltaSauve>0)
        {
            co++;
            if (co > physicalmodel.DeltaSauve)
            {
                co = 0;
                if (!contact)
                    SauveASCII(NSauve++, newnumagg);
            }
        }

        if (contact && physicalmodel.DeltaSauve >= 0)
            SauveASCII(NSauve++, newnumagg);
    }

    printf("Nombre total d'aggregats : %d\nNombre d'iterations sans contacts' : %d\n",NAgg,it_without_contact);
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
        const array<double, 4> pos = Aggregates[i].Position();
        for (j = 1; j <= 3; j++)
            printf("%e\t", pos[j]*1E9);
        printf("\t%e\t%e\n",OldAggregates[i][7]*1E25,OldAggregates[i][9]*1E25);
    }
    printf("\n\n");

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

__attribute__((noreturn)) void MainWindow::BoutonQuitter()
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
    QString pathSuiviTempo = tmp3.absolutePath(); //Cette variable ne retient que le chemin du fichier Suivi Tempo

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

    ui->progressBar->setMaximum(physicalmodel.N);
    ui->progressBar->setVisible(true);
    ui->progressBar->setValue(0);

    physicalmodel.coeffB = ui->doubleSpinBoxCoeffCroissance->value();
    physicalmodel.xsurfgrowth = ui->doubleSpinBoxPuissancex->value();
    physicalmodel.dfe = ui->doubleSpinBoxDimFractale->value();
    physicalmodel.kfe = ui->doubleSpinBoxPrefFractal->value();

    if (ui->ActModulePhysique->isChecked())     physicalmodel.ActiveModulephysique = 1;
    else    physicalmodel.ActiveModulephysique = 0;

    if (ui->SuiviTempo->isChecked())    physicalmodel.ActiveVariationTempo = 1;
    else    physicalmodel.ActiveVariationTempo = 0;

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

    physicalmodel.coeffB = 1.0;
    physicalmodel.xsurfgrowth = 2.0;
    physicalmodel.dfe = 1.80;
    physicalmodel.kfe = 1.40;
    physicalmodel.ActiveModulephysique = 1;
    physicalmodel.ActiveVariationTempo = 0;

    Calcul();

    return 0;
}

