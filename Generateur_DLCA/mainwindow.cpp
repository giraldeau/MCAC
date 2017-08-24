#include <libgen.h>
#include <sys/stat.h>
#include "mainwindow.h"
#include <Sphere.h>
#include <aggregat.h>
#include <physical_model.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <string>
#include <list>
#include <cmath>
#include <limits.h> /* PATH_MAX */

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

#ifdef WITH_GUI
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <qdir.h>
#include <QString>

MainWindow* GUI;
#endif

PhysicalModel physicalmodel;

const double PI = atan(1.0)*4; // 3.14159
array<double,4> Vectdir; // Direction of the translation of an aggregate
int NumAgg; // Number of the aggregate in translation
time_t secondes; // CPU Time variable
int NSauve; // last file number


// Nagg
int NAgg; // Nombre d'aggrégats (1 à l'initialisation)
ListAggregat Aggregates; // Position of the center of gravity
int** AggLabels; // 2 dimensionnal array. It stocks in its first dimensions the Ids of the different aggregates, The second dimensions stocks at the index 0, the number of spheres in the
                // Aggregate, and then , the labels of saif spheres
double* TriCum; // Cumulated probabilities
vector<int> IndexPourTri;
double RayonAggMax=0.; // Value of the radius of the bigger aggregate in the box


// Collision Agg
int** IdPossible;  // Array of the Ids of the aggregates that could be in contact with the one moving
int max_npossible;

// suivit tempo
string FichierParam = "___";
string pathParam = "___";
string FichierSuiviTempo = "___";
double** tab; // generic array
int iValTab=0;
int nb_line;


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
    double m = Aggregates[1][colonne];
    for (int i = 1;i < nmax; i++)
        m =max(Aggregates[i][colonne], m);

    return m;
}


void SupprimeLigne(int ligne)
{// This functions deletes a line in the arrays Aggregates Agglabels, it is called in Reunit(), when 2 aggregates are in contact and merge into one aggregate
    int i, j;
    //printf("SupLigne  : ");

    for (i = ligne + 1; i<= NAgg; i++)
    {

        for (j = 0; j <= 11; j++)
            Aggregates[i-1][j] = Aggregates[i][j];
        const array<double, 4> newpos = Aggregates[i].GetPosition();

        Aggregates[i-1].SetPosition(newpos);

    }
    Aggregates.verlet.Remove(NAgg,Aggregates[NAgg].VerletIndex());

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
    ListSphere SpheresToReLabel(Aggregates.spheres, AggLabels,ligne,NAgg);
    SpheresToReLabel.DecreaseLabel();

    for (i=ligne;i<=NAgg;i++)
        Aggregates[i].UpdatesSpheres(Aggregates.spheres, AggLabels[i]);
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
        double* TpT = new double [physicalmodel.N+1];

    //$ Sort the timesteps
        for (i=1; i <= NAgg; i++)
            TpT[i] = max/Aggregates[i][5];

        MonTri(NAgg, TpT, IndexPourTri); //$

    //$ Accumulate the timesteps
        TriCum[1] = TpT[1];
        for (i=2; i <= NAgg; i++)
            TriCum[i] = TriCum[i-1]+TpT[i];
        delete[] TpT;
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

//########################################## Determination of the contacts between agrgates ##########################################
void CalculDistance(int id, double &distmin, int &aggcontact)
{   _List_iterator<int> p;
    double lpm,dist;
    int bornei1,bornei2,bornej1,bornej2,bornek1,bornek2;
    int i,s,j,k;
    int npossible;
    int dx,dy,dz;
    lpm = Aggregates[id][4]; // mean free path of the agregate labeled id
    npossible = 0;
    aggcontact = 0;
    distmin = 1.0;

    //+++++++++++++++++++++++++++++++++ Determination of the potential contacts (= in this part, we're considering aggreghates as spheres  with a diameter of Aggregate[6],+++++++++++++++++
    //+++++++++++++++++++++++++++++++++ the distance between the center of gravity of the aggregate and the furthest sphere in the agregate +++++++++++++++++++++++++++++++++++++++++++++++++++

    //$ Find potential contacts between agregates

    Sphere s1(physicalmodel, Aggregates[id].GetPosition(), Aggregates[id][6]); // Represents the sphere containing the agregate we're testing
    //$ [3 imbricated loops on dx,dy,dz to look into the 27 boxes]

    if (! physicalmodel.use_verlet)
    {
        //$ [For all other agregate]
        for (i = 1;i <= NAgg; i++)
        {
            if (i != id)
            {
                double inix,iniy,iniz,inir;
                const array<double, 4> pos = Aggregates[i].GetPosition();
                inix = pos[1];
                iniy = pos[2];
                iniz = pos[3];
                inir = Aggregates[i][6]; //represents the different agregates

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
                                if (npossible == max_npossible)
                                {
                                    cout << "Too many possible collisions" << endl;
                                    exit(2);
                                }
                                }
                        }
                    }
                }
            }
        }    
    }
    else
    {


        double mindist = (Aggregates[id][6]+RayonAggMax);
        const array<double, 4> posid = Aggregates[id].GetPosition();

        double xp(posid[1]+mindist+MAX(lpm*Vectdir[1],0));
        double xm(posid[1]-mindist+MIN(lpm*Vectdir[1],0));
        double yp(posid[2]+mindist+MAX(lpm*Vectdir[2],0));
        double ym(posid[2]-mindist+MIN(lpm*Vectdir[2],0));
        double zp(posid[3]+mindist+MAX(lpm*Vectdir[3],0));
        double zm(posid[3]-mindist+MIN(lpm*Vectdir[3],0));

        bornei1 = int(floor(xm*physicalmodel.GridDiv/physicalmodel.L));
        bornei2 = int(floor(xp*physicalmodel.GridDiv/physicalmodel.L)+1);
        bornej1 = int(floor(ym*physicalmodel.GridDiv/physicalmodel.L));
        bornej2 = int(floor(yp*physicalmodel.GridDiv/physicalmodel.L)+1);
        bornek1 = int(floor(zm*physicalmodel.GridDiv/physicalmodel.L));
        bornek2 = int(floor(zp*physicalmodel.GridDiv/physicalmodel.L)+1);

        if (bornei2-bornei1>=physicalmodel.GridDiv)
            {bornei1=0 ; bornei2=physicalmodel.GridDiv-1;}
        if (bornej2-bornej1>=physicalmodel.GridDiv)
            {bornej1=0 ; bornej2=physicalmodel.GridDiv-1;}
        if (bornek2-bornek1>=physicalmodel.GridDiv)
            {bornek1=0 ; bornek2=physicalmodel.GridDiv-1;}


        // ///////
        for (i=bornei1;i<=bornei2;i++)
        {
            for (j=bornej1;j<=bornej2;j++)
            {
                for (k=bornek1;k<=bornek2;k++)
                {
                    int ii(i),jj(j),kk(k);

                    // periodic
                    while (ii<0)
                        ii += physicalmodel.GridDiv;
                    while (jj<0)
                        jj += physicalmodel.GridDiv;
                    while (kk<0)
                        kk += physicalmodel.GridDiv;
                    while (ii>=physicalmodel.GridDiv)
                        ii -= physicalmodel.GridDiv;
                    while (jj>=physicalmodel.GridDiv)
                        jj -= physicalmodel.GridDiv;
                    while (kk>=physicalmodel.GridDiv)
                        kk -= physicalmodel.GridDiv;


                    list<int>* cell = Aggregates.verlet.GetCell(ii,jj,kk);

                    for(p=cell->begin();p!=cell->end();p++)
                    {
                        if (*p != id)
                        {
                            double inix,iniy,iniz,inir;
                            const array<double, 4> posp = Aggregates[*p].GetPosition();
                            inix = posp[1];
                            iniy = posp[2];
                            iniz = posp[3];
                            inir = Aggregates[*p][6]; //represents the different agregates

                            Sphere s2(physicalmodel,inix,iniy,iniz,inir);

                            dist = s1.Collision(s2, Vectdir, lpm);
                            if (dist <= lpm)
                            {
                                npossible++; // Number of aggregates that could be hit
                                IdPossible[npossible][1] = *p;  //Label of an aggregate that could be in contact with the one moving
                                IdPossible[npossible][2] = 0; //X coordinate of the "box" where this agregate was
                                IdPossible[npossible][3] = 0; //Y coordinate of the "box" where this agregate was
                                IdPossible[npossible][4] = 0; //Z coordinate of the "box" where this agregate was
                                if (npossible == max_npossible)
                                {
                                    cout << "Too many possible collisions" << endl;
                                    exit(2);
                                }
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
            int agg = IdPossible[s][1];

            array<double,4> trans;
            trans[1] = IdPossible[s][2]*physicalmodel.L;
            trans[2] = IdPossible[s][3]*physicalmodel.L;
            trans[3] = IdPossible[s][4]*physicalmodel.L;

            dist = Aggregates[id].Distance_Aggregate(Aggregates[agg],trans,Vectdir);
            if (dist >= 0 && dist < distmin)
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

    ListSphere SpheresToDelete(Aggregates.spheres, AggLabels[numreject]);
    int nselect = SpheresToDelete.size;
    //$ Update of the labels of the spheres that were in the deleted aggregate
    for (i = 1; i <= nselect; i++)
    {
        SpheresToDelete[i].SetLabel(numstudy);
    }

    //$ Deletionn of the aggregate that was absorbed, using SupprimeLigne()

    SupprimeLigne(numreject);
    delete[] TamponValeurs;

    Aggregates[numstudy].UpdatesSpheres(Aggregates.spheres, AggLabels[numstudy]);

    //$ Index of the Reunited aggregate is returned
    return numstudy;
}
//###############################################################################################################################

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
        print(commentaires);
    }
    else
    {
        sprintf(commentaires, "Le module physique n'est pas activé.\n");
        print(commentaires);
    }

    if (physicalmodel.ActiveVariationTempo)
    {
        LectureSuiviTempo();
        sprintf(commentaires, "Le fichier de données de suivi temporel est lu.\n");
        print(commentaires);
    }
    else
    {
        sprintf(commentaires, "Le fichier de données sélectionné est le fichier 'params.txt'.\n");
        print(commentaires);
    }

    sprintf(commentaires, "\nDimension fractale : %1.2f\nPréfacteur fractal : %1.2f\nB = %1.2f\nx = %1.2f\n",
                          physicalmodel.dfe, physicalmodel.kfe, physicalmodel.coeffB, physicalmodel.xsurfgrowth);
    print(commentaires);

    double Asurfgrowth = physicalmodel.coeffB*1E-3;
    sprintf(commentaires, "Coefficient de croissance de surface : %e\n", Asurfgrowth);
    print(commentaires);

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
#ifdef WITH_GUI
        qApp->processEvents(); //Permet de rafraichir la fenêtre Qt
#endif
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
                    print(commentaires);
                }
                finmem = finfichiersuivitempo;
            }

            //$ Surface Growth of all spheres
            Aggregates.spheres.CroissanceSurface(deltatemps);

            //$ Aggregates parameter update

            for (i = 1; i<= NAgg; i++)
            {
                Aggregates[i].Update();
                double rmax = Aggregates[i][6];
                if(rmax>RayonAggMax)
                    RayonAggMax=rmax;
            }
            lpm = Aggregates[NumAgg][4];
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

            double trans[4];
            trans[1] = IdPossible[aggcontact][2]*physicalmodel.L;
            trans[2] = IdPossible[aggcontact][3]*physicalmodel.L;
            trans[3] = IdPossible[aggcontact][4]*physicalmodel.L;

            Aggregates[IdPossible[aggcontact][1]].Translate(trans);

            //$ Aggregates in contact are reunited
            newnumagg = Reunit(NumAgg, IdPossible[aggcontact][1], err);

            Aggregates[newnumagg].Update();
            double rmax = Aggregates[newnumagg][6];
            if(rmax>RayonAggMax)
                RayonAggMax=rmax;

            physicalmodel.temps = physicalmodel.temps-deltatemps*(1-distmin/lpm);
            /*
            //+++++ Test de superposition des monomères dans un agrégat lors d'un contact +++++
            tmp = CalculSuperposition(newnumagg);
            if (tmp == 1)     {superpo++;}
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            */
             printf("NAgg=%d  temps=%5.1f E-6 s     CPU=%d sec    after %d it\n", NAgg, physicalmodel.temps*1E6, int(secondes),it_without_contact);
            it_without_contact = 0;
#ifdef WITH_GUI
             if (!(GUI == nullptr)            )
                GUI->progress(physicalmodel.N-NAgg+1);
#endif
        }
        else
        {
            newnumagg = NumAgg;
            it_without_contact++;
        }



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
        const array<double, 4> pos = Aggregates[i].GetPosition();
        for (j = 1; j <= 3; j++)
            printf("%e\t", pos[j]*1E9);
        printf("\t%e\t%e\n",Aggregates[i][7]*1E25,Aggregates[i][9]*1E25);
    }
    printf("\n\n");

    Fermeture();
    print(CheminSauve);

    sprintf(commentaires,"\nFin du calcul  ...\n");
    print(commentaires);
}


void Init()
{

    int testmem = 0;
    NAgg = physicalmodel.N;
    max_npossible = NAgg;
    TriCum = new double[NAgg+1];
    IdPossible = new int* [max_npossible+1];
    IndexPourTri.assign(NAgg+1,0);
    Aggregates.Init(physicalmodel, NAgg);

    // Agglabels
    AggLabels= new int*[NAgg+1];// Array containing the labels of the spheres in each aggregate

    for (int i=1;i<=NAgg;i++)
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

    for (int i = 1; i <= max_npossible; i++)
    {
        IdPossible[i] = new int[5];
    }

    for (int i = 1; i <= NAgg; i++)
    {          

        array<double, 4> newpos;
        //random position
        for (int j = 1; j<= 3 ; j++)
            newpos[j] = Random()*physicalmodel.L;

        //random size
        double x = Random();
        double Dp =0;
        if (physicalmodel.Mode == 1)
            Dp = physicalmodel.Dpm+sqrt(2.0)*physicalmodel.sigmaDpm*inverf(2*x-1); //Loi normale
        else
            Dp = exp(log(physicalmodel.Dpm)+sqrt(2.0)*log(physicalmodel.sigmaDpm)*inverf(2*x-1)); //Loi log-normale

        Dp = Dp/1E9;
        if (Dp <= 0)  Dp = physicalmodel.Dpm*1E-9;

        //++++++++++++ Test de superposition des sphérules lors de leur génération aléatoire ++++++++++++
        int test=0;
        for (int k = 1; k <= i-1; k++)
        {
            double dist = Aggregates.spheres[k].Distance(newpos); // Calcule la distance centre à centre entre le monomère k et tous les autres
            if (dist <= Aggregates.spheres[k].Radius()+Dp/2)
                test++;
        }
        testmem = testmem + test; //Comptabilise le nombre d'échecs à positionner une sphère sans superposition

        if (test > 0)
            i--;
        else
        {
            Aggregates[i].Init(physicalmodel,Aggregates.verlet,newpos,i,Aggregates.spheres,Dp);

            if(Dp/2>RayonAggMax)
            {
                RayonAggMax=Dp/2;
            }
        }
        if (testmem > NAgg)
        {
            printf("Impossible de générer tous les monomères sans superposition.\n");
            printf("La fraction volumique doit être diminuée.\n");
            exit(0);
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
}

void Fermeture()
{

    for (int i = 1; i <= physicalmodel.N; i++)
    {
        delete[] IdPossible[i];

    }
    delete[] TriCum;
    delete[] IdPossible;
    delete[] tab;
}



//################################ Lecture des données physiques d'entrée variables dans le temps ###############################
int LectureSuiviTempo()
{
    FILE* f = nullptr;
    f = fopen(FichierSuiviTempo.c_str(), "rt");
    int i, j, k;
    char skip_line[500], data[500];

    nb_line=0;

    if (f != nullptr)
    {
        //On lit la première ligne pour la sauter
        if( fgets(skip_line, 500, f) == nullptr)
        {
            cout << "le fichier d'entrée est vide !" << endl;
            exit(1);
        }

        while(fgets(data, 500, f) != nullptr)
        {
            nb_line++; //On compte le nombre de lignes ...
        }

        printf("\nLe fichier contient %d lignes.\n", nb_line); //... et on l'affiche à l'écran
        printf("\n");

        //On se replace au début du fichier ...
        fseek(f, 0, SEEK_SET);

        //On lit la première ligne pour la sauter
        if( fgets(skip_line, 500, f) == nullptr)
        {
            cout << "le fichier d'entrée est vide !" << endl;
            exit(1);
        }

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
            if (fgets(data, 500, f) != nullptr)
            {
                //On crée une chaine qui contiendra le token
                char *token;

                //On découpe selon les tabulations
                token = strtok(data, "\t");

                while(token != nullptr)
                {
                    for (j = 0; j < 2; j++)
                    {
                        //On enregistre la valeur convertie en double dans le tableau
                        tab[i][j] = atof(token);

                        //On affiche le tableau
                        //printf("ligne %d : tableau[%d][%d] = %.6f\n",i, i, j, tab[i][j]);

                        token = strtok(nullptr, "\t");
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
        fprintf(f,"%s\n", Aggregates.spheres[i].str(1e9).c_str());

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
        fprintf(f, "%10.3f\t", Aggregates[i][0]*1E9);
        fprintf(f, "%d\t", int(Aggregates[i][1]));
        fprintf(f, "%d\t", int(Aggregates[i][2])); //Nombre de coordination
        fprintf(f, "%10.3f\t", Aggregates[i][3]*1E9);
        fprintf(f, "%10.3f\t", Aggregates[i][4]*1E9);
        fprintf(f, "%10.3f\t", Aggregates[i][5]*1E6);
        fprintf(f, "%10.3f\t", Aggregates[i][6]*1E9);

        const array<double, 4> pos = Aggregates[i].GetPosition();
        for (j = 1; j <= 3; j++)
            fprintf(f, "%10.3f\t", pos[j]*1E9);

        fprintf(f, "%10.3f\t", Aggregates[i][7]*1E25);  //Estimation du volume d'un agrégat
        fprintf(f, "%10.3f\t", Aggregates[i][8]*1E16);  //Estimation de la surface d'un agrégat
        fprintf(f, "%10.3f\t", Aggregates[i][9]*1E25);       //Volume d'un agrégat sans recouvrement (Avant c'était Taux de recouvrement volumique)
        fprintf(f, "%10.3f\t", Aggregates[i][10]);      //Coefficient de pénétration (paramètre de recouvrement Cov)
        fprintf(f, "%10.3f\n", Aggregates[i][11]*1E16); //Surface d'un agrégat sans recouvrement (Avant c'était Rapport surface/volume)
    }
    fclose(f);

    sprintf(NomComplet, "%s/AFractalPlot.txt", CheminSauve);
    if (Aggregates[id][1] > 10)
    {
        f=fopen(NomComplet, "a");
        fprintf(f, "%e\t%e\n", Aggregates[id][0]/(physicalmodel.Dpm*1E-9), Aggregates[id][1]);
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
        size_t f = mystring.find(".");
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
    f = fopen(FichierParam.c_str(), "rt");

    if (f == nullptr)
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

        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the N parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.N=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the T parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.T=latof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the Dpm parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.Dpm=latof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the sigmaDpm parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.sigmaDpm=latof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the FV parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.FV=latof(commentaires)*1E-6;
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the P parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.P=latof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the Rho parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.Rho=latof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the Mode parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.Mode=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the DeltaSauve parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.DeltaSauve=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the sauve parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",sauve,com);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the coeffB parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.coeffB=atof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the xsurfgrowth parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.xsurfgrowth=atof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the dfe parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.dfe=atof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the kfe parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.kfe=atof(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the ActiveModulephysique parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.ActiveModulephysique=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the ActiveVariationTempo parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.ActiveVariationTempo=atoi(commentaires);
        }
        fclose(f);
    }
    if (physicalmodel.Mode == 1)
        physicalmodel.X = pow(physicalmodel.N*PI/6.0/physicalmodel.FV*(1.0+3.0*physicalmodel.sigmaDpm*physicalmodel.sigmaDpm/physicalmodel.Dpm/physicalmodel.Dpm),1.0/3.0); //Loi normale
    else
        physicalmodel.X = pow(physicalmodel.N*PI/6.0/physicalmodel.FV*exp(9.0/2.0*log(physicalmodel.sigmaDpm)*log(physicalmodel.sigmaDpm)),1.0/3.0); //Loi log-normale

    strcat(CheminSauve,pathParam.c_str());
    strcat(CheminSauve,"/");
    strcat(CheminSauve,sauve);

    sprintf(commentaires, "N=%d \nT=%1.3f \nDpm=%1.3f \nsigmaDpm=%1.3f \nFV=%1.3e\nX=%1.3f \nP=%1.3f\nMode=%d\nRho=%1.3f \nDeltaSauve=%d\nCheminSauve=%s\n", physicalmodel.N, physicalmodel.T, physicalmodel.Dpm, physicalmodel.sigmaDpm, physicalmodel.FV, physicalmodel.X, physicalmodel.P, physicalmodel.Mode, physicalmodel.Rho, physicalmodel.DeltaSauve, CheminSauve);

    cout <<  commentaires << endl;
    if( !dirExists(CheminSauve))
    {
        const int dir_err = mkdir(CheminSauve, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err)
        {
            printf("Error creating directory\n");
            exit(1);
        }
    }
}



int No_GUI(int argc, char *argv[]){
    FichierParam = argv[1];

    char tmp[PATH_MAX + 1]; /* not sure about the "+ 1" */
    char* err = realpath(dirname(argv[1]), tmp);
    pathParam = tmp; //Cette variable ne retient que le chemin du fichier param

    if (err==nullptr)
    {
        cout << "Param file does not exist\n" << endl;
        exit(1);
    }

    LectureParams();

    Calcul();

    return 0;
}


void print(char* str)
{
#ifdef WITH_GUI
    if (GUI == nullptr)
        printf("%s",str);
    else
        GUI->print(str);
#else
    printf("%s",str);
#endif
}

int dirExists(const char *path)
{
    struct stat info;

    if(stat( path, &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 1;
    else
        return 0;
}



#ifdef WITH_GUI



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

    ui->ActModulePhysique->setChecked(physicalmodel.ActiveModulephysique);
    ui->SuiviTempo->setChecked(physicalmodel.ActiveVariationTempo);
    ui->doubleSpinBoxCoeffCroissance->setValue(physicalmodel.coeffB);
    ui->doubleSpinBoxPuissancex->setValue(physicalmodel.xsurfgrowth);
    ui->doubleSpinBoxDimFractale->setValue(physicalmodel.dfe);
    ui->doubleSpinBoxPrefFractal->setValue(physicalmodel.kfe);

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

#endif


