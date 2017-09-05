#include "mainwindow.h"
#include <libgen.h>
#include <sys/stat.h>
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

void Calcul() //Coeur du programme
{
    double deltatemps, distmin, lpm;
    double thetarandom, phirandom;
    int aggcontact, newnumagg, finfichiersuivitempo, finmem = 0;
    //int tmp,superpo;
    int i, j, co;
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
            //NumAgg = Probabilite(contact, deltatemps);// Choice of an agrgegate, the probability of said agrgegate to be chosen proportionnal to his lpm
            double max = Aggregates.GetMaxTimeStep();
            if (contact)
                Aggregates.SortTimeSteps(max);
            NumAgg = Aggregates.RandomPick(deltatemps,Random());
            deltatemps = max/deltatemps;

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
            for (Aggregate* Agg : Aggregates)
            {
                Agg->Update();
            }
            lpm = Aggregates[NumAgg][4];
        }
        else
        {
            //$ Random Choice of an aggregate
            NumAgg = int(Random()*double(NAgg));
            deltatemps = 0.0;
            physicalmodel.temps = physicalmodel.temps + 1E-9;
            lpm = physicalmodel.Dpm*1E-9; //On fixe le lpm au Dpm
        }

        double distmove = lpm;

        //$ looking for potential contacts
        aggcontact = Aggregates.DistanceToNextContact(NumAgg, Vectdir, distmin);

        contact = (aggcontact > 0);
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
            //$ Aggregates in contact are reunited;
            newnumagg = Aggregates.Merge(NumAgg,aggcontact);
            NAgg--;

            Aggregates[newnumagg].Update();

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

    for (i = 0; i < NAgg; i++)
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
    Aggregates.Init(physicalmodel, NAgg);

    for (int i = 0; i < NAgg; i++)
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
        for (int k = 0; k <= i-1; k++)
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

    for (i=0;i<physicalmodel.N;i++)
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

    for (i = 0; i < NAgg; i++)
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


