#include "mainwindow.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libgen.h>
#include <list>
#include <sstream>
#include <string>
#include <sys/stat.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

namespace fs = std::experimental::filesystem;
using namespace std;

#ifdef WITH_GUI
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <qdir.h>
#include <QString>
#endif


namespace DLCA{

#ifdef WITH_GUI
MainWindow* GUI;
#endif

const double PI = atan(1.0)*4; // 3.14159

/*
// suivit tempo
string FichierSuiviTempo = "___";
double** tab; // generic array
int iValTab=0;
int nb_line;


char CheminSauve[500];
char commentaires[500];
*/

void InitRandom()
{
    time_t t;
    time(&t);
    srand(uint(t));
    //srand(0);
}

double Random()
{
    double v = rand();
    v = v/RAND_MAX;
    return v;
}

void Calcul(PhysicalModel& physicalmodel) //Coeur du programme
{
    ListAggregat Aggregates;
    StatisticStorage Stats(physicalmodel);

    Init(physicalmodel, Stats, Aggregates);

    // Contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool contact(true);

    //$ Loop on the N monomeres
    while (! physicalmodel.Finished(Aggregates.size()))
    {
#ifdef WITH_GUI
        qApp->processEvents(); //Permet de rafraichir la fenêtre Qt
#endif


        if(contact)
        {
            Aggregates.spheres.save();
            Aggregates.save();
            Stats.Analyze(Aggregates);
            //Stats.save();
        }


        // -- Generating a random direction --
        double thetarandom = Random()*PI*2;
        double phirandom = acos(1-2*Random());
        array<double,3> Vectdir{{sin(phirandom)*cos(thetarandom),
                                 sin(phirandom)*sin(thetarandom),
                                 cos(phirandom)}};

        // -- Pick an aggregate and it's corresponding timestep --
        double deltatemps(0);
        size_t NumAgg(0);
        if (physicalmodel.ActiveModulephysique)
        {
            //$ Choice of an aggregate according to his MFP
            double max = Aggregates.GetMaxTimeStep();
            if (contact)
            {
                Aggregates.SortTimeSteps(max);
            }
            NumAgg = Aggregates.RandomPick(deltatemps,Random());
            deltatemps = max/deltatemps;
        }
        else
        {
            //$ Random Choice of an aggregate
            NumAgg = size_t(Random()*double(Aggregates.size()));
            deltatemps = 1e-9;
        }

        //$ looking for potential contacts
        double distmove;
        int aggcontact = Aggregates.DistanceToNextContact(NumAgg, Vectdir, distmove);
        contact = (aggcontact >= 0);

        //$ Translation of the aggregate
        for (size_t i = 0; i < 3; i++)
        {
            Vectdir[i] = Vectdir[i]*distmove;
        }

        Aggregates[NumAgg].Translate(Vectdir);

        double lpm = Aggregates[NumAgg].GetLpm();
        if (contact)
        {
            //$ Aggregates in contact are reunited;
            NumAgg = Aggregates.Merge(NumAgg,size_t(aggcontact));
        }

        // adjust time step
        deltatemps = deltatemps*distmove/lpm;

        //$ Time incrementation
        physicalmodel.Time = physicalmodel.Time + deltatemps;

        //$ Update of the Aggregates/Spheres
        if (physicalmodel.ActiveModulephysique)
        {
/*
            if (physicalmodel.ActiveVariationTempo)
            {
                //$ Get Temperature in file
                finfichiersuivitempo = rechercheValTab();
                if (finmem != finfichiersuivitempo)
                {
                    sprintf(commentaires, "Attention le suivi temporel est plus long que celui du fichier lu.\n");
                    print(commentaires);
                }
                finmem = finfichiersuivitempo;
            }
*/
            if (physicalmodel.Asurfgrowth > 0.)
            {
	        //$ Growth of all spheres
                Aggregates.spheres.CroissanceSurface(deltatemps);
 
                //$ Aggregates update
                for (Aggregate* Agg : Aggregates)
                {
                    Agg->Update();
                }
            }

/*
            for (int i = 0; i<Aggregates.size(); i++)
            {
                for (int j = i+1; j<Aggregates.size(); j++)
                {
                    if(Aggregates[i].Contact(Aggregates[j]))
                    {
                        cout << "New contact !"<<endl;
                        Aggregates.Merge(i,j);
                    }
                }
            }
*/

        }

        //$ Show progress
        if(contact)
        {
            time_t t;
            time(&t);
            int secondes = int(round(t-physicalmodel.CPUStart));
            cout << "NAgg=" << Aggregates.size() << "    "
                 << "Time=" << physicalmodel.Time*1E6 << " E-6 s    "
                 << "CPU=" << secondes << "    "
                 << "after " << physicalmodel.Wait << " it " << endl;
            physicalmodel.Wait = 0;
#ifdef WITH_GUI
             if (!(GUI == nullptr)            )
             {
                GUI->progress(physicalmodel.N-Aggregates.size()+1);
             }
#endif
        }
        else
        {
            physicalmodel.Wait++;
        }
    }

    Aggregates.spheres.save(true);
    Aggregates.save(true);
    Stats.Analyze(Aggregates);
    //Stats.save(true);

    cout << "Final number of aggregates : " << Aggregates.size() << endl;

    /*
    cout << "Aggregats" << endl;
    for (size_t i = 0; i < Aggregates.size(); i++)
    {
        cout << i << "\t";
        const array<double, 3> pos = Aggregates[i].GetPosition();
        for (size_t j = 0; j < 3; j++)
        {
            cout << pos[j]*1E9 << "\t";
        }
        cout << Aggregates[i].GetVolAgregat()*1E25 << endl;
    }
    */
    cout << endl;
    Stats.print();
    cout << endl;

    //Fermeture();

    print("\nThe End\n");
}


void Init(PhysicalModel& physicalmodel, StatisticStorage& Stats, ListAggregat& Aggregates)
{

    //Initialize physical model


    physicalmodel.Init();
    if (physicalmodel.ActiveModulephysique)
    {
        print("Physical module activated.");
    }
    else
    {
        print("Physical module not activated.");
    }

    if (physicalmodel.ActiveVariationTempo)
    {
        //LectureSuiviTempo();
        //print("Le fichier de données de suivi temporel est lu.");
        print("Temporal variation of physical parameters not implemented.\n");
        exit(50);
    }
    else
    {
        print("No temporal variation of physical parameters.\n");
    }

    //Initialize randomness
    InitRandom();

    //Initialize the aggregates

    size_t testmem = 0;
    Aggregates.Init(physicalmodel, physicalmodel.N);


    for (size_t i = 0; i < Aggregates.size(); i++)
    {          

        //random position
        array<double, 3> newpos{{Random()*physicalmodel.L,
                                 Random()*physicalmodel.L,
                                 Random()*physicalmodel.L}};

        //random size
        double x = Random();
        double Dp =0;
        if (physicalmodel.Mode == 1)
        {
            Dp = physicalmodel.Dpm+sqrt(2.0)*physicalmodel.sigmaDpm*inverf(2*x-1); //Loi normale
        }
        else
        {
            Dp = exp(log(physicalmodel.Dpm)+sqrt(2.0)*log(physicalmodel.sigmaDpm)*inverf(2*x-1)); //Loi log-normale
        }

        Dp = Dp/1E9;
        if (Dp <= 0)
        {
            Dp = physicalmodel.Dpm*1E-9;
        }

        //++++++++++++ Test de superposition des sphérules lors de leur génération aléatoire ++++++++++++
        size_t test=0;
        if (i>0)
        {
            for (size_t k = 0; k <= i-1; k++)
            {
                double dist = Aggregates.spheres[k].Distance(newpos); // Calcule la distance centre à centre entre le monomère k et tous les autres
                if (dist <= Aggregates.spheres[k].Radius()+Dp/2)
                {
                    test++;
                }
            }
        }
        testmem = testmem + test; //Comptabilise le nombre d'échecs à positionner une sphère sans superposition

        if (test > 0)
        {
            i--;
        }
        else
        {
            Aggregates[i].Init(physicalmodel,Aggregates.verlet,newpos,i,Aggregates.spheres,Dp);
        }
        if (testmem > Aggregates.size())
        {
            cout << "Impossible de générer tous les monomères sans superposition." << endl;
            cout << "La fraction volumique doit être diminuée." << endl;
            exit(0);
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    }


    Stats.Init();
}


/*

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

        //for (i = 0; i < nb_line; i++)
        //{
        //    for (j = 0; j < 2; j++)
        //    {
        //        cout << tab[i][j] << "  ";
        //    }
        //    cout << endl;
        //}

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

    while (tab[iValTab][0] < physicalmodel.Time && iValTab < (nb_line-1))
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

        const array<double, 3> pos = Aggregates[i].GetPosition();
        for (j = 0; j < 3; j++)
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
    */

bool locale_with_dots()
{
    static bool tested = false;
    static bool with_dots;

    if (tested)
    {
        return with_dots;
    }
    double testfloat = 1.5;
    string teststr1 = "1.5";
    string teststr2 = "1,5";
    double test1=atof(teststr1.c_str());
    double test2=atof(teststr2.c_str());

    if (fabs(test1-testfloat)<1e-3)
    {
        with_dots = true;
    }
    else if (fabs(test2-testfloat)<1e-3)
    {
        with_dots = false;
    }
    else
    {
        cout << "What locale are you using ?" << endl;
        exit(1);
    }
    return with_dots;
}

double latof(const char* _char)
{
    string mystring = _char;
    if (!locale_with_dots())
    {
        size_t f = mystring.find('.');
        if (f>0)
        {
            mystring.replace(f, 1, ",");
        }
    }
    return atof(mystring.c_str());
}


PhysicalModel LectureParams(const string& FichierParam)
{
    PhysicalModel physicalmodel;

    FILE* f;
    char sauve[500];

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
        char commentaires[500];
        char com[500]; // Char array used in the ASCII Save

        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the N parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.N=size_t(atoi(commentaires));
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
            cout << "I need the ActiveVariationTempo parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.ActiveVariationTempo=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the AggMin parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.AggMin=size_t(atoi(commentaires));
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the WaitLimit parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.WaitLimit=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the CPULimit parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.CPULimit=atoi(commentaires);
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the DeltaSauve parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",commentaires,com);
            physicalmodel.DeltaSauve=size_t(atoi(commentaires));
        }
        if( fgets(com, 500, f) == nullptr)
        {
            cout << "I need the output_dir parameter" << endl;
            exit(1);
        }
        else
        {
            sscanf(com,"%s  %s",sauve,com);
        }
        fclose(f);
    }
    if (physicalmodel.Mode == 1)
    {
        physicalmodel.X = pow(double(physicalmodel.N)*PI/6.0/physicalmodel.FV*(1.0+3.0*physicalmodel.sigmaDpm*physicalmodel.sigmaDpm/physicalmodel.Dpm/physicalmodel.Dpm),1.0/3.0); //Loi normale
    }
    else
    {
        physicalmodel.X = pow(double(physicalmodel.N)*PI/6.0/physicalmodel.FV*exp(9.0/2.0*log(physicalmodel.sigmaDpm)*log(physicalmodel.sigmaDpm)),1.0/3.0); //Loi log-normale
    }

    fs::path pathParam = extractPath(FichierParam);
    physicalmodel.CheminSauve = pathParam / sauve;

    if( ! fs::exists(physicalmodel.CheminSauve))
    {
        if( ! fs::create_directory(physicalmodel.CheminSauve))
        {
            cout << "Error creating directory " << physicalmodel.CheminSauve << endl;
            exit(1);
        }
    }
    else
    {
        if( ! fs::is_directory(physicalmodel.CheminSauve))
        {
            cout << "Error not a directory " << physicalmodel.CheminSauve << endl;
            exit(1);
        }
    }

    return physicalmodel;
}

fs::path extractPath(const string& filename)
{
    fs::path path = filename;
    fs::path fullpath = fs::absolute(path);

    if(! fs::exists(fullpath))
    {
        cout << "File does not exist\n" << endl;
        exit(1);
    }

    fs::path parentpath = fullpath.parent_path();

    return parentpath; //Cette variable ne retient que le chemin du fichier param
}

int No_GUI(int argc, char *argv[]){
    if(argc <=1)
    {
        cout << "Missing argument : param file" << endl;
        return 1;
    }

    string FichierParam = argv[1];

    PhysicalModel physicalmodel = LectureParams(FichierParam);

    Calcul(physicalmodel);

    return 0;
}


void print(const string str)
{
#ifdef WITH_GUI
    if (GUI == nullptr)
        cout << str << endl;
    else
        GUI->print(str.c_str());
#else
    cout << str << endl;
#endif
}

bool dirExists(const char *path)
{
    struct stat info{};

    if(stat( path, &info ) != 0)
    {
        return false;
    }
    if(info.st_mode & S_IFDIR)
    {
        return true;
    }
    return false;
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
    string FichierParam = QFileDialog::getOpenFileName(this,"Sélectionner le fichier de données DLCA",
                                                        "C:/Users/dlca/Desktop/DLCA_sous_Qt/",
                                                        "Fichier de paramètres (*.txt)");
    //Affiche le chemin dans la zone de texte
    ui->AfficheurRep->append(FichierParam);
    LectureParams(FichierParam);
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
}  // namespace DLCA


