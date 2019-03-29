#include "calcul.hpp"
#include <cmath>
#include <iomanip>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

namespace DLCA{

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
    //srand(uint(t));
    srand(0);
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

    size_t multiply_threshold = Aggregates.size() / 8;

    //$ Loop on the N monomeres
    while (! physicalmodel.Finished(Aggregates.size()))
    {
        if(contact)
        {
            Aggregates.spheres.save();
            Aggregates.save();
            Stats.Analyze(Aggregates);
            //Stats.save();

            if (Aggregates.size() <= multiply_threshold)
            {
                cout << "Duplication : " << Aggregates.spheres.size()
                     << " spheres in " << Aggregates.size() << " aggregates";

                Aggregates.Multiply();
                multiply_threshold = Aggregates.size() / 8;

                cout << " duplicated into " << Aggregates.spheres.size()
                     << " spheres in " << Aggregates.size() << " aggregates" << endl;
            }

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
                    cout << commentaires << endl;
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
            clock_t now = clock();
            double elapse = double(now - physicalmodel.CPUStart) / CLOCKS_PER_SEC;
            cout.precision(3);
            cout << scientific;
            cout << "  NAgg=" << setw(5) << Aggregates.size()
                 << "  Time=" << setw(5) << physicalmodel.Time << " s"
                 << "   CPU=" << setw(5) << elapse << " s"
                 << " after " << setw(5) << physicalmodel.Wait << " it --- ";

            auto InstantaneousFractalLaw = Aggregates.getInstantaneousFractalLaw();
            if(get<0>(InstantaneousFractalLaw))
            {
                cout << "  "
                     << exp(get<2>(InstantaneousFractalLaw))
                     << " * x^ "
                     << get<1>(InstantaneousFractalLaw)
                     << "  --- r= "
                     << get<3>(InstantaneousFractalLaw) << endl;
                /*
                physicalmodel.dfe = get<1>(InstantaneousFractalLaw);
                physicalmodel.kfe = get<2>(InstantaneousFractalLaw);
                */
            }
            else {
                 cout << endl;
            }

            physicalmodel.Wait = 0;
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

    cout << "\nThe End\n" << endl;
}


void Init(PhysicalModel& physicalmodel, StatisticStorage& Stats, ListAggregat& Aggregates)
{

    //Initialize physical model


    physicalmodel.Init();
    if (physicalmodel.ActiveModulephysique)
    {
        cout << "Physical module activated." << endl;
    }
    else
    {
        cout << "Physical module not activated." << endl;
    }

    if (physicalmodel.ActiveVariationTempo)
    {
        //LectureSuiviTempo();
        //cout << "Le fichier de données de suivi temporel est lu." << endl;
        cout << "Temporal variation of physical parameters not implemented.\n" << endl;
        exit(50);
    }
    else
    {
        cout << "No temporal variation of physical parameters.\n" << endl;
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

    */

}  // namespace DLCA


