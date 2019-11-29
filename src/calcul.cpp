#include "calcul.hpp"
#include <cmath>
#include <iomanip>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

namespace MCAC{

const double PI = atan(1.0)*4; // 3.14159

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
//    StatisticStorage Stats(physicalmodel);

    Init(physicalmodel,
//        Stats,
        Aggregates);

    // Contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool contact(true);

    size_t multiply_threshold = Aggregates.size() / 8;

    //$ Loop on the N monomeres
    while (! physicalmodel.Finished(Aggregates.size(), Aggregates.GetAvg_npp()))
    {
        if(contact)
        {
          Aggregates.spheres.save();
          Aggregates.save();
//            Stats.Analyze(Aggregates);
            //Stats.Save();

            if (Aggregates.size() <= multiply_threshold)
            {
                cout << "Duplication : "    << Aggregates.spheres.size()
                     << " spheres in "      << Aggregates.size() << " aggregates";

                Aggregates.Duplication();
                multiply_threshold = Aggregates.size() / 8;

                cout << " duplicated into "     << Aggregates.spheres.size()
                     << " spheres in "          << Aggregates.size() << " aggregates" << endl;
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
            /*
            //$ Choice of an aggregate according to his MFP
            double max = Aggregates.GetMaxTimeStep();
            if (contact)
            {
                Aggregates.SortTimeSteps(max);
            }
            NumAgg = Aggregates.RandomPick(Random());
            deltatemps = Aggregates.GetTimeStep(max);
            */
            NumAgg = Aggregates.PickLast();
            deltatemps = Aggregates[NumAgg].GetTimeStep();
        }
        else
        {
            //$ Random Choice of an aggregate
            NumAgg = size_t(Random()*double(Aggregates.size()));
            deltatemps = 1e-9;
        }
        double lpm = Aggregates[NumAgg].GetLpm();

        //$ looking for potential contacts
        double distmove;
        int aggcontact = Aggregates.DistanceToNextContact(NumAgg, Vectdir, distmove);

        // adjust time step
        deltatemps = deltatemps*distmove/lpm;

        //$ Translation of the aggregate
        for (size_t i = 0; i < 3; i++)
        {
            Vectdir[i] = Vectdir[i]*distmove;
        }
        Aggregates[NumAgg].Translate(Vectdir);
        Aggregates[NumAgg].TimeForward(deltatemps);

        //$ Time incrementation
        deltatemps = deltatemps / double(Aggregates.size());
        physicalmodel.Time = physicalmodel.Time + deltatemps;

        contact = (aggcontact >= 0);
        if (contact)
        {
            //$ Aggregates in contact are reunited;
            NumAgg = Aggregates.Merge(NumAgg, size_t(aggcontact));
        }

        //$ Update of the Aggregates/Spheres
        if (physicalmodel.ActiveModulephysique)
        {
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
            cout << "  Npp_avg=" << setw(4) << Aggregates.GetAvg_npp()
                 << "  NAgg=" << setw(4) << Aggregates.size()
                 << "  Time=" << setw(4) << physicalmodel.Time << " s"
                 << "   CPU=" << setw(4) << elapse << " s"
                 << " after " << setw(4) << physicalmodel.Wait << " it --- ";

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
                cout << "1.000e+00 * x^ 1.000e+00  --- r= 0" << endl;
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
//    Stats.Analyze(Aggregates);
    //Stats.Save(true);

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
//    Stats.print();
//    cout << endl;

    //Fermeture();

    cout << "\nThe End\n" << endl;
}


void Init(PhysicalModel& physicalmodel,
//    StatisticStorage& Stats,
    ListAggregat& Aggregates)
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

        //random size
        double x = Random();
        double Dp =0;
        if (physicalmodel.Mode == 1)
        {
            Dp = physicalmodel.Dpm+sqrt(2.0)*physicalmodel.sigmaDpm*inverf(2.0*x-1.0); //Loi normale
        }
        else
        {
            Dp = exp(log(physicalmodel.Dpm)+sqrt(2.0)*log(physicalmodel.sigmaDpm)*inverf(2.0*x-1.0)); //Loi log-normale
        }

        Dp = Dp/1E9;
        if (Dp <= 0)
        {
            Dp = physicalmodel.Dpm*1E-9;
        }

        bool placed = false;
        while (!placed)
        {
            //random position
            array<double, 3> newpos{{Random()*physicalmodel.L,
                                     Random()*physicalmodel.L,
                                     Random()*physicalmodel.L}};

            //++++++++++++ Test de superposition des sphérules lors de leur génération aléatoire ++++++++++++
            if (Aggregates.TestFreeSpace(newpos, Dp))
            {
                i--;
                testmem++;
            }
            else
            {
                Aggregates[i].Init(physicalmodel,Aggregates.verlet,newpos,i,Aggregates.spheres,Dp);
                placed = true;
            }
            if (testmem > Aggregates.size())
            {
                cout << "Impossible de générer tous les monomères sans superposition." << endl;
                cout << "La fraction volumique doit être diminuée." << endl;
                exit(0);
            }
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }
    }


//    Stats.Init();
}

}  // namespace MCAC


