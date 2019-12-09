#include "calcul.hpp"
#include "tools/tools.hpp"
#include <iomanip>
#include <iostream>


using namespace std;
namespace MCAC {
void Calcul(PhysicalModel &physicalmodel, AggregatList &Aggregates) {//Coeur du programme

    // contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool contact(true);
    size_t multiply_threshold = Aggregates.size() / 8;

    //$ Loop on the N monomeres
    while (!physicalmodel.finished(Aggregates.size(), Aggregates.get_avg_npp())) {
        if (contact) {
            Aggregates.refresh();
            Aggregates.spheres.save();
            Aggregates.save();
            if (Aggregates.size() <= multiply_threshold) {
                cout << "Duplication : " << Aggregates.spheres.size()
                     << " spheres in " << Aggregates.size() << " aggregates";
                Aggregates.duplication();
                multiply_threshold = Aggregates.size() / 8;
                cout << " duplicated into " << Aggregates.spheres.size()
                     << " spheres in " << Aggregates.size() << " aggregates" << endl;
            }
        }

        // -- Generating a random direction --
        double thetarandom = random() * 2 * _pi;
        double phirandom = acos(1 - 2 * random());
        array<double, 3> Vectdir{{sin(phirandom) * cos(thetarandom),
                                  sin(phirandom) * sin(thetarandom),
                                  cos(phirandom)}};

        // -- Pick an aggregate and it's corresponding timestep --
        double deltatemps(0);
        size_t NumAgg(0);
        /*
        //$ Choice of an aggregate according to his MFP
        double max = Aggregates.get_max_time_step();
        if (contact)
        {
            Aggregates.sort_time_steps(max);
        }
        NumAgg = Aggregates.pick_random();
        deltatemps = Aggregates.get_time_step(max);
        */
        NumAgg = Aggregates.pick_last();
        deltatemps = Aggregates[NumAgg].get_time_step();
        double lpm = Aggregates[NumAgg].get_lpm();

        //$ looking for potential contacts
        auto contact_distance = Aggregates.distance_to_next_contact(NumAgg, Vectdir);
        int aggcontact = contact_distance.first;
        double distmove = contact_distance.second;

        // adjust time step
        deltatemps = deltatemps * distmove / lpm;

        //$ Translation of the aggregate
        for (size_t i = 0; i < 3; i++) {
            Vectdir[i] = Vectdir[i] * distmove;
        }
        Aggregates[NumAgg].translate(Vectdir);
        Aggregates[NumAgg].time_forward(deltatemps);

        //$ Time incrementation
        deltatemps = deltatemps / double(Aggregates.size());
        physicalmodel.time = physicalmodel.time + deltatemps;
        contact = (aggcontact >= 0);
        if (contact) {
            //$ Aggregates in contact are reunited;
            NumAgg = Aggregates.merge(NumAgg, size_t(aggcontact));
        }

        //$ update of the Aggregates/Spheres
        if (physicalmodel.a_surfgrowth > 0.) {
            //$ Growth of all spheres
            Aggregates.spheres.croissance_surface(deltatemps);

            //$ Aggregates update
            for (Aggregate *Agg : Aggregates) {
                Agg->update();
            }
        }
        //            for (int i = 0; i < Aggregates.size(); i++) {
        //                for (int j = i + 1; j < Aggregates.size(); j++) {
        //                    if (Aggregates[i].contact(Aggregates[j])) {
        //                        cout << "New contact !" << endl;
        //                        Aggregates.merge(i, j);
        //                    }
        //                }
        //            }

        //$ Show progress
        if (contact) {
            clock_t now = clock();
            double elapse = double(now - physicalmodel.cpu_start) / CLOCKS_PER_SEC;
            cout.precision(3);
            cout << scientific;
            cout << "  Npp_avg=" << setw(4) << Aggregates.get_avg_npp()
                 << "  NAgg=" << setw(4) << Aggregates.size()
                 << "  Time=" << setw(4) << physicalmodel.time << " s"
                 << "   CPU=" << setw(4) << elapse << " s"
                 << " after " << setw(4) << physicalmodel.n_iter_without_contact << " it --- ";
            auto InstantaneousFractalLaw = Aggregates.get_instantaneous_fractal_law();
            if (get<0>(InstantaneousFractalLaw)) {
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
            } else {
                cout << "1.000e+00 * x^ 1.000e+00  --- r= 0" << endl;
            }
            physicalmodel.n_iter_without_contact = 0;
        } else {
            physicalmodel.n_iter_without_contact++;
        }
    }
    Aggregates.spheres.save(true);
    Aggregates.save(true);
    cout << "Final number of aggregates : " << Aggregates.size() << endl;
    cout << endl;
    cout << "\nThe End\n" << endl;
}
}  // namespace MCAC


