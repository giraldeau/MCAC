#include "calcul.hpp"
#include "tools/tools.hpp"
#include <iomanip>
#include <iostream>


namespace mcac {
void calcul(PhysicalModel *physicalmodel, AggregatList *aggregates) {//Coeur du programme

    // contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool contact(true);
    size_t multiply_threshold = aggregates->size() / 8;

    //$ Loop on the N monomeres
    while (!physicalmodel->finished(aggregates->size(), aggregates->get_avg_npp())) {
        if (contact) {
            aggregates->refresh();
            aggregates->spheres.save();
            aggregates->save();
            if (aggregates->size() <= multiply_threshold) {
                std::cout << "Duplication : " << aggregates->spheres.size()
                          << " spheres in " << aggregates->size() << " aggregates";
                aggregates->duplication();
                multiply_threshold = aggregates->size() / 8;
                std::cout << " duplicated into " << aggregates->spheres.size()
                          << " spheres in " << aggregates->size() << " aggregates" << std::endl;
            }
        }

        // -- Generating a random direction --
        double thetarandom = random() * 2 * _pi;
        double phirandom = acos(1 - 2 * random());
        std::array<double, 3> vectdir{{sin(phirandom) * cos(thetarandom),
                                       sin(phirandom) * sin(thetarandom),
                                       cos(phirandom)}};

        // -- Pick an aggregate and it's corresponding timestep --
        double deltatemps(0);
        size_t num_agg(0);
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
        num_agg = aggregates->pick_last();
        deltatemps = (*aggregates)[num_agg].get_time_step();
        double lpm = (*aggregates)[num_agg].get_lpm();

        //$ looking for potential contacts
        auto contact_distance = aggregates->distance_to_next_contact(num_agg, vectdir);
        int aggcontact = contact_distance.first;
        double distmove = contact_distance.second;

        // adjust time step
        deltatemps = deltatemps * distmove / lpm;

        //$ Translation of the aggregate
        for (size_t i = 0; i < 3; i++) {
            vectdir[i] = vectdir[i] * distmove;
        }
        (*aggregates)[num_agg].translate(vectdir);
        (*aggregates)[num_agg].time_forward(deltatemps);

        //$ Time incrementation
        deltatemps = deltatemps / double(aggregates->size());
        physicalmodel->time = physicalmodel->time + deltatemps;
        contact = (aggcontact >= 0);
        if (contact) {
            //$ Aggregates in contact are reunited;
            num_agg = aggregates->merge(num_agg, size_t(aggcontact));
        }

        //$ update of the Aggregates/Spheres
        if (physicalmodel->a_surfgrowth > 0.) {
            //$ Growth of all spheres
            aggregates->spheres.croissance_surface(deltatemps);

            //$ Aggregates update
            for (Aggregate *agg : *aggregates) {
                agg->update();
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
            double elapse = double(now - physicalmodel->cpu_start) / CLOCKS_PER_SEC;
            std::cout.precision(3);
            std::cout << std::scientific;
            std::cout << "  Npp_avg=" << std::setw(4) << aggregates->get_avg_npp()
                      << "  NAgg=" << std::setw(4) << aggregates->size()
                      << "  Time=" << std::setw(4) << physicalmodel->time << " s"
                      << "   CPU=" << std::setw(4) << elapse << " s"
                      << " after " << std::setw(4) << physicalmodel->n_iter_without_contact << " it --- ";
            auto [success, dfe, kfe, error] = aggregates->get_instantaneous_fractal_law();
            if (success) {
                std::cout << "  " << kfe << " * x^ " << dfe << "  --- r= " << error << std::endl;
                /*
                physicalmodel.dfe = dfe;
                physicalmodel.kfe = kfe;
                */
            } else {
                std::cout << "1.000e+00 * x^ 1.000e+00  --- r= 0" << std::endl;
            }
            physicalmodel->n_iter_without_contact = 0;
        } else {
            physicalmodel->n_iter_without_contact++;
        }
    }
    aggregates->spheres.save(true);
    aggregates->save(true);
    std::cout << "Final number of aggregates : " << aggregates->size() << std::endl;
    std::cout << std::endl;
    std::cout << "\nThe End\n" << std::endl;
}
}  // namespace mcac


