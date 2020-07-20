/*
 * MCAC
 * Copyright (C) 2020 CORIA
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "calcul.hpp"
#include "tools/tools.hpp"
#include "exceptions.hpp"
#include "tools/contact_info.hpp"
#include <iomanip>
#include <iostream>


namespace mcac {
void calcul(PhysicalModel &physicalmodel, AggregatList &aggregates) {//Coeur du programme

    // contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool event(true);
    size_t multiply_threshold = aggregates.size() / 8;

    //$ Loop on the N monomeres
    while (!physicalmodel.finished(aggregates.size(), aggregates.get_avg_npp())) {
        if (physicalmodel.n_iter_without_event % physicalmodel.write_between_event_every == 0) {
            aggregates.spheres.save();
            aggregates.save();
        }
        if (event) {
            if (aggregates.size() <= multiply_threshold) {
                std::cout << "Duplication : " << aggregates.spheres.size()
                          << " spheres in " << aggregates.size() << " aggregates";
                aggregates.duplication();
                multiply_threshold = aggregates.size() / 8;
                std::cout << " duplicated into " << aggregates.spheres.size()
                          << " spheres in " << aggregates.size() << " aggregates" << std::endl;
            }
        }

        // -- Generating a random direction --
        double thetarandom = random() * 2 * _pi;
        double phirandom = std::acos(1 - 2 * random());
        std::array<double, 3> vectdir{{std::sin(phirandom) * std::cos(thetarandom),
                                       std::sin(phirandom) * std::sin(thetarandom),
                                       std::cos(phirandom)}};

        // -- Pick an aggregate and it's corresponding timestep --
        double deltatemps(0);
        size_t num_agg(0);
        if (physicalmodel.pick_method == PickMethods::PICK_RANDOM) {
            //$ Choice of an aggregate according to his MFP
            double max = aggregates.get_max_time_step();
            if (event) {
                aggregates.sort_time_steps(max);
            }
            num_agg = aggregates.pick_random();
            deltatemps = aggregates.get_time_step(max);
        } else if (physicalmodel.pick_method == PickMethods::PICK_LAST) {
            num_agg = aggregates.pick_last();
            deltatemps = aggregates[num_agg].get_time_step();
        }
        double full_distance = aggregates[num_agg].get_lpm();
        double move_distance = full_distance;

        //$ looking for potential contacts
        AggregateContactInfo next_contact = aggregates.distance_to_next_contact(num_agg, vectdir, full_distance);

        bool contact = next_contact <= full_distance;
        if (contact){
            move_distance = next_contact.distance;
            deltatemps = deltatemps * move_distance / full_distance;
        }

        //$ Translation of the aggregate
        aggregates[num_agg].translate(vectdir * move_distance);
        aggregates[num_agg].time_forward(deltatemps);

        //$ Time incrementation
        if (physicalmodel.pick_method == PickMethods::PICK_LAST) {
            deltatemps = deltatemps / double(aggregates.size());
        }
        physicalmodel.time = physicalmodel.time + deltatemps;

        if (contact) {
            //$ Aggregates in contact are reunited;
            num_agg = aggregates.merge(next_contact);
        }

        //$ update of the Aggregates/Spheres
        bool split = false;
        if (std::abs(physicalmodel.a_surfgrowth) > 0.) {
            //$ Growth of all spheres
            aggregates.croissance_surface(deltatemps);
            split = aggregates.split();

            //$ Aggregates update
            for (const auto& agg : aggregates) {
                agg->update();
            }
        }
        event = split || contact;

        auto[success, fractal_dimension, fractal_prefactor, error] = aggregates.get_instantaneous_fractal_law();
        if (!success) {
            fractal_dimension = physicalmodel.fractal_dimension;
            fractal_prefactor = physicalmodel.fractal_prefactor;
        }
        //$ Show progress
        if (event) {
            aggregates.refresh();
            clock_t now = clock();
            double elapse = double(now - physicalmodel.cpu_start) / CLOCKS_PER_SEC;
            std::cout.precision(3);
            std::cout << std::scientific << std::boolalpha;
            std::cout << " Npp_avg=" << std::setw(4) << aggregates.get_avg_npp()
                      << " NAgg="    << std::setw(4) << aggregates.size()
                      << " Time="    << std::setw(4) << physicalmodel.time << "s"
                      << " CPU="     << std::setw(4) << elapse << "s"
                      << " contact=" << std::setw(4) << contact
                      << " split="   << std::setw(4) << split
                      << " after "   << std::setw(4) << physicalmodel.n_iter_without_event << " it"
                      << " --- "     << fractal_prefactor << " * x^ " << fractal_dimension << "  --- r= " << error
                      << std::endl;
            physicalmodel.n_iter_without_event = 0;
        } else {
            physicalmodel.n_iter_without_event++;
        }
        if (contact) {
            aggregates.remove_sphere(0);
            if(aggregates.split()) {
                for (Aggregate *agg : aggregates) {
                    agg->update();
                }
            }
            if (aggregates.spheres.size() >= 101) {
                try {
                    aggregates.add(1);
                }
                catch (TooDenseError const &e) {
                    std::cout << e.what() << std::endl;
                    break;
                }
            }
        }
    }
    aggregates.spheres.save(true);
    aggregates.save(true);
    std::cout << "Final number of aggregates : " << aggregates.size() << std::endl;
    std::cout << "Output files saved on: " << physicalmodel.output_dir << std::endl;
    std::cout << std::endl;
    std::cout << "\nThe End\n" << std::endl;
}
}  // namespace mcac


