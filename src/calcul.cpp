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
#include "tools/contact_info.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>


namespace mcac {
/********************************************************************************
* write time-resolved data to a text file
********************************************************************************/
static void save_advancement(PhysicalModel &physicalmodel, AggregatList &aggregates) noexcept {
    std::ofstream outfile;
    outfile.open(physicalmodel.output_dir / "advancement.dat", std::ios_base::app);
    outfile << physicalmodel.time
            << " " << physicalmodel.aggregate_concentration
            << " " << physicalmodel.volume_fraction
            << " " << aggregates.get_avg_npp()
            << " " << physicalmodel.temperature
            << std::endl;
    outfile.close();
}
/********************************************************************************
* Calcul: The main function of MCAC
********************************************************************************/
void calcul(PhysicalModel &physicalmodel, AggregatList &aggregates) {
    // contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool event(true);
    size_t multiply_threshold = aggregates.size() / 8;
    size_t total_events(0);

    // load flame-coupling info.
    FlameCoupling flame;
    if (physicalmodel.with_flame_coupling) {
        flame = FlameCoupling(physicalmodel.flame_file);
        physicalmodel.update_from_flame(flame);
    }
    physicalmodel.print();

    //$ Loop on the N monomeres
    while (!physicalmodel.finished(aggregates.size(), aggregates.get_avg_npp())) {
        if (physicalmodel.n_iter_without_event % physicalmodel.write_between_event_frequency == 0) {
            aggregates.spheres.save();
            aggregates.save();
            save_advancement(physicalmodel, aggregates);
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
        std::array<double, 3> vectdir = random_direction();

        // -- Pick an aggregate and it's corresponding timestep --
        double deltatemps(0);
        size_t num_agg(0);
        if (physicalmodel.pick_method == PickMethods::PICK_RANDOM) {
            //$ Choice of an aggregate according to his MFP
            double max = aggregates.get_max_time_step();
            if (event || physicalmodel.with_surface_reactions) {
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

        bool contact = false;
        AggregateContactInfo next_contact;
        if (physicalmodel.with_collisions) {
            //$ looking for potential contacts
            next_contact = aggregates.distance_to_next_contact(num_agg, vectdir, full_distance);
            contact = next_contact <= full_distance;
            if (contact) {
                move_distance = next_contact.distance;
                deltatemps = deltatemps * move_distance / full_distance;
            }
        }

        //$ Translation of the aggregate
        aggregates[num_agg].translate(vectdir * move_distance);
        aggregates[num_agg].time_forward(deltatemps);

        //$ Time incrementation
        if (physicalmodel.pick_method == PickMethods::PICK_LAST) {
            deltatemps = deltatemps / double(aggregates.size());
        }
        physicalmodel.time = physicalmodel.time + deltatemps;

        //$ Event management
        if (contact) {
            //$ Aggregates in contact are reunited;
            num_agg = aggregates.merge(next_contact);
        }
        bool split = false;
        if (physicalmodel.with_surface_reactions) {
            aggregates.croissance_surface(deltatemps);
            // Maybe will be modified later -> spliting only happens when u_sg is negative
            if(physicalmodel.u_sg < 0.0) {
                split = aggregates.split();
            }
        }
        event = split || contact;
        size_t current_n_iter_without_event = physicalmodel.n_iter_without_event;
        if (event) {
            physicalmodel.n_iter_without_event = 0;
            total_events++;
        } else {
            physicalmodel.n_iter_without_event++;
        }

        //$ Update aggregates
        if (event) {
            aggregates.refresh();
        }
        if (physicalmodel.with_surface_reactions 
            || physicalmodel.with_flame_coupling) {
            if (physicalmodel.n_iter_without_event % physicalmodel.full_aggregate_update_frequency == 0) {
                for (const auto &agg : aggregates) {
                    agg->update();
                }
            } else {
                for (const auto &agg : aggregates) {
                    agg->update_partial();
                }
            }
        }

        /* no calculation of fractal dimension in time
        auto[success, fractal_dimension, fractal_prefactor, error] = aggregates.get_instantaneous_fractal_law();
        if (!success) {
            fractal_dimension = physicalmodel.fractal_dimension;
            fractal_prefactor = physicalmodel.fractal_prefactor;
        }
        */

        //$ Show progress
        if (event) {
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
                      << " after "   << std::setw(4) << current_n_iter_without_event << " it"
                      << " total_events " << std::setw(4) <<  total_events
                      // << " --- "     << fractal_prefactor << " * x^ " << fractal_dimension << "  --- r= " << error
                      << std::endl;
        }
        //$ Update physical model
        if (event
            || physicalmodel.with_surface_reactions) {
            physicalmodel.update(aggregates.size(), aggregates.get_total_volume());
        }
        if (physicalmodel.with_flame_coupling) {
            physicalmodel.update_from_flame(flame);
        }
    }
    save_advancement(physicalmodel, aggregates);
    aggregates.spheres.save();
    aggregates.save();
    std::cout << " Final residence time="    << std::setw(4) << physicalmodel.time << "s" << std::endl;
    std::cout << "Final number of aggregates : " << aggregates.size() << std::endl;
    std::cout << "Output files saved on: " << physicalmodel.output_dir << std::endl;
    std::cout << std::endl;
    std::cout << "\nThe End\n" << std::endl;
}
}  // namespace mcac
