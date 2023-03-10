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
#include "tools/contact_info.hpp"
#include "tools/tools.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace mcac {
/********************************************************************************
 * write time-resolved data to a text file
 ********************************************************************************/
static void save_advancement(PhysicalModel &physicalmodel, AggregatList &aggregates) noexcept {
    std::ofstream outfile;
    outfile.open(physicalmodel.output_dir / "advancement.dat", std::ios_base::app);
    outfile << physicalmodel.time;
    outfile << " " << physicalmodel.aggregate_concentration;
    outfile << " " << physicalmodel.volume_fraction;
    outfile << " " << aggregates.get_avg_npp();
    outfile << " " << physicalmodel.temperature;
    outfile << " " << physicalmodel.box_volume;
    outfile << " " << physicalmodel.monomer_concentration;
    outfile << " " << physicalmodel.u_sg;
    outfile << " " << physicalmodel.flux_nucleation;
    outfile << std::endl;
    outfile.close();
}
void print_bool(bool the_bool, int width, std::string txt) {
    if (the_bool) {
        std::cout << std::setw(width / 2 + width % 2) << "X" << std::setw(width / 2 + 3) << " | ";
        // std::cout << std::setw(width) << txt << " | ";
    } else {
        std::cout << std::setw(width + 3) << " | ";
    }
}
/********************************************************************************
 * Calcul: The main function of MCAC
 ********************************************************************************/
void calcul(PhysicalModel &physicalmodel, AggregatList &aggregates) {
    // contact is initialized to true for saving the initial set of monomeres and to sort the timesteps
    bool event(true);
    size_t duplication_threshold = aggregates.size() / 8;
    size_t reduction_threshold = aggregates.size() * 8;
    size_t total_events(0);
    size_t write_phys_time_int(0);

    physicalmodel.print();

    //$ Loop on the N monomeres
    while (!physicalmodel.finished(aggregates.size(), aggregates.get_avg_npp())) {
        if (physicalmodel.time_to_write(total_events)) {
            aggregates.spheres.save();
            aggregates.save();
            save_advancement(physicalmodel, aggregates);
        }
        if (event) {
            if (physicalmodel.with_domain_duplication && aggregates.size() <= duplication_threshold
                && !(physicalmodel.u_sg < 0.0)) {
                std::cout << "Duplication : " << aggregates.spheres.size() << " spheres in " << aggregates.size()
                          << " aggregates";
                aggregates.duplication();
                // duplication_threshold = aggregates.size() / 8;
                std::cout << " duplicated into " << aggregates.spheres.size() << " spheres in " << aggregates.size()
                          << " aggregates" << std::endl;
            }

            if (physicalmodel.with_domain_reduction && aggregates.size() >= reduction_threshold) {
                std::cout << "Reduction : " << aggregates.spheres.size() << " spheres in " << aggregates.size()
                          << " aggregates";
                aggregates.reduction();
                // reduction_threshold = aggregates.size() * 8;
                std::cout << " reduced to " << aggregates.spheres.size() << " spheres in " << aggregates.size()
                          << " aggregates" << std::endl;
            }
        }

        // -- Pick an aggregate and it's corresponding timestep --
        double deltatemps(0);
        size_t num_agg(0);
        if (physicalmodel.pick_method == PickMethods::PICK_RANDOM) {
            //$ Choice of an aggregate according to his MFP
            double max = aggregates.get_max_time_step();
            if (event || physicalmodel.with_surface_reactions || physicalmodel.with_flame_coupling) {
                aggregates.sort_time_steps(max);
            }
            num_agg = aggregates.pick_random();
            deltatemps = aggregates.get_time_step(max);
        } else if (physicalmodel.pick_method == PickMethods::PICK_LAST) {
            num_agg = aggregates.pick_last();
            deltatemps = aggregates[num_agg]->get_time_step();
        }
        double deltatemps_indiv = aggregates[num_agg]->get_time_step();
        double full_distance = aggregates[num_agg]->get_lpm();

        //$ L2 - Loop on orientations
        std::array<double, 3> vectdir = {{0, 0, 0}};
        bool effective_move = false;
        double move_distance = full_distance;
        int n_try(0);
        bool contact = false;
        AggregateContactInfo next_contact;

        while (!effective_move) {
            //$ -- Generating a random direction --
            vectdir = random_direction();
            effective_move = true;
            move_distance = full_distance;
            n_try++;
            if (physicalmodel.with_collisions) {
                //$ looking for potential geometric contacts
                next_contact = aggregates.distance_to_next_contact(num_agg, vectdir, full_distance);
                contact = next_contact <= full_distance;
                if (contact) {
                    move_distance = next_contact.distance;
                    if (physicalmodel.with_potentials) {
                        //$ Check repulsion/rebound (Coulomb or Pauli)
                        InterPotentialRegime regime_potential = aggregates.check_InterPotentialRegime(next_contact);
                        //$ Update contact
                        if (regime_potential != InterPotentialRegime::STICKING) {
                            effective_move = false;
                        }
                    }
                }
            }
        }

        //$ Translation of the aggregate
        aggregates[num_agg]->translate(vectdir * move_distance);

        //$ Time incrementation
        deltatemps = deltatemps * (move_distance / full_distance + static_cast<double>(n_try - 1));
        deltatemps_indiv = deltatemps_indiv * (move_distance / full_distance + static_cast<double>(n_try - 1));
        aggregates[num_agg]->time_forward(deltatemps);
        if (physicalmodel.pick_method == PickMethods::PICK_LAST) {
            deltatemps = deltatemps / double(aggregates.size());
        }
        physicalmodel.time = physicalmodel.time + deltatemps;

        bool split = false;
        bool disappear = false;
        //$ Surface reactions and splitting
        if (physicalmodel.with_surface_reactions) {
            if (physicalmodel.individual_surf_reactions) {
                disappear = aggregates.croissance_surface(deltatemps_indiv, num_agg);
                //$ Maybe will be modified later -> splitting only happens when u_sg is negative
                if ((physicalmodel.u_sg < 0.0) && !disappear) {
                    split = aggregates.split(num_agg);
                }
            } else {
                disappear = aggregates.croissance_surface(deltatemps);
                // Maybe will be modified later -> splitting only happens when u_sg is negative
                if (physicalmodel.u_sg < 0.0) {
                    split = aggregates.split();
                }
            }
        }

        //$ Merge
        bool merge = false;
        if (contact) {
            if (!physicalmodel.individual_surf_reactions || !(split || disappear)) {
                //$ Aggregates in contact are reunited;
                merge = aggregates.merge(next_contact);
            }
        }

        //$ Update
        if (physicalmodel.with_surface_reactions || physicalmodel.with_flame_coupling) {
            if (physicalmodel.n_iter_without_event % physicalmodel.full_aggregate_update_frequency == 0) {
                if (physicalmodel.individual_surf_reactions && !merge && !split && !disappear) {
                    //$ Update already done if merge or split (and not necessary if disappear)
                    aggregates[num_agg]->update();
                } else {
                    //$ Updating twice num_agg if merge or split (already done)
                    for (const auto &agg : aggregates) {
                        agg->update();
                    }
                }
            } else {
                if (physicalmodel.individual_surf_reactions && !merge && !split && !disappear) {
                    //$ Update already done if merge or split (and not necessary if disappear)
                    aggregates[num_agg]->update_partial();
                } else {
                    //$ Updating twice num_agg if merge or split (already done)
                    for (const auto &agg : aggregates) {
                        agg->update_partial();
                    }
                }
            }
        }

        //$ Nucleation
        bool nucleation(false);
        int monomers_to_add(0);
        if (physicalmodel.with_nucleation) {
            physicalmodel.nucleation(deltatemps);
            if (physicalmodel.nucleation_accum > 1.0) {
                nucleation = true;
                monomers_to_add = static_cast<int>(std::floor(physicalmodel.nucleation_accum));
                physicalmodel.nucleation_accum -= static_cast<double>(monomers_to_add);
                aggregates.add(monomers_to_add);
                aggregates.sort_time_steps(aggregates.get_max_time_step());
            }
        }

        event = split || merge || disappear || nucleation;
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

        //$ Show progress
        if (event) {
            physicalmodel.cpu_last_event = clock();
            if (total_events % 20 == 1) {
                std::cout << std::setw(8) << "#"
                          << " | " << std::setw(9) << "Npp_avg"
                          << " | " << std::setw(8) << "NAgg"
                          << " | " << std::setw(10) << "Time"
                          << " | " << std::setw(10) << "CPU"
                          << " | " << std::setw(7) << "contact"
                          << " | " << std::setw(5) << "merge"
                          << " | " << std::setw(5) << "split"
                          << " | " << std::setw(9) << "disappear"
                          << " | " << std::setw(10) << "nucleation"
                          << " | " << std::setw(8) << "after" << std::endl;
            }
            clock_t now = clock();
            double elapse = double(now - physicalmodel.cpu_start) / CLOCKS_PER_SEC;
            std::cout.precision(3);
            std::cout << std::scientific;
            std::cout << std::setw(8) << total_events << " | ";
            std::cout << std::setw(8) << aggregates.get_avg_npp() << " | ";
            std::cout << std::setw(8) << aggregates.size() << " | ";
            std::cout << std::setw(8) << physicalmodel.time << "s"
                      << " | ";
            std::cout << std::setw(8) << elapse << "s"
                      << " | ";
            print_bool(contact, 7, "contact");
            print_bool(merge, 5, "merge");
            print_bool(split, 5, "split");
            print_bool(disappear, 9, "disappear");

            std::cout << std::setw(10) << monomers_to_add << " | " << std::setw(8) << current_n_iter_without_event
                      << std::endl;
        }
        //$ Update physical model
        if (event || physicalmodel.with_surface_reactions) {
            physicalmodel.update(aggregates.size(),
                                 aggregates.spheres.size(),
                                 aggregates.get_total_volume(),
                                 aggregates.get_total_surface());
        }
        if (physicalmodel.with_flame_coupling) {
            physicalmodel.update_from_flame();
        }
    }
    save_advancement(physicalmodel, aggregates);
    aggregates.spheres.save();
    aggregates.save();
    std::cout << " Final residence time=" << std::setw(4) << physicalmodel.time << "s" << std::endl;
    std::cout << "Final number of aggregates : " << aggregates.size() << std::endl;
    std::cout << "Output files saved on: " << physicalmodel.output_dir << std::endl;
    std::cout << std::endl;
    std::cout << "\nThe End\n" << std::endl;
}
} // namespace mcac
