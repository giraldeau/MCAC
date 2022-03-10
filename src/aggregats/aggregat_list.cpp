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
#include "aggregats/aggregat_distance.hpp"
#include "aggregats/aggregat_list.hpp"
#include "spheres/sphere_contact.hpp"
#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include "exceptions.hpp"
#include <iostream>

namespace mcac
{
    [[gnu::pure]] double AggregatList::get_total_volume() const
    {
        double total_volume(0.0);
        for (const auto &agg : list)
        {
            total_volume += (*agg->agregat_volume);
        }
        return total_volume;
    }
    [[gnu::pure]] double AggregatList::get_total_surface() const
    {
        double total_surface(0.0);
        for (const auto &agg : list)
        {
            total_surface += (*agg->agregat_surface);
        }
        return total_surface;
    }
    [[gnu::pure]] double AggregatList::get_avg_npp() const
    {
        return avg_npp;
    }
    [[gnu::pure]] double AggregatList::get_max_time_step() const
    {
        return max_time_step;
    }
    [[gnu::pure]] double AggregatList::get_time_step(double max) const
    {
        double deltatemps = max / cumulative_time_steps[size() - 1];
        return deltatemps;
    }
    [[gnu::pure]] size_t AggregatList::pick_random() const
    {
        //$ Pick a random sphere
        double val_alea = random() * cumulative_time_steps[size() - 1];
        long n = std::lower_bound(cumulative_time_steps.begin(), cumulative_time_steps.end(), val_alea) - cumulative_time_steps.begin();
        size_t num_agg = index_sorted_time_steps[static_cast<size_t>(n)];
        return num_agg;
    }
    [[gnu::pure]] size_t AggregatList::pick_last() const
    {
        // TODO(pouxa): cache result
        double time = *list[0]->proper_time;
        size_t latest = 0;
        for (const auto &agg : list)
        {
            if (*agg->proper_time < time)
            {
                time = *agg->proper_time;
                latest = agg->label;
            }
        }
        return latest;
    }
    void AggregatList::add(size_t n)
    {
        size_t initial_n_agg = size();
        size_t initial_n_sph = spheres.size();
        spheres.add(n, spheres);
        spheres.setpointers();
        add(n, *this);
        setpointers();

        // initialize data
        for (size_t i = 0; i < n; i++)
        {
            size_t i_agg = initial_n_agg + i;
            size_t i_sph = initial_n_sph + i;
            list[i_agg]->init(i_agg, i_sph, true);
        }
        refresh();
    }
    void AggregatList::refresh()
    {
        max_time_step = *list[0]->time_step;
        for (const auto &agg : list)
        {
            max_time_step = std::max(*agg->time_step, max_time_step);
        }
        avg_npp = static_cast<double>(spheres.size()) / static_cast<double>(size());
    }
    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v)
    {
        // initialize original index locations
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        std::sort(idx.begin(), idx.end(),
                  [&v](size_t i_1, size_t i_2)
                  {
                      return v[i_1] < v[i_2];
                  });
        return idx;
    }
    void AggregatList::sort_time_steps(double factor)
    {
        std::vector<double> tp_t(size());
#pragma omp simd
        for (size_t i = 0; i < size(); i++)
        {
            tp_t[i] = factor / (*list[i]->time_step);
        }
        index_sorted_time_steps = sort_indexes(tp_t);
        cumulative_time_steps.resize(size());

        //$ Accumulate the timesteps
        cumulative_time_steps[0] = tp_t[index_sorted_time_steps[0]];
        for (size_t i = 1; i < size(); i++)
        {
            cumulative_time_steps[i] = cumulative_time_steps[i - 1] + tp_t[index_sorted_time_steps[i]];
        }
    }
    void AggregatList::duplication()
    {
        size_t old_n_agg = size();
        double old_l = physicalmodel->box_length;
        physicalmodel->box_length *= 2;
        physicalmodel->n_monomeres *= 8;
        physicalmodel->box_volume = std::pow(physicalmodel->box_length, 3);
        for (const auto &agg : list)
        {
            agg->unset_verlet();
        }
        // TODO(pouxa): Rework this in order not to to it aggregate by aggregate

        for (size_t iagg = 0; iagg < old_n_agg; iagg++)
        {
            for (int i = 0; i <= 1; i++)
            {
                for (int j = 0; j <= 1; j++)
                {
                    for (int k = 0; k <= 1; k++)
                    {
                        if (i != 0 or j != 0 or k != 0)
                        {
                            auto new_agg = add(*list[iagg], *this);
                            new_agg->label = size() - 1;
                            setpointers();
                            spheres.setpointers();
                            std::array<double, 3> vec_move = {i * old_l,
                                                              j * old_l,
                                                              k * old_l};
                            new_agg->translate(vec_move);
                        }
                    }
                }
            }
        }
        //$ update Verlet
        verlet = Verlet(physicalmodel->n_verlet_divisions, physicalmodel->box_length);
        for (const auto &agg : list)
        {
            agg->set_verlet(&verlet);
        }

        //$ Update Physical model
        physicalmodel->update(size(),
                              spheres.size(),
                              get_total_volume(),
                              get_total_surface());
    }
    void AggregatList::reduction()
    {

        //$ update box
        physicalmodel->box_length /= 2;
        physicalmodel->box_volume = std::pow(physicalmodel->box_length, 3);

        //$ remove aggregates outside the box
        auto aggregate = list.begin();
        while (aggregate != list.end())
        {
            auto index = std::distance(list.begin(), aggregate);

            auto [x, y, z] = (*aggregate)->get_position();

            if (x > physicalmodel->box_length ||
                y > physicalmodel->box_length ||
                z > physicalmodel->box_length)
            {

                auto nspheres = (*aggregate)->size();
                while (nspheres > 0)
                {
                    auto sphere = (*aggregate)->myspheres.begin() + nspheres - 1;
                    remove_sphere((*sphere)->get_index());
                    nspheres--;
                }
                spheres.setpointers();
                setpointers();

                aggregate = list.begin() + index;
            }
            else
            {
                aggregate = list.begin() + index + 1;
            }
        }

        //$ update Verlet
        for (const auto &agg : list)
        {
            agg->unset_verlet();
        }
        verlet = Verlet(physicalmodel->n_verlet_divisions, physicalmodel->box_length);
        for (const auto &agg : list)
        {
            agg->set_verlet(&verlet);
        }

        //$ remove aggregates that collide with others due to new periodicity
        aggregate = list.begin();
        while (aggregate != list.end())
        {
            auto index = std::distance(list.begin(), aggregate);

            auto [x, y, z] = (*aggregate)->get_position();

            double dist_x = physicalmodel->box_length - x;
            double dist_y = physicalmodel->box_length - y;
            double dist_z = physicalmodel->box_length - z;

            double rmax = (*aggregate)->get_rmax();
            bool new_contact = false;
            // only if the aggregate is straddling the new periodicity limits
            if (dist_x < rmax || dist_y < rmax || dist_z < rmax)
            {

                // Use Verlet to reduce search of potential collision
                std::vector<size_t> neighborhood(get_neighborhood(index, {0., 0., 0.}, 0.));

                // Assimilate Aggregate as sphere to drasticly speed-up search
                std::multimap<double, size_t> filtered_neighborhood(filter_neighborhood(index, {0., 0., 0.}, neighborhood, 0.));

                //$ loop on the aggregates potentially in contact
                for (auto suspect : filtered_neighborhood)
                {
                    auto [suspect_distance, id] = suspect;

                    if (contact(**aggregate, *list[id]))
                    {
                        auto [other_x, other_y, other_z] = list[id]->get_position();
                        bool over_x = (other_x - physicalmodel->box_length / 2) * (x - physicalmodel->box_length / 2) < 0;
                        bool over_y = (other_y - physicalmodel->box_length / 2) * (y - physicalmodel->box_length / 2) < 0;
                        bool over_z = (other_z - physicalmodel->box_length / 2) * (z - physicalmodel->box_length / 2) < 0;

                        new_contact = over_x || over_y || over_z;
                    }
                    if (new_contact)
                    {
                        break;
                    }
                }
            }

            if (new_contact)
            {
                std::cout << " WARNING: removing aggregate during reduction due to new overlap" << std::endl;

                auto nspheres = (*aggregate)->size();
                while (nspheres > 0)
                {
                    auto sphere = (*aggregate)->myspheres.begin() + nspheres - 1;
                    remove_sphere((*sphere)->get_index());
                    nspheres--;
                }
                spheres.setpointers();
                setpointers();

                aggregate = list.begin() + index;
            }
            else
            {
                aggregate = list.begin() + index + 1;
            }
        }

        //$ Update Physical model
        physicalmodel->update(size(),
                              spheres.size(),
                              get_total_volume(),
                              get_total_surface());
    }
    InterPotentialRegime AggregatList::check_InterPotentialRegime(AggregateContactInfo contact_info)
    {
        auto moving_sphere = contact_info.moving_sphere.lock();
        auto other_sphere = contact_info.other_sphere.lock();

        //$ 1. Colliding primary spheres properties
        double D_moving = 2.0 * moving_sphere->get_radius();
        double D_other = 2.0 * other_sphere->get_radius();

        //$ 2. Electric charges
        /*    if (physicalmodel->with_dynamic_random_charges) {
                list[contact_info.moving_aggregate]->set_random_charge(kbT, intpotential_info);
                list[contact_info.other_aggregate]->set_random_charge(kbT, intpotential_info);
            }
        */

        //$ 3. Energy barrier/well depth
        double P_stick(1.0), P_coll(1.0);
        if (physicalmodel->with_external_potentials)
        {
            auto moving_aggregate = contact_info.moving_aggregate.lock();
            auto other_aggregate = contact_info.other_aggregate.lock();
            int q_moving = moving_aggregate->get_electric_charge();
            int q_other = other_aggregate->get_electric_charge();
            auto [E_bar, E_well] = physicalmodel->intpotential_info.get_Ebar_Ewell(D_moving, D_other, q_moving, q_other);
            double E_stick = std::abs(E_well) + std::abs(E_bar);
            P_stick = std::erf(std::sqrt(E_stick)) - std::sqrt(E_stick) * std::exp(-E_stick);
            P_coll = 1.0 - std::erf(std::sqrt(E_bar)) + std::sqrt(E_bar) * std::exp(-E_bar);
        }
        else
        {
            // Hou et al., J. Aerosol Sci. (2020) 105478.
            double kbT = _boltzmann * (physicalmodel->temperature);
            double D = D_moving * D_other / (D_moving + D_other);
            D = D * (1e+09);
            double E_well = (-6.6891e-23) * std::pow(D, 3) +
                            (1.1244e-21) * std::pow(D, 2) +
                            (1.1394e-20) * D - 5.5373e-21;
            P_stick = 1.0 - (1.0 + std::abs(E_well) / kbT) * std::exp(-std::abs(E_well) / kbT);
        }

        //$ 4. Check - Coulomb repulsion
        if (random() > P_coll)
        {
            return InterPotentialRegime::REPULSION;
        }

        //$ 5. Check - Pauli repulsion
        if (random() > P_stick)
        {
            return InterPotentialRegime::BOUNCING;
        }
        return InterPotentialRegime::STICKING;
    }
    bool AggregatList::merge(AggregateContactInfo contact_info)
    {
        auto moving_sphere = contact_info.moving_sphere.lock();
        auto other_sphere = contact_info.other_sphere.lock();
        if (!moving_sphere || !other_sphere)
        {
            return false;
        }
        if (!contact(*moving_sphere, *other_sphere))
        {
            return false;
        }

        auto _keeped = static_cast<const size_t>(std::min(moving_sphere->agg_label,
                                                          other_sphere->agg_label));
        auto _removed = static_cast<const size_t>(std::max(moving_sphere->agg_label,
                                                           other_sphere->agg_label));

        // compute proper time of the final aggregate
        // keeping global time constant
        double newtime = (*list[_keeped]->proper_time) + (*list[_removed]->proper_time) - physicalmodel->time;
        int total_charge = list[_keeped]->get_electric_charge() +
                           list[_removed]->get_electric_charge();

        // merge the two aggregate but do not remove the deleted one
        if (!list[_keeped]->merge(list[_removed], contact_info))
        {
            throw MergeError("ListAggregate want to merge but the aggregate refuses");
        }
        remove(_removed);
        setpointers();
        *list[_keeped]->proper_time = newtime;

        if (physicalmodel->with_dynamic_random_charges)
        {
            list[_keeped]->electric_charge = physicalmodel->get_random_charge(*(list[_keeped]->d_m));
        }
        else
        { // electric charges preservation
            list[_keeped]->electric_charge = total_charge;
        }

        return true;
    }
    bool AggregatList::split()
    {
        bool has_split = false;
        auto suspect = list.begin();
        while (suspect != list.end())
        {
            auto index = std::distance(list.begin(), suspect);
            // the split function will create the new aggregates
            if ((*suspect)->split())
            {
                has_split = true;
                // but we still have to destroy the current one
                remove(size_t(index));
                setpointers();
                suspect = list.begin() + index;
            }
            else
            {
                suspect = list.begin() + index + 1;
            }
        }
        return has_split;
    }
    bool AggregatList::split(const size_t numagg)
    {
        // The split function will create the new aggregates
        if (list[numagg]->split())
        {
            // but we still have to destroy the current one
            remove(numagg);
            setpointers();
            return true;
        }
        return false;
    }
    //################################# Determination of the contacts between aggregates ####################################
    AggregateContactInfo AggregatList::distance_to_next_contact(const size_t source,
                                                                const std::array<double, 3> &direction,
                                                                const double distance) const
    {
        // Use Verlet to reduce search
        std::vector<size_t> neighborhood(get_neighborhood(source, direction, distance));

        // Assimilate Aggregate as sphere to drasticly speed-up search
        std::multimap<double, size_t> filtered_neighborhood(filter_neighborhood(source, direction, neighborhood, distance));
        AggregateContactInfo closest_contact; // infinity by default

        //$ loop on the aggregates potentially in contact
        for (auto suspect : filtered_neighborhood)
        { // For every aggregate that could be in contact
            auto [suspect_distance, id] = suspect;

            if (closest_contact.distance <= 0.)
            {
                // cannot be closest than that
                break;
            }

            if (closest_contact.distance < suspect_distance)
            {
                // We already found the closest one
                break;
            }
            AggregateContactInfo potential_contact = distance_to_contact(list[source],
                                                                         list[id],
                                                                         direction,
                                                                         distance);
            if (potential_contact < closest_contact)
            {
                closest_contact = potential_contact;
            }
        }
        return closest_contact;
    }
    std::vector<size_t> AggregatList::get_neighborhood(const size_t source,
                                                       const std::array<double, 3> &direction,
                                                       const double distance) const
    {
        std::vector<size_t> neighborhood;

        // Extract from verlet
        double mindist(*list[source]->rmax + maxradius);
        std::array<double, 3> sourceposition = list[source]->get_position();
        std::array<double, 3> vector{distance * direction};
        neighborhood = verlet.get_neighborhood(sourceposition, vector, mindist);

        // Remove me
        for (size_t i = 0; i < neighborhood.size(); i++)
        {
            if (neighborhood[i] == source)
            {
                neighborhood.erase(neighborhood.begin() + long(i));
                return neighborhood;
            }
        }
        throw VerletError("Aggregate not on the verlet list ???");
        // return neighborhood;
    }
    //############################## Determination of the contacts between agrgates #######################################
    std::multimap<double, size_t> AggregatList::filter_neighborhood(const size_t moving_aggregate,
                                                                    const std::array<double, 3> &direction,
                                                                    const std::vector<size_t> &neighborhood,
                                                                    const double distance) const
    {
        std::multimap<double, size_t> sorted_neighborhood;
        Sphere sphere_me(*list[moving_aggregate]);
        Sphere sphere_other(*list[moving_aggregate]);

        //$ [For all other aggregate]
        for (const size_t &agg_other : neighborhood)
        {
            sphere_other.init_val(list[agg_other]->get_position(),
                                  *list[agg_other]->rmax);
            SphereContactInfo potential_contact = distance_to_contact(sphere_me, sphere_other, direction, distance);
            if (potential_contact.distance < distance)
            {
                // not too far
                sorted_neighborhood.insert(std::pair<double, size_t>(potential_contact.distance, agg_other));
            }
        }
        return sorted_neighborhood;
    }
    //################################### Determination of the contacts between agrgates ###################################
    bool AggregatList::test_free_space(std::array<double, 3> pos, double radius) const
    {
        // Use Verlet to reduce search
        std::vector<size_t> neighborhood = verlet.get_neighborhood(pos, radius + maxradius);
        Sphere sphere(*physicalmodel, pos, radius);
        //$ loop on the aggregates potentially in contact
        for (const size_t &suspect : neighborhood) // For every aggregate that could be in contact
        {
            if (contact(sphere, *list[suspect]))
            {
                return false;
            }
        }
        return true;
    }
    bool AggregatList::croissance_surface(double dt)
    {
        bool removed_aggregate(false);
        for (auto aggregate = list.begin(); aggregate != list.end();)
        {
            auto index = std::distance(list.begin(), aggregate);
            if (croissance_surface(dt, size_t(index)))
            {
                removed_aggregate = true;
                aggregate = list.begin() + index;
            }
            else
            {
                aggregate++;
            }
        }
        return removed_aggregate;
    }
    bool AggregatList::croissance_surface(const double dt, const size_t index)
    {
        if (list[index]->croissance_surface(dt))
        {
            remove(index);
            setpointers();
            return true;
        }
        else
        {
            return false;
        }
    }
} // namespace mcac
