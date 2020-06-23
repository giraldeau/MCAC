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
#include "spheres/sphere_collision.hpp"
#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include "exceptions.hpp"
#include <iostream>


namespace mcac {
[[gnu::pure]] double AggregatList::get_avg_npp() const {
    return avg_npp;
}
[[gnu::pure]] double AggregatList::get_max_time_step() const {
    return max_time_step;
}
[[gnu::pure]] double AggregatList::get_time_step(double max) const {
    double deltatemps = max / cumulative_time_steps[size() - 1];
    return deltatemps;
}
[[gnu::pure]] size_t AggregatList::pick_random() const {
    //$ Pick a random sphere
    double val_alea = random() * cumulative_time_steps[size() - 1];
    long n = lower_bound(cumulative_time_steps.begin(), cumulative_time_steps.end(), val_alea)
             - cumulative_time_steps.begin();
    size_t num_agg = index_sorted_time_steps[static_cast<size_t >(n)];
    return num_agg;
}
[[gnu::pure]]  size_t AggregatList::pick_last() const {
    // TODO(pouxa): cache result
    double time = *list[0]->time;
    size_t latest = 0;
    for (const Aggregate *agg : list) {
        if (*agg->time < time) {
            time = *agg->time;
            latest = agg->label;
        }
    }
    return latest;
}
void AggregatList::add(size_t n) {
    size_t initial_n_agg = size();
    size_t initial_n_sph = spheres.size();
    spheres.add(n, spheres);
    spheres.setpointers();
    add(n, *this);
    setpointers();

    //initialize data
    for (size_t i = 0; i < n; i++) {
        size_t i_agg = initial_n_agg + i;
        size_t i_sph = initial_n_sph + i;
        list[i_agg]->init(i_agg, i_sph);
    }
    refresh();
}
void AggregatList::refresh() {
    max_time_step = *list[0]->time_step;
    for (const Aggregate *agg : list) {
        max_time_step = MAX(*agg->time_step, max_time_step);
    }
    avg_npp = static_cast<double>(spheres.size()) / static_cast<double>(size());
}
template<typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i_1, size_t i_2) {
             return v[i_1] < v[i_2];
         });
    return idx;
}
void AggregatList::sort_time_steps(double factor) {
    std::vector<double> tp_t(size());
#pragma omp simd
    for (size_t i = 0; i < size(); i++) {
        tp_t[i] = factor / (*list[i]->time_step);
    }
    index_sorted_time_steps = sort_indexes(tp_t);
    cumulative_time_steps.resize(size());

    //$ Accumulate the timesteps
    cumulative_time_steps[0] = tp_t[index_sorted_time_steps[0]];
    for (size_t i = 1; i < size(); i++) {
        cumulative_time_steps[i] = cumulative_time_steps[i - 1] + tp_t[index_sorted_time_steps[i]];
    }
}
void AggregatList::duplication() {
    size_t old_n_agg = size();
    double old_l = physicalmodel->box_lenght;
    physicalmodel->box_lenght *= 2;
    physicalmodel->n_monomeres *= 8;

    // TODO(pouxa): Rework this in order not to to it aggregate by aggregate

    for (size_t iagg = 0; iagg < old_n_agg; iagg++) {
        for (int i = 0; i <= 1; i++) {
            for (int j = 0; j <= 1; j++) {
                for (int k = 0; k <= 1; k++) {
                    if (i != 0 or j != 0 or k != 0) {
                        Aggregate *new_agg =
                            ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::add(*list[iagg], *this);
                        new_agg->label = size() - 1;
                        new_agg->set_verlet(&verlet);
                        setpointers();
                        spheres.setpointers();
                        new_agg->set_verlet(&verlet); // TODO(pouxa): remove?
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
    verlet = Verlet(physicalmodel->n_verlet_divisions, physicalmodel->box_lenght);
    for (Aggregate *agg : list) {
        agg->update_verlet_index();
    }
}
size_t AggregatList::merge(size_t first, size_t second) {
    const size_t _keeped(MIN(first, second));
    const size_t _removed(MAX(first, second));

    // compute proper time of the final aggregate
    // keeping global time constant
    double newtime = double(size() - 1) * (*list[_keeped]->time + *list[_removed]->time) / double(size())
                     - physicalmodel->time;

    // merge the two aggregate but do not remove the deleted one
    list[_keeped]->merge(list[_removed]);
    remove(_removed);
    setpointers();
    *list[_keeped]->time = newtime;
    return _keeped;
}
bool AggregatList::split() {
    bool has_splitted = false;
    auto suspect = list.begin();
    while (suspect != list.end()) {
        auto index = std::distance(list.begin(), suspect);
        // the split function will create the new aggregates
        if ((*suspect)->split()) {
            has_splitted = true;
            // but we still have to destroy the current one
            remove(size_t(index));
            setpointers();
            suspect = list.begin() + index;
        } else {
            suspect = list.begin() + index + 1;
        }
    }
    return has_splitted;
}
//################################# Determination of the contacts between agrgates ####################################
std::pair<int, double> AggregatList::distance_to_next_contact(size_t source, std::array<double, 3> direction) const {
    // Use Verlet to reduce search
    std::vector<size_t> search_space(get_search_space(source, direction));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    std::vector<std::pair<size_t, double>> potential(sort_search_space(source, direction, search_space));
    int aggcontact(-1);
    double distmin = list[source]->get_lpm();
    double mindistagg(0.);

    //$ loop on the agregates potentially in contact
    for (const std::pair<size_t, double> &suspect : potential) //For every aggregate that could be in contact
    {
        size_t agg = suspect.first;
        double distagg = suspect.second;

        // We already found the closest one
        if (aggcontact >= 0) {
            double secu = 2 * (*list[size_t(aggcontact)]->rmax + *list[source]->rmax);
            if (distagg - mindistagg > secu) {
                break;
            }
        }
        double dist = distance(*list[source], *list[agg], direction);
        if (dist >= 0 && dist <= distmin) {
            distmin = dist;
            aggcontact = int(agg);
            mindistagg = distagg;

            // If two aggragates are already in contact (due to surface growing)
            if (distmin <= 1e-28) {
                break;
            }
        }
    }
    return {aggcontact, distmin};
}
std::vector<size_t> AggregatList::get_search_space(size_t source, std::array<double, 3> direction) const {
    std::vector<size_t> search_space;

    // Extract from verlet

    double lpm(*list[source]->lpm);
    double mindist(*list[source]->rmax + maxradius);
    std::array<double, 3> sourceposition = list[source]->get_position();
    std::array<double, 3> vector{lpm * direction};
    search_space = verlet.get_search_space(sourceposition, mindist, vector);

    // Remove me
    for (size_t i = 0; i < search_space.size(); i++) {
        if (search_space[i] == source) {
            search_space.erase(search_space.begin() + long(i));
            return search_space;
        }
    }
    throw VerletError("Aggregate not on the verlet list ???");
    //return SearchSpace;
}
//############################## Determination of the contacts between agrgates #######################################
std::vector<std::pair<size_t, double>> AggregatList::sort_search_space(size_t moving_aggregate,
                                                                       std::array<double, 3> direction,
                                                                       const std::vector<size_t> &search_space) const {
    std::vector<std::pair<size_t, double>> sorted_search_space;
    Sphere sphere_me(*list[moving_aggregate]);
    Sphere sphere_other(*list[moving_aggregate]);

    //$ [For all other agregate]
    for (const size_t &agg_other : search_space) {
        sphere_other.init_val(list[agg_other]->get_position(),
                              *list[agg_other]->rmax);
        auto pos = sorted_search_space.begin();
        AggregateCollisionInfo collision;
        if (contact(sphere_me, sphere_other)) {
            collision = AggregateCollisionInfo(moving_aggregate, agg_other, 0.);
            pos = sorted_search_space.begin();
        } else {
            SphereCollisionInfo potential_collision = sphere_collision(sphere_me, sphere_other, direction);
            if (potential_collision.is_collision && potential_collision.distance <= *list[moving_aggregate]->lpm) {
                collision = AggregateCollisionInfo(moving_aggregate, agg_other, potential_collision.distance);
                pos = lower_bound(sorted_search_space.begin(), sorted_search_space.end(), potential_collision);
            }
        }
        if (collision.is_collision) {
            sorted_search_space.insert(pos, collision);
        }
    }
    return sorted_search_space;
}
//################################### Determination of the contacts between agrgates ###################################
bool AggregatList::test_free_space(std::array<double, 3> pos, double radius) const {
    // Use Verlet to reduce search
    std::vector<size_t> search_space = verlet.get_search_space(pos, radius + maxradius);
    Sphere sphere(*physicalmodel, pos, radius);
    //$ loop on the agregates potentially in contact
    for (const size_t &suspect : search_space) //For every aggregate that could be in contact
    {
        if (contact(sphere, *list[suspect])) {
            return false;
        }
    }
    return true;
}
void AggregatList::croissance_surface(double dt) {
    spheres.croissance_surface(dt);
    for (Sphere *sphere : spheres) {
        if (sphere->get_radius() <= 0) {
            remove_sphere(sphere->get_index());
        }
    }
}
}// namespace mcac
