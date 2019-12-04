#include "aggregats/aggregat_list.hpp"
#include "aggregats/aggregat_distance.hpp"
#include "spheres/sphere_collision.hpp"
#include "spheres/sphere_distance.hpp"
#include "tools/tools.hpp"
#include <iostream>

using namespace std;
namespace MCAC {
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
    double valAlea = Random() * cumulative_time_steps[size() - 1];
    size_t n = lower_bound(cumulative_time_steps.begin(), cumulative_time_steps.end(), valAlea)
               - cumulative_time_steps.begin();
    size_t NumAgg = index_sorted_time_steps[n];
    return NumAgg;
}
[[gnu::pure]]  size_t AggregatList::pick_last() const {
    // TODO cache result
    double time = *list[0]->time;
    size_t latest = 0;
    for (const Aggregate *Agg : list) {
        if (*Agg->time < time) {
            time = *Agg->time;
            latest = Agg->label;
        }
    }
    return latest;
}
void AggregatList::refresh() {
    max_time_step = *list[0]->time_step;
    for (const Aggregate *Agg : list) {
        max_time_step = MAX(*Agg->time_step, max_time_step);
    }
}
template<typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {
             return v[i1] < v[i2];
         });
    return idx;
}
void AggregatList::sort_time_steps(double factor) {
    vector<double> TpT(size());
#pragma omp simd
    for (size_t i = 0; i < size(); i++) {
        TpT[i] = factor / (*list[i]->time_step);
    }
    index_sorted_time_steps = sort_indexes(TpT);
    cumulative_time_steps.resize(size());

    //$ Accumulate the timesteps
    cumulative_time_steps[0] = TpT[index_sorted_time_steps[0]];
    for (size_t i = 1; i < size(); i++) {
        cumulative_time_steps[i] = cumulative_time_steps[i - 1] + TpT[index_sorted_time_steps[i]];
    }
}
void AggregatList::duplication() {
    size_t oldNAgg = size();
    double oldL = physicalmodel->L;
    physicalmodel->L *= 2;
    physicalmodel->N *= 8;

    // TODO Rework this in order not to to it aggregate by aggregate

    for (size_t iagg = 0; iagg < oldNAgg; iagg++) {
        for (int i = 0; i <= 1; i++) {
            for (int j = 0; j <= 1; j++) {
                for (int k = 0; k <= 1; k++) {
                    if (i != 0 or j != 0 or k != 0) {
                        Aggregate *newAgg =
                            storage_list<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::add(*list[iagg], *this);
                        newAgg->label = size() - 1;
                        newAgg->set_verlet(verlet);
                        setpointers();
                        for (Aggregate *Agg : list) {
                            Agg->setpointers();
                            for (Sphere *Sph : Agg->myspheres) {
                                Sph->setpointers();
                            }
                        }
                        newAgg->set_verlet(verlet); // TODO remove?
                        array<double, 3> vec_move = {i * oldL,
                                                     j * oldL,
                                                     k * oldL};
                        newAgg->translate(vec_move);
                    }
                }
            }
        }
    }
    //$ update Verlet
    verlet.Init(physicalmodel->GridDiv, physicalmodel->L);
    for (Aggregate *Agg : list) {
        Agg->update_verlet_index();
    }
}
size_t AggregatList::merge(size_t first, size_t second) {
    const size_t keeped(MIN(first, second));
    const size_t removed(MAX(first, second));

    // compute proper time of the final aggregate
    // keeping global time constant
    double newtime = double(size() - 1) * (*list[keeped]->time + *list[removed]->time) / double(size())
                     - physicalmodel->Time;

    // compute the new average of npp
    avg_npp = avg_npp * size() / (static_cast<double>(size()) - 1);

    // merge the two aggregate but do not remove the deleted one
    list[keeped]->merge(*list[removed]);

    remove(removed);
    setpointers();
    *list[keeped]->time = newtime;
    return keeped;
}
//################################# Determination of the contacts between agrgates ####################################
pair<size_t, double> AggregatList::distance_to_next_contact(size_t source, array<double, 3> direction) const {
    // Use Verlet to reduce search
    vector<size_t> SearchSpace(get_search_space(source, direction));

    // Assimilate Aggregate as sphere to drasticly speed-up search
    vector<pair<size_t, double> > Potential(sort_search_space(source, direction, SearchSpace));
    int aggcontact(-1);
    double distmin = list[source]->get_lpm();
    double mindistagg(0.);

    //$ loop on the agregates potentially in contact
    for (const pair<size_t, double> &suspect : Potential) //For every aggregate that could be in contact
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
            if (distmin <= 1e-15){
                break;
            }
        }
    }
    return {aggcontact, distmin};
}
vector<size_t> AggregatList::get_search_space(size_t source, array<double, 3> direction) const {
    vector<size_t> SearchSpace;
    if (physicalmodel->use_verlet) {
        // Extract from verlet

        double lpm(*list[source]->lpm);
        double mindist(*list[source]->rmax + maxradius);
        array<double, 3> sourceposition = list[source]->get_position();
        array<double, 3> Vector{lpm * direction};
        SearchSpace = verlet.GetSearchSpace(sourceposition, mindist, Vector);

        // Remove me
        for (size_t i = 0; i < SearchSpace.size(); i++) {
            if (SearchSpace[i] == source) {
                SearchSpace.erase(SearchSpace.begin() + long(i));
                return SearchSpace;
            }
        }
        cout << "I'm not on the verlet list ???" << endl;
        cout << "This is an error" << endl;
        exit(ErrorCodes::VERLET_ERROR);
        //return SearchSpace;
    }

    // The full aggregat list index
    SearchSpace.resize(size());
    iota(SearchSpace.begin(), SearchSpace.end(), 0);

    // Except me
    SearchSpace.erase(SearchSpace.begin() + long(source));
    return SearchSpace;
}
//############################## Determination of the contacts between agrgates #######################################
vector<pair<size_t, double> > AggregatList::sort_search_space(size_t moving_aggregate,
                                                              array<double, 3> direction,
                                                              const vector<size_t> &search_space) const {
    vector<pair<size_t, double> > SortedSearchSpace;
    Sphere SphereMe(*list[moving_aggregate]);
    Sphere SphereOther(*list[moving_aggregate]);

    //$ [For all other agregate]
    for (const size_t &AggOther : search_space) {
        SphereOther.init_val(list[AggOther]->get_position(),
                             *list[AggOther]->rmax);
        auto pos = SortedSearchSpace.begin();
        double dist = 0.;
        bool iscollision = contact(SphereMe, SphereOther);
        if (!iscollision) {
            pair<bool, double> Collision = collision(SphereMe, SphereOther, direction);
            if (Collision.first && Collision.second <= *list[moving_aggregate]->lpm) {
                iscollision = Collision.first;
                dist = Collision.second;
                pos = lower_bound(SortedSearchSpace.begin(), SortedSearchSpace.end(), dist,
                                  [](pair<size_t, double> i1, double d) {
                                      return i1.second < d;
                                  });
            }
        }
        if (iscollision) {
            pair<size_t, double> suspect = {AggOther, dist};
            SortedSearchSpace.insert(pos, suspect);
        }
    }
    return SortedSearchSpace;
}
//################################### Determination of the contacts between agrgates ###################################
bool AggregatList::test_free_space(array<double, 3> pos, double diameter) const {
    // Use Verlet to reduce search
    double mindist(diameter * 0.5 + maxradius);
    vector<size_t> SearchSpace = verlet.GetSearchSpace(pos, mindist);
    double mindist2 = POW_2(mindist);

    //$ loop on the agregates potentially in contact
    for (const size_t &suspect : SearchSpace) //For every aggregate that could be in contact
    {
        //$ Loop on all the spheres of the aggregate
        for (const Sphere *sphere : list[suspect]->myspheres) {
            if (distance_2(sphere->get_position(), pos, sphere->physicalmodel->L) <= mindist2) {
                return false;
            }
        }
    }
    return true;
}
}// namespace MCAC

