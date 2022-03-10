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
#include "aggregats/aggregat_list.hpp"
#include "tools/tools.hpp"
#include "spheres/sphere.hpp"
#include <iostream>


namespace mcac {
void AggregatList::setpointers() {
    auto newdeb((*storage)[0].begin());
    auto newfin((*storage)[0].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin)) {
        return;
    }
    for (const auto& aggregate : list) {
        aggregate->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}
void AggregatList::remove(const size_t &id) noexcept {
    ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::remove(id);
    // keep index and label in sync
    for (size_t i = id; i < list.size(); i++) {
        list[i]->decrease_label();
    }
    setpointers();
}
void AggregatList::remove_sphere(const size_t &id) noexcept {
    auto agg_num = size_t(spheres[id]->agg_label);
    list[agg_num]->remove_sphere(id);
    if (list[agg_num]->size() == 0) {
        remove(agg_num);
    }
}
AggregatList::AggregatList(PhysicalModel *the_physical_model):
    ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>(),
    physicalmodel(the_physical_model),
    maxradius(0.),
    avg_npp(1),
    max_time_step(0.),
    index_sorted_time_steps(),
    cumulative_time_steps(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(std::make_unique<ThreadedIO>(physicalmodel->output_dir / "Aggregats", *physicalmodel, physicalmodel->n_monomeres)),
    last_saved(0),
    spheres(*the_physical_model, physicalmodel->n_monomeres),
    verlet(the_physical_model->n_verlet_divisions, the_physical_model->box_length) {
    ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::init(physicalmodel->n_monomeres, *this);
    setpointers();

    //Initialize the aggregates
    for (size_t i = 0; i < size(); i++) {
        list[i]->init(i, i, false);
    }
    refresh();

    if (physicalmodel->enforce_volume_fraction){
        double current_total_volume = get_total_volume();
        double prescribed_total_volume = physicalmodel->volume_fraction *  std::pow(physicalmodel->box_length, 3);
        double correction = std::pow(prescribed_total_volume / current_total_volume, 1./3.);

        for (const auto & sphere: spheres) {
            sphere->set_radius(sphere->get_radius() * correction);
            sphere->update_vol_and_surf();
        }
        for (const auto & aggregate: list) {
            aggregate->compute_volume_surface();
        }
    }
}
AggregatList::~AggregatList() noexcept {
    //#pragma omp simd
    for (const auto& aggregate : list) {
        aggregate->unset_verlet();
    }
}
}// namespace mcac
