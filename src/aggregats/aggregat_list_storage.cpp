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
#include <iostream>


namespace mcac {
void AggregatList::setpointers() {
    auto newdeb((*storage)[0].begin());
    auto newfin((*storage)[0].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin)) {
        return;
    }
    for (Aggregate *aggregate : list) {
        aggregate->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}
AggregatList::AggregatList(PhysicalModel *the_physical_model) noexcept:
    ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>(),
    physicalmodel(the_physical_model),
    maxradius(0.),
    avg_npp(1),
    max_time_step(0.),
    index_sorted_time_steps(),
    cumulative_time_steps(),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    writer(new ThreadedIO(*physicalmodel, physicalmodel->n_monomeres)),
    last_saved(0),
    spheres(*the_physical_model, physicalmodel->n_monomeres),
    verlet(the_physical_model->n_verlet_divisions, the_physical_model->box_lenght) {
    ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::init(physicalmodel->n_monomeres, *this);
    setpointers();

    //Initialize the aggregates
    for (size_t i = 0; i < size(); i++) {
        list[i]->init(i, i);
    }
    refresh();
}
AggregatList::~AggregatList() noexcept {
    delete writer;
    writer = nullptr;

    //#pragma omp simd
    for (Aggregate *aggregate : list) {
        aggregate->unset_verlet();
    }
}
}// namespace mcac

