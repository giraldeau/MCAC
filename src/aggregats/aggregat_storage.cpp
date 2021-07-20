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
#include "aggregats/aggregat.hpp"
#include "aggregats/aggregat_list.hpp"
#include "spheres/sphere.hpp"


namespace mcac {
void Aggregate::setpointers() noexcept {
    rg = &(*storage)[AggregatesFields::AGGREGAT_RG][index_in_storage];
    f_agg = &(*storage)[AggregatesFields::AGGREGAT_F_AGG][index_in_storage];
    lpm = &(*storage)[AggregatesFields::AGGREGAT_LPM][index_in_storage];
    time_step = &(*storage)[AggregatesFields::AGGREGAT_TIME_STEP][index_in_storage];
    rmax = &(*storage)[AggregatesFields::AGGREGAT_RMAX][index_in_storage];
    agregat_volume = &(*storage)[AggregatesFields::AGGREGAT_VOLUME][index_in_storage];
    agregat_surface = &(*storage)[AggregatesFields::AGGREGAT_SURFACE][index_in_storage];
    x = &(*storage)[AggregatesFields::AGGREGAT_X][index_in_storage];
    y = &(*storage)[AggregatesFields::AGGREGAT_Y][index_in_storage];
    z = &(*storage)[AggregatesFields::AGGREGAT_Z][index_in_storage];
    rx = &(*storage)[AggregatesFields::AGGREGAT_RX][index_in_storage];
    ry = &(*storage)[AggregatesFields::AGGREGAT_RY][index_in_storage];
    rz = &(*storage)[AggregatesFields::AGGREGAT_RZ][index_in_storage];
    proper_time = &(*storage)[AggregatesFields::AGGREGAT_TIME][index_in_storage];
    dp = &(*storage)[AggregatesFields::AGGREGAT_DP][index_in_storage];
    dg_over_dp = &(*storage)[AggregatesFields::AGGREGAT_DG_OVER_DP][index_in_storage];
    overlapping = &(*storage)[AggregatesFields::AGGREGAT_OVERLAPPING][index_in_storage];
    coordination_number = &(*storage)[AggregatesFields::AGGREGAT_COORDINATION_NUMBER][index_in_storage];
    d_m = &(*storage)[AggregatesFields::AGGREGAT_D_M][index_in_storage];
    CH_ratio = &(*storage)[AggregatesFields::AGGREGAT_CH_RATIO][index_in_storage];
}
Aggregate::Aggregate() noexcept:
    ElemStorage<AggregatesFields::AGGREGAT_NFIELDS, AggregatList>(),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    agregat_volume(nullptr),
    agregat_surface(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    proper_time(nullptr),
    dp(nullptr),
    dg_over_dp(nullptr),
    overlapping(nullptr),
    coordination_number(nullptr),
    d_m(nullptr),
    CH_ratio(nullptr),
    electric_charge(0),
    n_spheres(0),
    alpha_vs_extreme(0),
    label(0),
    bulk_density(0.),
    distances(),
    distances_center(),
    volumes(),
    surfaces(),
    physicalmodel(nullptr),
    verlet(nullptr),
    index_verlet({{0, 0, 0}}),
    myspheres() {
}
Aggregate::Aggregate(AggregatList *aggregat_list, size_t newlabel) noexcept:
    ElemStorage<AggregatesFields::AGGREGAT_NFIELDS, AggregatList>(*aggregat_list, newlabel),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    agregat_volume(nullptr),
    agregat_surface(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    proper_time(nullptr),
    dp(nullptr),
    dg_over_dp(nullptr),
    overlapping(nullptr),
    coordination_number(nullptr),
    d_m(nullptr),
    CH_ratio(nullptr),
    electric_charge(0),
    n_spheres(0),
    alpha_vs_extreme(0),
    label(newlabel),
    bulk_density(0.),
    distances(),
    distances_center(),
    volumes(),
    surfaces(),
    physicalmodel(aggregat_list->physicalmodel),
    verlet(nullptr),
    index_verlet({{0, 0, 0}}),
    myspheres(*aggregat_list->physicalmodel, 1) {
}
Aggregate::~Aggregate() noexcept {
    unset_verlet();
}
/** Copy constructor */
Aggregate::Aggregate(const Aggregate &other, AggregatList *aggregat_list, size_t newlabel) noexcept:
    ElemStorage<AggregatesFields::AGGREGAT_NFIELDS, AggregatList>(*aggregat_list, newlabel),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    agregat_volume(nullptr),
    agregat_surface(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    proper_time(nullptr),
    dp(nullptr),
    dg_over_dp(nullptr),
    overlapping(nullptr),
    coordination_number(nullptr),
    d_m(nullptr),
    CH_ratio(nullptr),
    electric_charge(0),
    n_spheres(other.n_spheres),
    label(other.label),
    bulk_density(other.bulk_density),
    distances(other.distances),
    distances_center(other.distances_center),
    volumes(other.volumes),
    surfaces(other.surfaces),
    physicalmodel(other.physicalmodel),
    verlet(nullptr),
    index_verlet({{0, 0, 0}}),
    myspheres(&aggregat_list->spheres, {}) {
    for (const auto& sphere : other.myspheres) {
        myspheres.list.push_back(aggregat_list->spheres.add(*sphere, aggregat_list->spheres));
    }
    setpointers();
    for (const auto& s : myspheres) {
        s->set_label(long(index_in_storage));
    }
}
//Aggregate::Aggregate(const Aggregate &other) :
//    storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat>(other),
//    StatisicsData(other),
//    physicalmodel(other.physicalmodel),
//    myspheres(other.myspheres),
//    verlet(nullptr),
//    IndexVerlet({{0, 0, 0}}),
//    _distances(other._distances),
//    distances_center(other.distances_center),
//    volumes(other.volumes),
//    surfaces(other.surfaces),
//    rg(nullptr),
//    f_agg(nullptr),
//    lpm(nullptr),
//    time_step(nullptr),
//    rmax(nullptr),
//    volAgregat(nullptr),
//    surfAgregat(nullptr),
//    x(nullptr),
//    y(nullptr),
//    z(nullptr),
//    rx(nullptr),
//    ry(nullptr),
//    rz(nullptr),
//    time(nullptr),
//    Np(other.Np),
//    Label(other.Label),
//    InVerlet(false) {
//    setpointers();
//}
/** Move constructor */
//Aggregate::Aggregate(Aggregate &&other)
//noexcept: /* noexcept needed to enable optimizations in containers */
//    storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat>(move(other)),
//    StatisicsData(move(other)),
//    myspheres(move(other.myspheres)),
//    verlet(nullptr),
//    index_verlet({{0, 0, 0}}),
//    distances(),
//    distances_center(),
//    volumes(),
//    surfaces(),
//    rg(nullptr),
//    f_agg(nullptr),
//    lpm(nullptr),
//    time_step(nullptr),
//    rmax(nullptr),
//    agregat_volume(nullptr),
//    agregat_surface(nullptr),
//    x(nullptr),
//    y(nullptr),
//    z(nullptr),
//    rx(nullptr),
//    ry(nullptr),
//    rz(nullptr),
//    time(nullptr),
//    n_spheres(other.n_spheres),
//    label(other.label),
//    physicalmodel(other.physicalmodel) {
//    swap(distances, other.distances);
//    swap(distances_center, other.distances_center);
//    swap(volumes, other.volumes);
//    swap(surfaces, other.surfaces);
//    setpointers();
//}


/** Copy assignment operator */
/*
Aggregate& Aggregate::operator= (const Aggregate& other)
{
    Aggregate tmp(other);      // re-use copy-constructor
    *this = std::move(tmp); // re-use move-assignment
    return *this;
}
*/
/** Move assignment operator */
//Aggregate &Aggregate::operator=(Aggregate &&other) noexcept {
//    if (physicalmodel && physicalmodel->toBeDestroyed) {
//        delete physicalmodel;
//    }
//    physicalmodel = other.physicalmodel;
//    other.physicalmodel = nullptr;
//    swap(myspheres, other.myspheres);
//    swap(_distances, other._distances);
//    swap(distances_center, other.distances_center);
//    swap(volumes, other.volumes);
//    swap(surfaces, other.surfaces);
//    swap(IndexVerlet, other.IndexVerlet);
//    swap(Label, other.Label);
//    swap(Np, other.Np);
//    StatisicsData::operator=(static_cast<StatisicsData>(other));
//    storage_elem<AggregatesFields::AGGREGAT_NFIELDS,
//                 ListAggregat>::operator=(move(static_cast<storage_elem<AggregatesFields::AGGREGAT_NFIELDS,
//                                                                        ListAggregat> >(other)));
//    setpointers();
//    return *this;
//}
}  // namespace mcac
