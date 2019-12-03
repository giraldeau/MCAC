#include "aggregats/aggregat.hpp"
#include "aggregats/aggregat_list.hpp"
#include "spheres/sphere.hpp"


using namespace std;
namespace MCAC {
void Aggregate::setpointers() noexcept {
    rg = &(*Storage)[AggregatesFields::AGGREGAT_RG][indexInStorage];
    f_agg = &(*Storage)[AggregatesFields::AGGREGAT_F_AGG][indexInStorage];
    lpm = &(*Storage)[AggregatesFields::AGGREGAT_LPM][indexInStorage];
    time_step = &(*Storage)[AggregatesFields::AGGREGAT_TIME_STEP][indexInStorage];
    rmax = &(*Storage)[AggregatesFields::AGGREGAT_RMAX][indexInStorage];
    agregat_volume = &(*Storage)[AggregatesFields::AGGREGAT_VOLUME][indexInStorage];
    agregat_surface = &(*Storage)[AggregatesFields::AGGREGAT_SURFACE][indexInStorage];
    x = &(*Storage)[AggregatesFields::AGGREGAT_X][indexInStorage];
    y = &(*Storage)[AggregatesFields::AGGREGAT_Y][indexInStorage];
    z = &(*Storage)[AggregatesFields::AGGREGAT_Z][indexInStorage];
    rx = &(*Storage)[AggregatesFields::AGGREGAT_RX][indexInStorage];
    ry = &(*Storage)[AggregatesFields::AGGREGAT_RY][indexInStorage];
    rz = &(*Storage)[AggregatesFields::AGGREGAT_RZ][indexInStorage];
    time = &(*Storage)[AggregatesFields::AGGREGAT_TIME][indexInStorage];
}
Aggregate::Aggregate() noexcept:
    storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat>(),
    StatisicsData(),
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
    time(nullptr),
    n_spheres(0),
    label(0),
    distances(),
    distances_center(),
    volumes(),
    surfaces(),
    physicalmodel(nullptr),
    verlet(nullptr),
    index_verlet({{0, 0, 0}}),
    myspheres() {
}
Aggregate::Aggregate(ListAggregat &_storage, const size_t newlabel) noexcept:
    storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat>(_storage, newlabel),
    StatisicsData(),
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
    time(nullptr),
    n_spheres(0),
    label(newlabel),
    distances(),
    distances_center(),
    volumes(),
    surfaces(),
    physicalmodel(_storage.physicalmodel),
    verlet(nullptr),
    index_verlet({{0, 0, 0}}),
    myspheres(*_storage.physicalmodel, 1) {
}
Aggregate::~Aggregate() noexcept {
    if (static_cast<bool>(verlet)) {
        verlet->Remove(get_label(), index_verlet);
    }
}
/** Copy constructor */
Aggregate::Aggregate(const Aggregate &other, ListAggregat &_Storage) noexcept:
    storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat>(other, *this, _Storage),
    StatisicsData(other),
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
    time(nullptr),
    n_spheres(other.n_spheres),
    label(other.label),
    distances(other.distances),
    distances_center(other.distances_center),
    volumes(other.volumes),
    surfaces(other.surfaces),
    physicalmodel(other.physicalmodel),
    verlet(nullptr),
    index_verlet({{0, 0, 0}}),
    myspheres(other.myspheres, _Storage.spheres) {
    setpointers();
    for (Sphere *s : myspheres) {
        s->set_label(long(indexInStorage));
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
//    StatisicsData::operator=(static_cast<StatisicsData &>(other));
//    storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat>::operator=(move(static_cast<storage_elem<AggregatesFields::AGGREGAT_NFIELDS, ListAggregat> &>(other)));
//    setpointers();
//    return *this;
//}
}  // namespace MCAC
