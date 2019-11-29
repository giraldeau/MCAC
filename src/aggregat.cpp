#include "aggregat.hpp"
#include "sblvolumesurface.hpp"
#include "spheres/sphere_collision.hpp"
#include "spheres/sphere_distance.hpp"
#include "spheres/sphere_intersection.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <functional>
#include <iostream>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))
using namespace std;
namespace MCAC {
const double PI = atan(1.0) * 4;
Aggregate::Aggregate() :
    storage_elem<15, ListAggregat>(),
    StatisicsData(),
    physicalmodel(nullptr),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0, 0, 0}}),
    _distances(),
    distances_center(),
    volumes(),
    surfaces(),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    time(nullptr),
    Np(0),
    Label(0),
    InVerlet(false) {
    Init();
}
Aggregate::Aggregate(PhysicalModel &_physicalmodel) :
    storage_elem<15, ListAggregat>(),
    StatisicsData(),
    physicalmodel(&_physicalmodel),
    myspheres(),
    verlet(nullptr),
    IndexVerlet({{0, 0, 0}}),
    _distances(),
    distances_center(),
    volumes(),
    surfaces(),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    time(nullptr),
    Np(0),
    Label(0),
    InVerlet(false) {
    Init();
}
Aggregate::Aggregate(ListAggregat &_storage, const size_t label) :
    storage_elem<15, ListAggregat>(_storage, label),
    StatisicsData(),
    physicalmodel(_storage.physicalmodel),
    myspheres(*_storage.physicalmodel, 1),
    verlet(nullptr),
    IndexVerlet({{0, 0, 0}}),
    _distances(),
    distances_center(),
    volumes(),
    surfaces(),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    time(nullptr),
    Np(0),
    Label(label),
    InVerlet(false) {
    Init();
}
void Aggregate::Init() {
    if (external_storage) {
        external_storage->setpointers();
    }
    setpointers();
    *rg = 0.;              //Gyration Radius
    *f_agg = 0.;              //Friction coeff
    *lpm = 0.;              //Aggregate apparent Mean Free Path
    *time_step = 0.;              //Time to move along lpm
    *rmax = 0.;              //Radius of the sphere containing Agg
    *volAgregat = 0.;              //Etimation of the Aggregate's volume
    *surfAgregat = 0.;              //Estimation of the sufrace of the aggregate

    *rx = 0;
    *ry = 0;
    *rz = 0;
    *time = 0;
    IndexVerlet = {{0, 0, 0}};
    SetPosition(0., 0., 0.);
    UpdateDistances();
}
void Aggregate::setpointers() {
    size_t myindex = indexInStorage;
    rg = &(*Storage)[0][myindex];               //Gyration Radius
    f_agg = &(*Storage)[1][myindex];               //Friction coeff
    lpm = &(*Storage)[2][myindex];               //Aggregate apparent Mean Free Path
    time_step = &(*Storage)[3][myindex];               //Time to move along lpm
    rmax = &(*Storage)[4][myindex];               //Radius of the sphere containing Agg
    volAgregat = &(*Storage)[5][myindex];               //Etimation of the Aggregate's volume
    surfAgregat = &(*Storage)[6][myindex];               //Estimation of the surface of the aggregate
    x = &(*Storage)[7][myindex];
    y = &(*Storage)[8][myindex];
    z = &(*Storage)[9][myindex];
    rx = &(*Storage)[10][myindex];
    ry = &(*Storage)[11][myindex];
    rz = &(*Storage)[12][myindex];
    time = &(*Storage)[13][myindex];
}
void Aggregate::Init(PhysicalModel &_physicalmodel,
                     Verlet &_verlet,
                     const array<double, 3> position,
                     const size_t _label,
                     ListSphere &spheres, const double D) {
    if (physicalmodel && physicalmodel->toBeDestroyed) {
        delete physicalmodel;
    }
    physicalmodel = &_physicalmodel;
    verlet = &_verlet;
    Label = _label;
    indexInStorage = _label;
    if (external_storage) {
        external_storage->setpointers();
    }
    setpointers();
    if (physicalmodel->use_verlet) {
        verlet->Add(_label, IndexVerlet);
        InVerlet = true;
    }
    SetPosition(position);
    spheres[_label].set_label(int(_label));
    spheres[_label].init_val(position, D / 2);
    myspheres = ListSphere(spheres, {_label});
    Np = myspheres.size();
    UpdateDistances();
    Update();
    fullStatistics();
}
//####################################### Calcul de la distance inter-agrégats ########################################
bool Aggregate::Contact(Aggregate &other) const noexcept {
    /*
        Sphere SphereMe(*this);
        Sphere SphereOther(other);

        if (! SphereMe.contact(SphereOther))
            return false;
    */
    //$ Loop on all the spheres of the other aggregate
    for (const Sphere *othersphere : other.myspheres) {
        //$ For every sphere in the aggregate :
        for (const Sphere *mysphere: myspheres) {
            if (contact(*mysphere, *othersphere)) {
                return true;
            }
        }
    }
    return false;
}
double Aggregate::Distance(Aggregate &other, array<double, 3> Vectdir) const {
    double mindist(*lpm);
    bool contact(false);

    //$ For every sphere in the aggregate :
    for (const Sphere *mysphere: myspheres) {
        vector<double> dists = collisions(*mysphere, other.myspheres, Vectdir);
        if (!dists.empty()) {
            double lmindist = *min_element(dists.begin(), dists.end());
            if (lmindist <= mindist) {
                mindist = lmindist;
                contact = true;
            }
        }
    }
    if (contact) {
        return mindist;
    }
    return -1.;
}

//###################################################################################################################

__attribute__((pure)) array<double, 3> Aggregate::GetPosition() const noexcept {
    array<double, 3> mypos{{*x, *y, *z}};
    return mypos;
}
void Aggregate::SetPosition(const double newx, const double newy, const double newz) noexcept {
    *x = periodicPosition(newx, physicalmodel->L);
    *y = periodicPosition(newy, physicalmodel->L);
    *z = periodicPosition(newz, physicalmodel->L);
    if (InVerlet && physicalmodel->use_verlet) {
        //$ Update Verlet
        array<size_t, 3> newindexVerlet = GetVerletIndex();
        if (newindexVerlet != IndexVerlet) {
            verlet->Remove(GetLabel(), IndexVerlet);
            IndexVerlet = newindexVerlet;
            verlet->Add(GetLabel(), IndexVerlet);
        }
    }
}
void Aggregate::SetPosition(const array<double, 3> position) noexcept {
    SetPosition(position[0], position[1], position[2]);
}
void Aggregate::Translate(const array<double, 3> vector) noexcept {
    // keep track of the previous position
    array<double, 3> oldpos = GetPosition();

    // move the aggregate
    SetPosition(*x + vector[0], *y + vector[1], *z + vector[2]);

    // Real movement of the aggregate (periodic)
    array<double, 3> newpos = GetPosition();
    newpos[0] -= oldpos[0];
    newpos[1] -= oldpos[1];
    newpos[2] -= oldpos[2];

    // move the first sphere taking care of the periodicity
    *myspheres[0].x = periodicPosition(*myspheres[0].x + newpos[0], physicalmodel->L);
    *myspheres[0].y = periodicPosition(*myspheres[0].y + newpos[1], physicalmodel->L);
    *myspheres[0].z = periodicPosition(*myspheres[0].z + newpos[2], physicalmodel->L);

    // move all the other sphere relatively to the first
    array<double, 3> refpos = myspheres[0].get_position();
    for (Sphere *mysphere : myspheres) {
        newpos = mysphere->get_relative_position();
        newpos[0] += refpos[0];
        newpos[1] += refpos[1];
        newpos[2] += refpos[2];
        mysphere->set_position(newpos);
    }
}
void Aggregate::TimeForward(double deltatemps) noexcept {
    *time += deltatemps;
}
Aggregate::~Aggregate() noexcept {
    if (InVerlet && physicalmodel->use_verlet) {
        verlet->Remove(GetLabel(), IndexVerlet);
    }
    if (physicalmodel && physicalmodel->toBeDestroyed) {
        delete physicalmodel;
    }
}
__attribute__((pure)) double Aggregate::GetLpm() const noexcept {
    return *lpm;
}
__attribute__((pure)) double Aggregate::GetVolAgregat() const noexcept {
    return *volAgregat;
}
__attribute__((pure)) size_t Aggregate::GetLabel() const noexcept {
    return Label;
}
__attribute__((pure)) double Aggregate::GetTimeStep() const noexcept {
    return *time_step;
}
//### Mise à jour des paramètres physiques d'un agrégat (rayon de giration, masse, nombre de sphérules primaires) #####
void Aggregate::Update() {
    // This function will update the parameter of Agg
    Volume();
    MassCenter();
    CalcRadius();
    RayonGiration();

    //$ Determination of the friction coefficient
    *f_agg = physicalmodel->friction_coeff(Np);
    //*f_agg = physicalmodel->friction_coeff2(*rg);

    if (physicalmodel->ActiveModulephysique) {
        double masse = physicalmodel->Rho * (*volAgregat);    //Determination of the real mass of Agg
        double diff = physicalmodel->diffusivity(*f_agg);
        double rel_time = physicalmodel->relax_time(masse, *f_agg);
        *time_step = 3.0 * rel_time;
        *lpm = sqrt(6.0 * diff * (*time_step));
    } else {
        *lpm = physicalmodel->Dpm * 1E-9;
        *time_step = 1E-6;
    }
    partialStatistics();
    if (external_storage) {
        if (*rmax > external_storage->maxradius) {
            external_storage->maxradius = *rmax;
        }
    }
}
//####### Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ########
void Aggregate::Volume() {
    *volAgregat = *surfAgregat = 0.0; // Volume and surface of Agg Id

    //$ Initialisation of the arrays of volume, surface of each sphere, and the center of mass
    volumes.resize(Np);
    surfaces.resize(Np);

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < Np; i++) {
        //$ Calculation of the volume and surface of monomere i of Agg id
        volumes[i] = *myspheres[i].volume;      //Calculation of the volume of i
        surfaces[i] = *myspheres[i].surface;    //Calculation of the surface of i
    }
#ifdef WITH_SBL
    pair<vector<double>, vector<double>> volume_surface(compute_volume_surface(myspheres));

    volumes = volume_surface.first;
    surfaces = volume_surface.second;

    for (size_t i = 0; i < Np; i++)
    {
        *volAgregat = *volAgregat + volumes[i];    //Total Volume of Agg id
        *surfAgregat = *surfAgregat + surfaces[i]; //Total Surface of Agg id
    }
#else
    for (size_t i = 0; i < Np; i++) {
        for (size_t j = i + 1; j < Np; j++) //for the j spheres composing Aggregate n°id
        {
            double voli, volj, surfi, surfj;
            //            voli = volj = surfi = surfj = 0.;

            //$ Calculation of the intersection between the spheres i and j if in contact
            if (SphereDistance(i, j) > 0) {
                //                intersection(myspheres[i], myspheres[j], SphereDistance(i, j), voli, volj, surfi, surfj);
                //            }
                Intersection intersection = Intersection(myspheres[i], myspheres[j], SphereDistance(i, j));
                voli = intersection.volume_1;
                volj = intersection.volume_2;
                surfi = intersection.surface_1;
                surfj = intersection.surface_2;
            } else {
                voli = volj = surfi = surfj = 0.;
            }

            //$ The volume and surface covered by j is substracted from those of i
            volumes[i] = volumes[i] - voli;     //Calcul du volume de la sphérule i moins le volume de
            //la calotte due à la sphérule j
            surfaces[i] = surfaces[i] - surfi;  //Calcul de la surface de la sphérule i moins la surface de
            //la calotte due à la sphérule j

            //$ The volume and surface covered by i is substracted from those of j
            volumes[j] = volumes[j] - volj;     //Calcul du volume de la sphérule j moins le volume de
            //la calotte due à la sphérule i
            surfaces[j] = surfaces[j] - surfj;  //Calcul de la surface de la sphérule j moins la surface de
            //la calotte due à la sphérule i
        }
        //$ Calculation of the total volume and surface of the aggregate
        *volAgregat = *volAgregat + volumes[i];    //Total Volume of Agg id
        *surfAgregat = *surfAgregat + surfaces[i]; //Total Surface of Agg id
    }
#endif
}
void Aggregate::MassCenter() {
    const size_t loopsize(Np);
    double _rx(0.);
    double _ry(0.);
    double _rz(0.);

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < loopsize; i++) {
        //$ Calculation of the position of the center of mass
        _rx += *myspheres[i].rx * volumes[i];
        _ry += *myspheres[i].ry * volumes[i];
        _rz += *myspheres[i].rz * volumes[i];
    }
    _rx /= *volAgregat;
    _ry /= *volAgregat;
    _rz /= *volAgregat;
    for (size_t i = 0; i < loopsize; i++) {
        double dx = *myspheres[i].rx - _rx;
        double dy = *myspheres[i].ry - _ry;
        double dz = *myspheres[i].rz - _rz;
        distances_center[i] = sqrt(POW_2(dx) + POW_2(dy) + POW_2(dz));
    }
    SetPosition(*myspheres[0].x + _rx,
                *myspheres[0].y + _ry,
                *myspheres[0].z + _rz);
    *rx = _rx;
    *ry = _ry;
    *rz = _rz;
}
void Aggregate::CalcRadius() {
    // Maximum radius of the aggregate,
    // this corresponds to the distance between
    //    the center of mass of the aggregate
    //    and the edge of the furthest ball from said center.
    // It is used to assimilate the aggregate to a sphere when checking for intersections
    *rmax = 0.0;
    const size_t loopsize(Np);
    for (size_t i = 0; i < loopsize; i++) {
        *rmax = MAX(*rmax, *myspheres[i].r + SphereDistance(i));
    }
}
void Aggregate::RayonGiration() {
    // This function determines the Gyration Radius of the Aggregate Id.

    // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient,
    // they are used  used in the final formula of the Radius of Gyration
    double Arg(0.);
    double Brg(0.);
    const size_t loopsize(Np);
    for (size_t i = 0; i < loopsize; i++) {
        //$ Calculation of Rg
        Arg = Arg + volumes[i] * POW_2(SphereDistance(i)); // distance to the gravity center
        Brg = Brg + volumes[i] * POW_2(*myspheres[i].r);
    }
    *rg = sqrt(fabs((Arg + 3.0 / 5.0 * Brg) / (*volAgregat)));
    *volAgregat = fabs(*volAgregat);
}
//#####################################################################################################################


array<size_t, 3> Aggregate::GetVerletIndex() noexcept {
    array<size_t, 3> index({{0, 0, 0}});
    double step = double(physicalmodel->GridDiv) / physicalmodel->L;
    index[0] = size_t(floor((*x) * step));
    index[1] = size_t(floor((*y) * step));
    index[2] = size_t(floor((*z) * step));
    return index;
}
Sphere::Sphere(const Aggregate &Agg) noexcept : Sphere() {
    if (physicalmodel && physicalmodel->toBeDestroyed) {
        delete physicalmodel;
    }
    physicalmodel = Agg.physicalmodel;
    init_val(Agg.GetPosition(), *Agg.rmax);
}
void Aggregate::Merge(Aggregate &other) {
    //$ Update of the labels of the spheres that were in the deleted aggregate
    //$ And their new relative position
    double dx = periodicDistance(*other.myspheres[0].x - *myspheres[0].x, physicalmodel->L);
    double dy = periodicDistance(*other.myspheres[0].y - *myspheres[0].y, physicalmodel->L);
    double dz = periodicDistance(*other.myspheres[0].z - *myspheres[0].z, physicalmodel->L);

    // For all the spheres that were in the deleted aggregate
    array<double, 3> refpos = myspheres[0].get_position();
    for (Sphere *othersphere : other.myspheres) {
        // change the Label to the new owner
        othersphere->set_label(long(GetLabel()));

        // change the relative position to the new aggregate
        othersphere->relative_translate({dx, dy, dz});

        // Move them accordingly (periodicity)
        array<double, 3> newpos = othersphere->get_relative_position();
        newpos[0] += refpos[0];
        newpos[1] += refpos[1];
        newpos[2] += refpos[2];
        othersphere->set_position(newpos);
    }

    // Merge the spheresLists
    myspheres.merge(other.myspheres);
    Np = myspheres.size();
    UpdateDistances();
    Update();
    fullStatistics();
}
void Aggregate::DecreaseLabel() noexcept {
    if (InVerlet) {
        verlet->Remove(GetLabel(), IndexVerlet);
    }


    // Keep index and label in sync
    DecreaseIndex();
    Label--;

    // Keep aggLabel of myspheres in sync
    myspheres.DecreaseLabel();
    if (InVerlet) {
        verlet->Add(GetLabel(), IndexVerlet);
    }
}
void Aggregate::UpdateDistances() noexcept {
    _distances.resize(Np);
    distances_center.resize(Np);
    const size_t loopsize(Np);
    for (size_t i = 0; i < loopsize; i++) {
        // The last index is the distance to the mass center
        _distances[i].clear();
    }
    for (size_t i = 0; i < loopsize; i++) {
        for (size_t j = i + 1; j < loopsize; j++) {
            // Compute the distance between sphere i and j without taking periodicity into account
            double dist = relative_distance(myspheres[i], myspheres[j]);

            // if in contact
            if (dist < myspheres[i].get_radius() + myspheres[j].get_radius()) {
                _distances[i].push_back(pair<size_t, double>(j, dist));
                // distances are symetric !
                _distances[j].push_back(pair<size_t, double>(i, dist));
            }
        }
    }
}
/** Copy constructor */
Aggregate::Aggregate(const Aggregate &other) :
    storage_elem<15, ListAggregat>(other),
    StatisicsData(other),
    physicalmodel(other.physicalmodel),
    myspheres(other.myspheres),
    verlet(nullptr),
    IndexVerlet({{0, 0, 0}}),
    _distances(other._distances),
    distances_center(other.distances_center),
    volumes(other.volumes),
    surfaces(other.surfaces),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    time(nullptr),
    Np(other.Np),
    Label(other.Label),
    InVerlet(false) {
    setpointers();
}
/** Copy constructor */
Aggregate::Aggregate(const Aggregate &other, ListAggregat &_Storage) :
    storage_elem<15, ListAggregat>(other, *this, _Storage),
    StatisicsData(other),
    physicalmodel(other.physicalmodel),
    myspheres(other.myspheres, _Storage.spheres),
    verlet(nullptr),
    IndexVerlet({{0, 0, 0}}),
    _distances(other._distances),
    distances_center(other.distances_center),
    volumes(other.volumes),
    surfaces(other.surfaces),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    time(nullptr),
    Np(other.Np),
    Label(other.Label),
    InVerlet(false) {
    setpointers();
    for (Sphere *s : myspheres) {
        s->set_label(long(indexInStorage));
    }
}
/** Move constructor */
Aggregate::Aggregate(Aggregate &&other) noexcept: /* noexcept needed to enable optimizations in containers */
    storage_elem<15, ListAggregat>(move(other)),
    StatisicsData(move(other)),
    physicalmodel(other.physicalmodel),
    myspheres(move(other.myspheres)),
    verlet(nullptr),
    IndexVerlet({{0, 0, 0}}),
    _distances(),
    distances_center(),
    volumes(),
    surfaces(),
    rg(nullptr),
    f_agg(nullptr),
    lpm(nullptr),
    time_step(nullptr),
    rmax(nullptr),
    volAgregat(nullptr),
    surfAgregat(nullptr),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    time(nullptr),
    Np(other.Np),
    Label(other.Label),
    InVerlet(false) {
    swap(_distances, other._distances);
    swap(distances_center, other.distances_center);
    swap(volumes, other.volumes);
    swap(surfaces, other.surfaces);
    setpointers();
}


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
Aggregate &Aggregate::operator=(Aggregate &&other) noexcept {
    if (physicalmodel && physicalmodel->toBeDestroyed) {
        delete physicalmodel;
    }
    physicalmodel = other.physicalmodel;
    other.physicalmodel = nullptr;
    swap(myspheres, other.myspheres);
    swap(_distances, other._distances);
    swap(distances_center, other.distances_center);
    swap(volumes, other.volumes);
    swap(surfaces, other.surfaces);
    swap(IndexVerlet, other.IndexVerlet);
    swap(Label, other.Label);
    swap(Np, other.Np);
    StatisicsData::operator=(static_cast<StatisicsData &>(other));
    storage_elem<15, ListAggregat>::operator=(move(static_cast<storage_elem<15, ListAggregat> &>(other)));
    setpointers();
    return *this;
}
double Aggregate::SphereDistance(size_t i, size_t j) const {
    size_t ii = i;
    size_t jj = j;
    if (_distances[i].size() > _distances[j].size()) {
        ii = j;
        jj = i;
    }
    std::function<bool(pair<size_t, double>)> IsJJ = [jj](pair<size_t, double> e) {
        return e.first == jj;
    };
    auto findIter = std::find_if(_distances[ii].begin(), _distances[ii].end(), IsJJ);
    if (findIter != _distances[ii].end()) {
        return findIter->second;
    }
    return -1;
}
__attribute__((pure)) double Aggregate::SphereDistance(size_t i) const {
    return distances_center[i];
}
void Aggregate::print() const {
    cout << "Printing details of Aggregat " << indexInStorage << " " << Label << endl;
    if (external_storage) {
        cout << "  With external Storage" << endl;
    } else {
        cout << "  Without external Storage" << endl;
    }
    if (InVerlet) {
        cout << "  In Verlet list : " << IndexVerlet[0] << " "
             << IndexVerlet[1] << " "
             << IndexVerlet[2] << endl;
    } else {
        cout << "  Not in Verlet list" << endl;
    }
    cout << "  Caracteristics" << endl;
    cout << "    Gyration radius   : " << *rg << endl;
    cout << "    Geometric radius  : " << *rmax << endl;
    cout << "    Friction coeff.   : " << *f_agg << endl;
    cout << "    Mean Free Path    : " << *lpm << endl;
    cout << "    Delta t           : " << *time_step << endl;
    cout << "    Volume            : " << *volAgregat << endl;
    cout << "    Surface           : " << *surfAgregat << endl;
    cout << "    Position          : " << *x << " " << *y << " " << *z << endl;
    cout << "    Proper time       : " << *time << endl;
    myspheres.print();

    /*

    std::vector<std::vector <double > > distances;
    std::vector<double> volumes;
    std::vector<double> surfaces;
    double *rx,*ry,*rz; // position of the gravity center

    */
}
}  // namespace MCAC
