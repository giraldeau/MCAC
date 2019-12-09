#include "aggregats/aggregat.hpp"
#include "aggregats/aggregat_list.hpp"
#include "spheres/sphere.hpp"
#include "spheres/sphere_intersection.hpp"
#include "spheres/sphere_distance.hpp"
#include "verlet/verlet.hpp"
#include "tools/tools.hpp"
#include <iostream>


using namespace std;
namespace MCAC {
/* getters */
[[gnu::pure]] double Aggregate::get_rg() const noexcept {
    return *rg;
}
[[gnu::pure]] double Aggregate::get_f_agg() const noexcept {
    return *f_agg;
}
[[gnu::pure]] double Aggregate::get_lpm() const noexcept {
    return *lpm;
}
[[gnu::pure]] double Aggregate::get_time_step() const noexcept {
    return *time_step;
}
[[gnu::pure]] double Aggregate::get_rmax() const noexcept {
    return *rmax;
}
[[gnu::pure]] double Aggregate::get_agregat_volume() const noexcept {
    return *agregat_volume;
}
[[gnu::pure]] double Aggregate::get_agregat_surface() const noexcept {
    return *agregat_surface;
}
[[gnu::pure]] array<double, 3> Aggregate::get_position() const noexcept {
    array<double, 3> mypos{{*x, *y, *z}};
    return mypos;
}
[[gnu::pure]] array<double, 3> Aggregate::get_relative_position() const noexcept {
    array<double, 3> mypos{{*rx, *ry, *rz}};
    return mypos;
}
[[gnu::pure]] std::array<size_t, 3> Aggregate::get_verlet_index() const noexcept {
    return index_verlet;
}
[[gnu::pure]] double Aggregate::get_time() const noexcept {
    return *time;
}
[[gnu::pure]] size_t Aggregate::size() const noexcept {
    return n_spheres;
}
[[gnu::pure]] size_t Aggregate::get_label() const noexcept {
    return label;
}
/* modifiers */
void Aggregate::decrease_label() noexcept {
    if (static_cast<bool>(verlet)) {
        verlet->Remove(get_label(), index_verlet);
    }
    // Keep index and label in sync
    decrease_index();
    label--;

    // Keep aggLabel of myspheres in sync
    myspheres.decrease_label();
    if (static_cast<bool>(verlet)) {
        verlet->Add(get_label(), index_verlet);
    }
}
void Aggregate::set_verlet(Verlet &newverlet) noexcept {
    verlet = &newverlet;
}
void Aggregate::unset_verlet() noexcept {
    verlet = nullptr;
}
void Aggregate::set_time(double newtime) noexcept {
    *time = newtime;
}
void Aggregate::time_forward(double deltatemps) noexcept {
    *time += deltatemps;
}
void Aggregate::set_position(const array<double, 3> &position) noexcept {
    array<double, 3> newposition{periodic_position(position, physicalmodel->box_lenght)};
    *x = newposition[0];
    *y = newposition[1];
    *z = newposition[2];
    if (static_cast<bool>(verlet)) {
        //$ update Verlet
        update_verlet_index();
    }
}
void Aggregate::translate(array<double, 3> vector) noexcept {
    // move the aggregate
    set_position(get_position() + vector);

    // move the first sphere taking care of the periodicity
    array<double, 3> refpos{periodic_position(myspheres[0].get_position() + vector, physicalmodel->box_lenght)};
    myspheres[0].set_position(refpos);

    // move all the other sphere relatively to the first
    for (Sphere *mysphere : myspheres) {
        array<double, 3> relpos = mysphere->get_relative_position();
        mysphere->set_position(refpos + relpos);
    }
}
void Aggregate::init(const PhysicalModel &new_physicalmodel,
                     Verlet &new_verlet,
                     const array<double, 3> &position,
                     size_t new_label,
                     SphereList &spheres,
                     double sphere_diameter) noexcept {
    physicalmodel = &new_physicalmodel;
    if (static_cast<bool>(verlet)) {
        verlet->Remove(label, index_verlet);
    }
    label = new_label;
    index_in_storage = new_label;
    if (static_cast<bool>(external_storage)) {
        external_storage->setpointers();
    }
    setpointers();
    set_position(position);
    verlet = &new_verlet;
    update_verlet_index();
    spheres[label].set_label(int(label));
    spheres[label].init_val(position, sphere_diameter / 2);
    myspheres = SphereList(spheres, {label});
    n_spheres = myspheres.size();
    update_distances();
    update();
    fullStatistics();
}

//###################################################################################################################

//### Mise à jour des paramètres physiques d'un agrégat (rayon de giration, masse, nombre de sphérules primaires) #####
void Aggregate::update() noexcept {
    // This function will update the parameter of Agg
    compute_volume();
    compute_mass_center();
    compute_max_radius();
    compute_giration_radius();

    //$ Determination of the friction coefficient
    *f_agg = physicalmodel->friction_coeff(n_spheres);
    double masse = physicalmodel->density * (*agregat_volume);
    double relax_time = physicalmodel->relax_time(masse, *f_agg);
    *time_step = 3.0 * relax_time;
    double diffusivity = physicalmodel->diffusivity(*f_agg);
    *lpm = sqrt(6.0 * diffusivity * (*time_step));
    partialStatistics();
    if (static_cast<bool>(external_storage)) {
        if (*rmax > external_storage->maxradius) {
            external_storage->maxradius = *rmax;
        }
    }
}
//####### Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ########
void Aggregate::compute_volume() noexcept {
    *agregat_volume = *agregat_surface = 0.0; // compute_volume and surface of Agg Id

    //$ Initialisation of the arrays of volume, surface of each sphere, and the center of mass
    volumes.resize(n_spheres);
    surfaces.resize(n_spheres);

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < n_spheres; i++) {
        //$ Calculation of the volume and surface of monomere i of Agg id
        volumes[i] = myspheres[i].get_volume();      //Calculation of the volume of i
        surfaces[i] = myspheres[i].get_surface();    //Calculation of the surface of i
    }
#ifdef WITH_SBL
    pair<vector<double>, vector<double>> volume_surface(compute_volume_surface(myspheres));

    volumes = volume_surface.first;
    surfaces = volume_surface.second;

    for (size_t i = 0; i < n_spheres; i++)
    {
        *agregat_volume = *agregat_volume + volumes[i];    //Total compute_volume of Agg id
        *agregat_surface = *agregat_surface + surfaces[i]; //Total Surface of Agg id
    }
#else
    for (size_t i = 0; i < n_spheres; i++) {
        for (size_t j = i + 1; j < n_spheres; j++) //for the j spheres composing Aggregate n°id
        {
            //$ Calculation of the intersection between the spheres i and j if in contact
            Intersection intersection(myspheres[i],
                                      myspheres[j],
                                      internal_sphere_distance(i, j));

            //$ The volume and surface covered by j is substracted from those of i
            volumes[i] = volumes[i] - intersection.volume_1;
            surfaces[i] = surfaces[i] - intersection.surface_1;

            //$ The volume and surface covered by i is substracted from those of j
            volumes[j] = volumes[j] - intersection.volume_2;
            surfaces[j] = surfaces[j] - intersection.surface_2;
        }
        //$ Calculation of the total volume and surface of the aggregate
        *agregat_volume = *agregat_volume + volumes[i];
        *agregat_surface = *agregat_surface + surfaces[i];
    }
#endif
}
void Aggregate::compute_mass_center() noexcept {
    const size_t loopsize{n_spheres};
    array<double, 3> _r{0., 0., 0.};

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < loopsize; i++) {
        //$ Calculation of the position of the center of mass
        _r += myspheres[i].get_relative_position() * volumes[i];
    }
    _r /= *agregat_volume;
    for (size_t i = 0; i < loopsize; i++) {
        array<double, 3> diff{myspheres[i].get_relative_position() - _r};
        distances_center[i] = sqrt(POW_2(diff[0]) + POW_2(diff[1]) + POW_2(diff[2]));
    }
    set_position(myspheres[0].get_position() + _r);
    *rx = _r[0];
    *ry = _r[1];
    *rz = _r[2];
}
void Aggregate::compute_max_radius() noexcept {
    // Maximum radius of the aggregate,
    // this corresponds to the distance between
    //    the center of mass of the aggregate
    //    and the edge of the furthest ball from said center.
    // It is used to assimilate the aggregate to a sphere when checking for intersections
    *rmax = 0.0;
    const size_t loopsize(n_spheres);
    for (size_t i = 0; i < loopsize; i++) {
        *rmax = MAX(*rmax, myspheres[i].get_radius() + sphere_distance_center(i));
    }
}
void Aggregate::compute_giration_radius() noexcept {
    // This function determines the Gyration Radius of the Aggregate Id.

    // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient,
    // they are used  used in the final formula of the Radius of Gyration
    double Arg(0.);
    double Brg(0.);
    const size_t loopsize(n_spheres);
    for (size_t i = 0; i < loopsize; i++) {
        //$ Calculation of Rg
        Arg = Arg + volumes[i] * POW_2(sphere_distance_center(i)); // distance to the gravity center
        Brg = Brg + volumes[i] * POW_2(myspheres[i].get_radius());
    }
    *rg = sqrt(fabs((Arg + 3.0 / 5.0 * Brg) / (*agregat_volume)));
    *agregat_volume = fabs(*agregat_volume);
}
//#####################################################################################################################

void Aggregate::merge(Aggregate &other) noexcept {
    //$ update of the labels of the spheres that were in the deleted aggregate
    //$ And their new relative position
    array<double, 3> refpos = myspheres[0].get_position();
    array<double, 3> diffpos = periodic_distance(other.myspheres[0].get_position() - refpos, physicalmodel->box_lenght);

    // For all the spheres that were in the deleted aggregate
    for (Sphere *othersphere : other.myspheres) {
        // change the Label to the new owner
        othersphere->set_label(long(get_label()));

        // change the relative position to the new aggregate
        othersphere->relative_translate(diffpos);

        // Move them accordingly (periodicity)
        array<double, 3> newpos = othersphere->get_relative_position();
        newpos += refpos;
        othersphere->set_position(newpos);
    }

    // merge the spheresLists
    myspheres.merge(other.myspheres);
    n_spheres = myspheres.size();
    update_distances();
    update();
    fullStatistics();
}
void Aggregate::print() const noexcept {
    cout << "Printing details of Aggregat " << index_in_storage << " " << label << endl;
    if (static_cast<bool>(external_storage)) {
        cout << "  With external Storage" << endl;
    } else {
        cout << "  Without external Storage" << endl;
    }
    if (static_cast<bool>(verlet)) {
        cout << "  In Verlet list : "
             << index_verlet[0] << " "
             << index_verlet[1] << " "
             << index_verlet[2] << endl;
    } else {
        cout << "  Not in Verlet list" << endl;
    }
    cout << "  Caracteristics" << endl;
    cout << "    Gyration radius   : " << *rg << endl;
    cout << "    Geometric radius  : " << *rmax << endl;
    cout << "    Friction coeff.   : " << *f_agg << endl;
    cout << "    Mean Free Path    : " << *lpm << endl;
    cout << "    Delta t           : " << *time_step << endl;
    cout << "    Volume            : " << *agregat_volume << endl;
    cout << "    Surface           : " << *agregat_surface << endl;
    cout << "    Position          : " << *x << " " << *y << " " << *z << endl;
    cout << "    Proper time       : " << *time << endl;
    myspheres.print();
}
void Aggregate::update_verlet_index() noexcept {
    double step = double(physicalmodel->n_verlet_divisions) / physicalmodel->box_lenght;
    std::array<size_t, 3> new_verlet_index{size_t(floor((*x) * step)),
                                           size_t(floor((*y) * step)),
                                           size_t(floor((*z) * step))};
    verlet->Remove(get_label(), index_verlet);
    index_verlet = new_verlet_index;
    verlet->Add(get_label(), index_verlet);
}
void Aggregate::update_distances() noexcept {
    distances.resize(n_spheres);
    distances_center.resize(n_spheres);
    const size_t loopsize(n_spheres);
    for (size_t i = 0; i < loopsize; i++) {
        // The last index is the distance to the mass center
        distances[i].clear();
    }
    for (size_t i = 0; i < loopsize; i++) {
        for (size_t j = i + 1; j < loopsize; j++) {
            // Compute the distance between sphere i and j without taking periodicity into account
            double dist = relative_distance(myspheres[i], myspheres[j]);

            // if in contact
            if (dist < myspheres[i].get_radius() + myspheres[j].get_radius()) {
                distances[i].push_back(pair<size_t, double>(j, dist));
                // distances are symetric !
                distances[j].push_back(pair<size_t, double>(i, dist));
            }
        }
    }
}
double Aggregate::internal_sphere_distance(size_t i, size_t j) const noexcept {
    size_t ii = i;
    size_t jj = j;
    if (distances[i].size() > distances[j].size()) {
        ii = j;
        jj = i;
    }
    std::function<bool(pair<size_t, double>)> IsJJ = [jj](pair<size_t, double> e) {
        return e.first == jj;
    };
    auto findIter = std::find_if(distances[ii].begin(), distances[ii].end(), IsJJ);
    if (findIter != distances[ii].end()) {
        return findIter->second;
    }
    return -1;
}
__attribute__((pure)) double Aggregate::sphere_distance_center(size_t i) const noexcept {
    return distances_center[i];
}
}  // namespace MCAC
