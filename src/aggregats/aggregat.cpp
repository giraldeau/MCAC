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
#include "spheres/sphere_distance.hpp"
#include "spheres/sphere_intersection.hpp"
#include "tools/tools.hpp"
#include "verlet/verlet.hpp"
#include "exceptions.hpp"
#include <iostream>


namespace mcac {
/* getters */
[[gnu::pure]] const double &Aggregate::get_rg() const noexcept {
    return *rg;
}
[[gnu::pure]] const double &Aggregate::get_f_agg() const noexcept {
    return *f_agg;
}
[[gnu::pure]] const double &Aggregate::get_lpm() const noexcept {
    return *lpm;
}
[[gnu::pure]] const double &Aggregate::get_time_step() const noexcept {
    return *time_step;
}
[[gnu::pure]] const double &Aggregate::get_rmax() const noexcept {
    return *rmax;
}
[[gnu::pure]] const double &Aggregate::get_agregat_volume() const noexcept {
    return *agregat_volume;
}
[[gnu::pure]] const double &Aggregate::get_agregat_surface() const noexcept {
    return *agregat_surface;
}
[[gnu::pure]] std::array<double, 3> Aggregate::get_position() const noexcept {
    std::array<double, 3> mypos{{*x, *y, *z}};
    return mypos;
}
[[gnu::pure]] std::array<double, 3> Aggregate::get_relative_position() const noexcept {
    std::array<double, 3> mypos{{*rx, *ry, *rz}};
    return mypos;
}
[[gnu::pure]] std::array<size_t, 3> Aggregate::get_verlet_index() const noexcept {
    return index_verlet;
}
[[gnu::pure]] const double &Aggregate::get_time() const noexcept {
    return *time;
}
[[gnu::pure]] const size_t &Aggregate::size() const noexcept {
    return n_spheres;
}
[[gnu::pure]] const size_t &Aggregate::get_label() const noexcept {
    return label;
}
/* modifiers */
void Aggregate::decrease_label() noexcept {
    if (static_cast<bool>(verlet)) {
        verlet->remove(get_label(), index_verlet);
    }
    label--;

    // Keep aggLabel of myspheres in sync
    myspheres.decrease_label();
    if (static_cast<bool>(verlet)) {
        verlet->add(get_label(), index_verlet);
    }
}
void Aggregate::set_verlet(Verlet *newverlet) noexcept {
    verlet = newverlet;
    index_verlet = compute_index_verlet();
    verlet->add(get_label(), index_verlet);
}
void Aggregate::unset_verlet() noexcept {
    if (static_cast<bool>(verlet)) {
        verlet->remove(label, index_verlet);
    }
    verlet = nullptr;
}
void Aggregate::set_time(double newtime) noexcept {
    *time = newtime;
}
void Aggregate::time_forward(double deltatemps) noexcept {
    *time += deltatemps;
}
void Aggregate::set_position(const std::array<double, 3> &position) noexcept {
    std::array<double, 3> newposition{periodic_position(position, physicalmodel->box_lenght)};
    *x = newposition[0];
    *y = newposition[1];
    *z = newposition[2];
    if (static_cast<bool>(verlet)) {
        //$ update Verlet
        update_verlet();
    }
}
void Aggregate::translate(std::array<double, 3> vector) noexcept {
    // move the aggregate
    set_position(get_position() + vector);

    // move the first sphere taking care of the periodicity
    std::array<double, 3> refpos{get_position() - get_relative_position()};
    myspheres[0].set_position(refpos);

    // move all the other sphere relatively to the first
    for (const auto& mysphere : myspheres) {
        std::array<double, 3> relpos = mysphere->get_relative_position();
        mysphere->set_position(refpos + relpos);
    }
}
void Aggregate::init(size_t new_label,
                     size_t sphere_index) {
    //initialize data

    //random size
    double diameter = 0;
    if (physicalmodel->monomeres_initialisation_type == MonomeresInitialisationMode::NORMAL_INITIALISATION) {
        diameter = (physicalmodel->mean_diameter)
                   + std::sqrt(2.) * physicalmodel->dispersion_diameter * inverf(2. * random() - 1.0);
    } else if (physicalmodel->monomeres_initialisation_type
               == MonomeresInitialisationMode::LOG_NORMAL_INITIALISATION) {
        diameter = std::pow(physicalmodel->mean_diameter,
                            std::sqrt(2.) * std::log(physicalmodel->dispersion_diameter) * inverf(2. * random() - 1.0));
    } else {
        throw InputError("Monomere initialisation mode unknown");
    }
    if (diameter <= 0) {
        diameter = physicalmodel->mean_diameter;
    }
    diameter *= 1E-9;
    for (size_t n_try = 0; n_try < external_storage->spheres.size(); n_try++) {
        //random position
        std::array<double, 3> newpos{{random() * physicalmodel->box_lenght,
                                      random() * physicalmodel->box_lenght,
                                      random() * physicalmodel->box_lenght}};
        if (external_storage->test_free_space(newpos, diameter * 0.5)) {
            init(*external_storage->physicalmodel,
                 &external_storage->spheres,
                 &external_storage->verlet,
                 external_storage->physicalmodel->time,
                 new_label, sphere_index, newpos, diameter);
            return;
        }
    }
    throw TooDenseError();
}
void Aggregate::init(const PhysicalModel &new_physicalmodel,
                     SphereList *spheres,
                     Verlet *new_verlet,
                     double new_time,
                     size_t new_label,
                     size_t sphere_index,
                     const std::array<double, 3> &position,
                     double sphere_diameter) noexcept {
    physicalmodel = &new_physicalmodel;
    label = new_label;
    index_in_storage = new_label;
    if (static_cast<bool>(external_storage)) {
        external_storage->setpointers();
    }
    setpointers();
    set_position(position);
    *time = new_time;
    if (static_cast<bool>(new_verlet)) {
        set_verlet(new_verlet);
    }
    (*spheres)[sphere_index].set_label(int(label));
    (*spheres)[sphere_index].init_val(position, sphere_diameter * 0.5);
    myspheres = SphereList(spheres, {sphere_index});
    n_spheres = myspheres.size();
    update_distances();
    update();
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
    double relax_time = mcac::PhysicalModel::relax_time(masse, *f_agg);
    *time_step = 3. * relax_time;
    double diffusivity = physicalmodel->diffusivity(*f_agg);
    *lpm = std::sqrt(6. * diffusivity * (*time_step));
    *dp = 0.;
    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < n_spheres; i++) {
        *dp += myspheres[i].get_radius();
    }
    *dp = 2 * (*dp) / double(n_spheres);
    *dg_over_dp = 2 * (*rg) / (*dp);
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
    std::pair<std::vector<double>, std::vector<double>> volume_surface(compute_volume_surface_sbl(myspheres));

    volumes = volume_surface.first;
    surfaces = volume_surface.second;

    for (size_t i = 0; i < n_spheres; i++)
    {
        *agregat_volume = *agregat_volume + volumes[i];    //Total compute_volume of Agg id
        *agregat_surface = *agregat_surface + surfaces[i]; //Total Surface of Agg id
    }
#else
    for (size_t i = 0; i < n_spheres; i++) {
#ifdef FULL_INTERNAL_DISTANCES
        for (size_t j = i + 1; j < n_spheres; j++) //for the j spheres composing Aggregate n°id
        {
            double dist = internal_sphere_distance(i, j);
#else
        for (const auto &[j, dist] : distances[i]) {
            if (j <= i) {
                continue;
            }
#endif
            //$ Calculation of the intersection between the spheres i and j if in contact
            Intersection intersection(myspheres[i],
                                      myspheres[j],
                                      dist);

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
    const size_t _loopsize{n_spheres};
    std::array<double, 3> r{0., 0., 0.};

    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < _loopsize; i++) {
        //$ Calculation of the position of the center of mass
        r += myspheres[i].get_relative_position() * volumes[i];
    }
    r /= *agregat_volume;
    for (size_t i = 0; i < _loopsize; i++) {
        std::array<double, 3> diff{myspheres[i].get_relative_position() - r};
        distances_center[i] = std::sqrt(std::pow(diff[0], 2) + std::pow(diff[1], 2) + std::pow(diff[2], 2));
    }
    set_position(myspheres[0].get_position() + r);
    *rx = r[0];
    *ry = r[1];
    *rz = r[2];
}
void Aggregate::compute_max_radius() noexcept {
    // Maximum radius of the aggregate,
    // this corresponds to the distance between
    //    the center of mass of the aggregate
    //    and the edge of the furthest ball from said center.
    // It is used to assimilate the aggregate to a sphere when checking for intersections
    *rmax = 0.0;
    const size_t _loopsize(n_spheres);
    for (size_t i = 0; i < _loopsize; i++) {
        *rmax = std::max(*rmax, myspheres[i].get_radius() + sphere_distance_center(i));
    }
}
void Aggregate::compute_giration_radius() noexcept {
    // This function determines the Gyration Radius of the Aggregate Id.

    // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient,
    // they are used  used in the final formula of the Radius of Gyration
    double arg(0.);
    double brg(0.);
    const size_t _loopsize(n_spheres);
    for (size_t i = 0; i < _loopsize; i++) {
        //$ Calculation of Rg
        arg = arg + volumes[i] * std::pow(sphere_distance_center(i), 2); // distance to the gravity center
        brg = brg + volumes[i] * std::pow(myspheres[i].get_radius(), 2);
    }
    *rg = std::sqrt(std::abs((arg + 3. / 5. * brg) / (*agregat_volume)));
    *agregat_volume = std::abs(*agregat_volume);
}
//#####################################################################################################################

void Aggregate::merge(Aggregate *other, AggregateContactInfo contact_info) noexcept {
    size_t i_mysphere, i_othersphere;
    if (contact_info.moving_aggregate == get_label()) {
        i_mysphere = contact_info.moving_sphere;
        i_othersphere = contact_info.other_sphere;
    } else {
        i_mysphere = contact_info.other_sphere;
        i_othersphere = contact_info.moving_sphere;
    }

    //$ update of the labels of the spheres that were in the deleted aggregate
    //$ And their new relative position
    std::array<double, 3> refpos = myspheres[0].get_position();

    //use periodic_distance at contact point
    std::array<double, 3> ref_root_to_contact = myspheres[i_mysphere].get_relative_position();
    std::array<double, 3> diffcontact =
        periodic_distance(other->myspheres[i_othersphere].get_position() - myspheres[i_mysphere].get_position(),
                          physicalmodel->box_lenght);
    std::array<double, 3> other_root_to_contact = other->myspheres[i_othersphere].get_relative_position();
    std::array<double, 3> diffpos = ref_root_to_contact + diffcontact - other_root_to_contact;

    // For all the spheres that were in the deleted aggregate
    for (const auto& othersphere : other->myspheres) {
        // change the Label to the new owner
        othersphere->set_label(long(get_label()));

        // change the relative position to the new aggregate
        othersphere->relative_translate(diffpos);

        // Move them accordingly (periodicity)
        std::array<double, 3> newpos = othersphere->get_relative_position();
        newpos += refpos;
        othersphere->set_position(newpos);
    }

    // merge the spheresLists
    myspheres.merge(other->myspheres);
    n_spheres = myspheres.size();
    update_distances();
    update();
}
bool Aggregate::split() {
    bool have_splitted = false;
    // first identify what to split
    std::vector<std::vector<size_t>> independant_components;
    // all spheres are unvisited
    std::vector<size_t> unvisisted;
    for (size_t i = 0; i < n_spheres; ++i) {
        unvisisted.insert(unvisisted.end(), i);
    }
    while (!unvisisted.empty()) {
        // the first unvisited sphere is the start of a new component
        std::vector<size_t> component({*unvisisted.begin()});
        unvisisted.erase(unvisisted.begin());
        // discovered are spheres of which we need to explore the neighborhood
        std::vector<size_t> discovered = component;
        while (!discovered.empty()) {
            size_t agg_1 = discovered.back();
            discovered.pop_back();
            auto suspect = unvisisted.begin();
            while (suspect != unvisisted.end()) {
                auto index = std::distance(unvisisted.begin(), suspect);
                if (contact(myspheres[agg_1], myspheres[*suspect])) {
                    // discovered are spheres of which we need to explore the neighborhood
                    component.push_back(*suspect);
                    discovered.push_back(*suspect);
                    unvisisted.erase(suspect);
                    suspect = unvisisted.begin() + index;
                } else {
                    suspect = unvisisted.begin() + index + 1;
                }
            }
        }
        independant_components.push_back(component);
    }
    if (independant_components.size() > 1) {
        // splitting occurs when we have at least 2 independant componnents
        have_splitted = true;
        auto initial_number_of_spheres = external_storage->spheres.size();
        for (auto &split: independant_components) {
            // duplicate the current aggregate
            auto *agg = external_storage->add(*this, *external_storage);
            agg->label = external_storage->size() - 1;
            external_storage->setpointers();
            // the verlet reference is not conserved by the duplication
            agg->set_verlet(verlet);
            // destroy the spheres of the duplication
            for (size_t i = 0; i < agg->n_spheres; i++) {
                external_storage->spheres.remove(initial_number_of_spheres);
            }
            // copy reference of the selection into the duplication
            agg->myspheres = SphereList(&myspheres, split);
            agg->n_spheres = agg->myspheres.size();
            for (const auto& sph : agg->myspheres) {
                sph->set_label(long(agg->get_label()));
            }
            // by creating and destroying spheres, this is important
            external_storage->spheres.setpointers();
            std::array<double, 3> refpos = agg->myspheres[0].get_position();
            for (const auto& sph : agg->myspheres) {
                // change the relative position of the new aggregate
                sph->set_relative_position(sph->get_position() - refpos);
            }
            // we have to recompute all the caracteristic of this new aggregate
            agg->update_distances();
            agg->update();
            // we have to move all the spheres (periodicity)
            refpos = agg->get_position() - agg->get_relative_position();
            for (const auto& sph : agg->myspheres) {
                sph->set_position(refpos + sph->get_relative_position());
            }
        }
    }
    return have_splitted;
}
void Aggregate::remove(const size_t &id) noexcept {
    for (size_t local_id = 0; local_id < n_spheres; local_id++) {
        if (myspheres[local_id].index_in_storage == id) {
            myspheres.list.erase(myspheres.list.begin() + long(local_id));
            break;
        }
    }
    external_storage->spheres.remove(id);
    n_spheres = myspheres.size();
    if (n_spheres > 0) {
        std::array<double, 3> refpos = myspheres[0].get_position();
        for (const auto& sph : myspheres) {
            // change the relative position of the new aggregate
            sph->set_relative_position(sph->get_position() - refpos);
        }
        // we have to recompute all the caracteristic of this new aggregate
        update_distances();
        update();
        // we have to move all the spheres (periodicity)
        refpos = get_position() - get_relative_position();
        for (const auto& sph : myspheres) {
            sph->set_position(refpos + sph->get_relative_position());
        }
    } // else it should be deleted ASAP
}
void Aggregate::print() const noexcept {
    std::cout << "Printing details of Aggregat " << index_in_storage << " " << label << std::endl;
    if (static_cast<bool>(external_storage)) {
        std::cout << "  With external Storage" << std::endl;
    } else {
        std::cout << "  Without external Storage" << std::endl;
    }
    if (static_cast<bool>(verlet)) {
        std::cout << "  In Verlet list : "
                  << index_verlet[0] << " "
                  << index_verlet[1] << " "
                  << index_verlet[2] << std::endl;
    } else {
        std::cout << "  Not in Verlet list" << std::endl;
    }
    std::cout << "  Caracteristics" << std::endl;
    std::cout << "    Gyration radius   : " << *rg << std::endl;
    std::cout << "    Geometric radius  : " << *rmax << std::endl;
    std::cout << "    Friction coeff.   : " << *f_agg << std::endl;
    std::cout << "    Mean Free Path    : " << *lpm << std::endl;
    std::cout << "    Delta t           : " << *time_step << std::endl;
    std::cout << "    Volume            : " << *agregat_volume << std::endl;
    std::cout << "    Surface           : " << *agregat_surface << std::endl;
    std::cout << "    Position          : " << *x << " " << *y << " " << *z << std::endl;
    std::cout << "    Proper time       : " << *time << std::endl;
    myspheres.print();
}
std::array<size_t, 3> Aggregate::compute_index_verlet() noexcept {
    double step = double(physicalmodel->n_verlet_divisions) / physicalmodel->box_lenght;
    std::array<size_t, 3> new_verlet_index{size_t(std::floor((*x) * step)),
                                           size_t(std::floor((*y) * step)),
                                           size_t(std::floor((*z) * step))};
    return new_verlet_index;
}
void Aggregate::update_verlet() noexcept {
    std::array<size_t, 3> new_verlet_index = compute_index_verlet();
    auto[a, b, c] = new_verlet_index;
    auto[i, j, k] = index_verlet;
    if (a != i ||
        b != j ||
        c != k) {
        verlet->remove(get_label(), index_verlet);
        index_verlet = new_verlet_index;
        verlet->add(get_label(), new_verlet_index);
    }
}
void Aggregate::update_distances() noexcept {
    distances.resize(n_spheres);
    distances_center.resize(n_spheres);
    const size_t _loopsize(n_spheres);
    for (size_t i = 0; i < _loopsize; i++) {
        // The last index is the distance to the mass center
#ifdef FULL_INTERNAL_DISTANCES
        distances[i].resize(n_spheres);
#else
        size_t old_size = distances[i].size();
        distances[i].clear();
        distances[i].reserve(old_size);
#endif
    }
    for (size_t i = 0; i < _loopsize; i++) {
        for (size_t j = i + 1; j < _loopsize; j++) {
            // Compute the distance between sphere i and j without taking periodicity into account
            double dist = relative_distance(myspheres[i], myspheres[j]);
#ifdef FULL_INTERNAL_DISTANCES
            distances[i][j] = dist;
            // distances are symetric !
            distances[j][i] = dist;
#else
            // if in contact (with some tolerence)
            if (dist <= (1. + 1e-10) * (myspheres[i].get_radius() + myspheres[j].get_radius())) {
                distances[i][j] = dist;
                // distances are symetric !
                distances[j][i] = dist;
            }
#endif
        }
    }
}
double Aggregate::internal_sphere_distance(size_t i, size_t j) const noexcept {
#ifdef FULL_INTERNAL_DISTANCES
    return distances[i][j];
#else
    size_t ii = i;
    size_t jj = j;
    if (distances[i].size() > distances[j].size()) {
        ii = j;
        jj = i;
    }
    auto it = distances[ii].find(jj);
    if (it != distances[ii].end()) {
        return it->second;
    }
    return std::numeric_limits<double>::infinity();
#endif
}
[[gnu::pure]] double Aggregate::sphere_distance_center(size_t i) const noexcept {
    return distances_center[i];
}
}  // namespace mcac
