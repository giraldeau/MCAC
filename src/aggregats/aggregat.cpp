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
#include "arvo/arvo_mcac_call.hpp"
#include <iostream>
#include <string.h>

namespace mcac {
/* getters */
[[gnu::pure]] const double &Aggregate::get_rg() const noexcept {
    return *rg;
}
[[gnu::pure]] const double &Aggregate::get_f_agg() const noexcept {
    return *f_agg;
}
[[gnu::pure]] const double &Aggregate::get_dp() const noexcept {
    return *dp;
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
[[gnu::pure]] const double &Aggregate::get_proper_time() const noexcept {
    return *proper_time;
}
[[gnu::pure]] const size_t &Aggregate::size() const noexcept {
    return n_spheres;
}
[[gnu::pure]] const size_t &Aggregate::get_label() const noexcept {
    return label;
}
[[gnu::pure]] const int &Aggregate::get_electric_charge() const noexcept {
    return electric_charge;
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
void Aggregate::set_proper_time(double newtime) noexcept {
    *proper_time = newtime;
}
void Aggregate::time_forward(double deltatemps) noexcept {
    *proper_time += deltatemps;
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
void Aggregate::set_bulk_density() {
    /*
    double dpp_nm = (*dp)*1e+09;
    if (physicalmodel->with_maturity){
        newdensity = (-4.30e-04)*std::pow(dpp_nm,4) +
                (9.12e-02)*std::pow(dpp_nm,3) +
                (-3.01e+00)*std::pow(dpp_nm,2) +
                (4.33e+01)*(dpp_nm) + 1.37e+03;
    } else {
        newdensity = physicalmodel->density;
    }
    */
    if (physicalmodel->with_maturity){
        set_CH_ratio();
        bulk_density = _bulk_density_young+
                (_bulk_density_mature-_bulk_density_young)/(_CH_mature-_CH_young) *
                ((*CH_ratio)-_CH_young);
        if (bulk_density < _bulk_density_young || bulk_density > _bulk_density_mature){
            throw InputError("Problem with bulk density: " + std::to_string(bulk_density));
        }
    } else {
        bulk_density = physicalmodel->density;
    }
}
void Aggregate::set_CH_ratio() noexcept {
    const double a(4.0), b(1.0);
    double dpp_nm = (*dp)*1e+09;
    *CH_ratio = 0.5*(std::erf((dpp_nm-a)/b)+1.0)*(_CH_mature-_CH_young)+_CH_young;
}
void Aggregate::translate(std::array<double, 3> vector) noexcept {
    // move the aggregate
    set_position(get_position() + vector);

    // move the first sphere taking care of the periodicity
    std::array<double, 3> refpos{get_position() - get_relative_position()};
    myspheres[0]->set_position(refpos);

    // move all the other sphere relatively to the first
    for (const auto &mysphere : myspheres) {
        std::array<double, 3> relpos = mysphere->get_relative_position();
        mysphere->set_position(refpos + relpos);
    }
}
void Aggregate::init(size_t new_label,
                     size_t sphere_index,
                     bool nucleation) {
    //initialize data
    set_proper_time(physicalmodel->time);

    //random size
    double diameter;
    if (nucleation) {
        diameter = physicalmodel->random_diameter(physicalmodel->mean_diameter_nucleation,
                                                  physicalmodel->dispersion_diameter_nucleation);
    } else {
        diameter = physicalmodel->random_diameter(physicalmodel->mean_diameter,
                                                  physicalmodel->dispersion_diameter);
    }

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
    *proper_time = new_time;
    if (static_cast<bool>(new_verlet)) {
        set_verlet(new_verlet);
    }
    (*spheres)[sphere_index]->set_label(int(label));
    (*spheres)[sphere_index]->init_val(position, sphere_diameter * 0.5);
    myspheres = SphereList(spheres, {sphere_index});
    n_spheres = myspheres.size();
    alpha_vs_extreme = 1.0/static_cast<double>(n_spheres);
    *d_m = sphere_diameter;
    // random initial charge
    if (physicalmodel->with_electric_charges) {
        electric_charge = physicalmodel->get_random_charge(*d_m);
    } else {
        electric_charge = 0;
    }
    (*spheres)[sphere_index]->set_sphere_charge(electric_charge);
    //update_distances_and_overlapping();
    update();
}
bool Aggregate::croissance_surface(double dt) {
    myspheres.croissance_surface(dt);

    for (auto sphere = myspheres.begin(); sphere != myspheres.end(); ) {
        if ((*sphere)->get_radius() <= physicalmodel->rp_min_oxid) {
            auto index = std::distance(myspheres.begin(), sphere);
            remove_sphere((*sphere)->get_index());
            sphere = myspheres.begin() + index;
        } else {
            sphere++ ;
        }
    }
    return myspheres.size() == 0;
}

//###################################################################################################################
//### Update all physical params of an aggregate except volume and surface (rayon de giration, masse, nombre de sphérules primaires) #####
void Aggregate::update_partial() noexcept {
    // This function will update the parameter of Agg
    compute_mass_center();
    compute_max_radius();
    compute_giration_radius();

    // Mean monomer diameter
    *dp = 0.;
    double vol_pp(0.0);
    //$ For the Spheres i in Agg Id
    for (size_t i = 0; i < n_spheres; i++) {
        vol_pp += myspheres[i]->get_volume();
        *dp += myspheres[i]->get_radius();
    }
    *dp = 2 * (*dp) / static_cast<double>(n_spheres);
    vol_pp = vol_pp / static_cast<double>(n_spheres);

    set_bulk_density();

    //$ Determination of the friction coefficient
    *f_agg = physicalmodel->friction_coeff(*agregat_volume, vol_pp, 0.5 * (*dp));
    *d_m = physicalmodel->mobility_diameter(*agregat_volume, vol_pp, 0.5 * (*dp));

    // Momentum relaxation time and lpm (persistent distance)
    double masse = bulk_density * (*agregat_volume);
    double relax_time = mcac::PhysicalModel::relax_time(masse, *f_agg);
    *time_step = 3. * relax_time;
    double diffusivity = physicalmodel->diffusivity(*f_agg);
    *lpm = sqrt(6. * diffusivity * (*time_step));
    *dg_over_dp = 2 * (*rg) / (*dp);
    if (static_cast<bool>(external_storage)) {
        if (*rmax > external_storage->maxradius) {
            external_storage->maxradius = *rmax;
        }
    }
}
//### Update all physical params of an aggregate (rayon de giration, masse, nombre de sphérules primaires) #####
void Aggregate::update() noexcept {
    update_distances_and_overlapping();
    compute_volume_surface();
    update_partial();
}
static double volume_alpha_correction(const double coordination_number,
                                      const double c_20,
                                      const double c_30,
                                      const double min_coordination_number,
                                      const double alpha_vs_extreme) noexcept {
    double difference_coordination = std::abs(coordination_number - min_coordination_number);
    double correction = 0.25 * (3.0 * c_20 - c_30) * coordination_number
                        - c_30 * difference_coordination * 0.62741833
                        - pow(difference_coordination, 1.5) * 0.00332425;
    // NOTE: correction parameters are obtained by post-processing fit
    if (correction < 0.0) {correction = 1.0;};
    correction = std::min(correction, 1.0);
    double alpha_v = 1.0 - correction;
    alpha_v = std::max(alpha_v, alpha_vs_extreme);
    return alpha_v;
}
static double surface_alpha_correction(const double coordination_number,
                                       const double c_10,
                                       const double min_coordination_number,
                                       const double alpha_vs_extreme) noexcept {
    double difference_coordination = std::abs(coordination_number - min_coordination_number);
    double correction = 0.5 * c_10 * coordination_number
                        - std::pow(c_10, 2) * difference_coordination * 0.29611
                        - std::pow(difference_coordination, 2) * 0.00155632;
    // NOTE: correction parameters are obtained by post-processing fit
    if (correction < 0.0) {correction = 1.0;};
    correction = std::min(correction, 1.0);
    double alpha_s = 1.0 - correction;
    alpha_s = std::max(alpha_s, alpha_vs_extreme);
    return alpha_s;
}
//####### Calculation of the volume, surface, center of mass and Giration radius of gyration of an aggregate ########
void Aggregate::compute_volume_surface() {
    *agregat_volume = *agregat_surface = 0.0; // compute_volume and surface of Agg Id

    //$ Initialisation of the arrays of volume, surface of each sphere, and the center of mass
    volumes.resize(n_spheres);
    surfaces.resize(n_spheres);

    // Compute volumes and surfaces using a library or ignoring overlap
    if (physicalmodel->volsurf_method == VolSurfMethods::EXACT_SBL) {
#ifdef WITH_SBL
        auto volume_surface(compute_volume_surface_sbl(myspheres));
        volumes = volume_surface.first;
        surfaces = volume_surface.second;
#else
        std::cout << "  ERROR SBL not available - activate it when doing cmake (cmake -DWITH_SBL=ON ..)" << std::endl;
        exit(ErrorCodes::SBL_ERROR);
#endif
    } else if (physicalmodel->volsurf_method == VolSurfMethods::EXACT_ARVO) {
#ifdef WITH_ARVO
        auto volume_surface(arvo_call(myspheres, distances));
        volumes = volume_surface.first;
        surfaces = volume_surface.second;
#else
        std::cout << "  ERROR ARVO not available - activate it when doing cmake (cmake -DWITH_ARVO=ON ..)" << std::endl;
        exit(ErrorCodes::ARVO_ERROR);
#endif

    } else {
        //$ For the Spheres i in Agg Id
        for (size_t i = 0; i < n_spheres; i++) {
            //$ Calculation of the volume and surface of monomere i of Agg id
            volumes[i] = myspheres[i]->get_volume();      //Calculation of the volume of i
            surfaces[i] = myspheres[i]->get_surface();    //Calculation of the surface of i
        }
    }

    // Use alpha or spherical caps to correct for the overlap
    if (physicalmodel->volsurf_method == VolSurfMethods::SPHERICAL_CAPS) {
        for (size_t i = 0; i < n_spheres; i++) {
#ifdef FULL_INTERNAL_DISTANCES
            for (size_t j = i + 1; j < n_spheres; j++) { //for the j spheres composing Aggregate n°id
                double dist = internal_sphere_distance(i, j);
#else
            for (const auto &[j, dist] : distances[i]) {
                if (j <= i) {
                    continue;
                }
#endif
                //$ Calculation of the intersection between the spheres i and j if in contact
                Intersection intersection(*myspheres[i],
                                          *myspheres[j],
                                          dist);

                //$ The volume and surface covered by j is substracted from those of i
                volumes[i] = volumes[i] - intersection.volume_1;
                surfaces[i] = surfaces[i] - intersection.surface_1;

                //$ The volume and surface covered by i is substracted from those of j
                volumes[j] = volumes[j] - intersection.volume_2;
                surfaces[j] = surfaces[j] - intersection.surface_2;
            }
            //$ For large overlapping, this function gives negative values, so in this case=0
            volumes[i] = std::max(volumes[i], 0.0);
            surfaces[i] = std::max(surfaces[i], 0.0);
        }
    }
    for (size_t i = 0; i < n_spheres; i++) {
        *agregat_volume = *agregat_volume + volumes[i];    //Total compute_volume of Agg id
        *agregat_surface = *agregat_surface + surfaces[i]; //Total Surface of Agg id
    }
    if (physicalmodel->volsurf_method == VolSurfMethods::ALPHAS) {
        *overlapping = 0.0; // average overlapping and coord. numbers updated here
        double c_v30(0.0), c_v20(0.0), c_s10(0.0);
        double vp_sum(0.0), sp_sum(0.0);
        size_t intersections(0);
        for (size_t i = 0; i < n_spheres; i++) {
            // loop over the neighbor list
            for (const auto &it : distances[i]) {
                double radius_1 = myspheres[i]->get_radius();
                double radius_2 = myspheres[it.first]->get_radius();
                double c_ij = (radius_1 + radius_2 - it.second) / (radius_1 + radius_2);
                double vp1 = std::pow(radius_1, 3);
                double vp2 = std::pow(radius_2, 3);
                double sp1 = std::pow(radius_1, 2);
                double sp2 = std::pow(radius_2, 2);
                vp_sum += (vp1 + vp2);
                sp_sum += (sp1 + sp2);
                *overlapping += c_ij;
                c_s10 += c_ij * (sp1 + sp2);
                c_v20 += std::pow(c_ij, 2) * (vp1 + vp2);
                c_v30 += std::pow(c_ij, 3) * (vp1 + vp2);
            }
            intersections += distances[i].size();
        }
        if (intersections > 0) {
            c_s10 /= sp_sum;
            c_v20 /= vp_sum;
            c_v30 /= vp_sum;
            double min_coordination_number = 2 * (1.0 - 1.0 / static_cast<double>(n_spheres));
            *overlapping /= static_cast<double>(intersections);
            *coordination_number = static_cast<double>(intersections) / static_cast<double>(n_spheres);

            *agregat_volume *= volume_alpha_correction(*coordination_number, c_v20, c_v30, min_coordination_number, alpha_vs_extreme);
            *agregat_surface *= surface_alpha_correction(*coordination_number, c_s10, min_coordination_number, alpha_vs_extreme);
        }
    }
    if (*agregat_volume <= 0 || *agregat_surface <= 0) {
        throw VolSurfError();
    }
}
void Aggregate::compute_mass_center() noexcept {
    const size_t _loopsize{n_spheres};
    std::array<double, 3> r{0., 0., 0.};

    //$ For the Spheres i in Agg Id
    // double sum_volume(0.0);
    for (size_t i = 0; i < _loopsize; i++) {
        //$ Calculation of the position of the center of mass
        r += myspheres[i]->get_relative_position() * volumes[i];
        // sum_volume += volumes[i];
    }
    r /= *agregat_volume;
    // r /= sum_volume;
    for (size_t i = 0; i < _loopsize; i++) {
        std::array<double, 3> diff{myspheres[i]->get_relative_position() - r};
        distances_center[i] = std::sqrt(std::pow(diff[0], 2) + std::pow(diff[1], 2) + std::pow(diff[2], 2));
    }
    set_position(myspheres[0]->get_position() + r);
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
        *rmax = std::max(*rmax, myspheres[i]->get_radius() + sphere_distance_center(i));
    }
}
void Aggregate::compute_giration_radius() noexcept {
    // This function determines the Gyration Radius of the Aggregate Id.

    // These correspond to the sum of the volumes of each spheres multiplied by their respective coefficient,
    // they are used  used in the final formula of the Radius of Gyration
    double arg(0.);
    double brg(0.);
    // double sum_volume(0.0);
    const size_t _loopsize(n_spheres);
    for (size_t i = 0; i < _loopsize; i++) {
        //$ Calculation of Rg
        arg = arg + volumes[i] * std::pow(sphere_distance_center(i), 2); // distance to the gravity center
        brg = brg + volumes[i] * std::pow(myspheres[i]->get_radius(), 2);
        // sum_volume += volumes[i];
    }
    // *rg = std::sqrt(std::abs((arg + 3. / 5. * brg) / sum_volume));
    *rg = std::sqrt(std::abs((arg + 3. / 5. * brg) / (*agregat_volume)));
    //*agregat_volume = std::abs(*agregat_volume);
}
//#####################################################################################################################

bool Aggregate::merge(std::shared_ptr<Aggregate> other, AggregateContactInfo contact_info) noexcept {

    std::shared_ptr<Aggregate> moving_aggregate = contact_info.moving_aggregate.lock();
    std::shared_ptr<Aggregate> other_aggregate = contact_info.other_aggregate.lock();
    if (!moving_aggregate || !other_aggregate) {
        // we lost one of the aggregates to be merged
        return false;
    }

    std::shared_ptr<Sphere> mysphere, othersphere;
    if (moving_aggregate.get() == this) {
        mysphere = contact_info.moving_sphere.lock();
        othersphere = contact_info.other_sphere.lock();
    } else if (other_aggregate.get() == this) {
        mysphere = contact_info.other_sphere.lock();
        othersphere = contact_info.moving_sphere.lock();
    } else {
        // ?
        return false;
    }
    if (!mysphere || !othersphere) {
        // we lost one of the spheres of the contact point
        return false;
    }

    //$ update of the labels of the spheres that were in the deleted aggregate
    //$ And their new relative position
    std::array<double, 3> refpos = myspheres[0]->get_position();

    //use periodic_distance at contact point
    std::array<double, 3> ref_root_to_contact = mysphere->get_relative_position();
    std::array<double, 3> diffcontact =
        periodic_distance(othersphere->get_position() - mysphere->get_position(),
                          physicalmodel->box_lenght);
    std::array<double, 3> other_root_to_contact = othersphere->get_relative_position();
    std::array<double, 3> diffpos = ref_root_to_contact + diffcontact - other_root_to_contact;

    // For all the spheres that were in the deleted aggregate
    for (const auto &sphere : other->myspheres) {
        // change the Label to the new owner
        sphere->set_label(long(get_label()));

        // change the relative position to the new aggregate
        sphere->relative_translate(diffpos);

        // Move them accordingly (periodicity)
        std::array<double, 3> newpos = sphere->get_relative_position();
        newpos += refpos;
        sphere->set_position(newpos);
    }

    // merge the spheresLists
    myspheres.merge(other->myspheres);
    n_spheres = myspheres.size();
    alpha_vs_extreme = 1.0/static_cast<double>(n_spheres);
    //update_distances_and_overlapping();
    update();
    return true;
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
                if (contact(*myspheres[agg_1], *myspheres[*suspect])) {
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
            auto agg = external_storage->add(*this, *external_storage);
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
            for (const auto &sph : agg->myspheres) {
                sph->set_label(long(agg->get_label()));
            }
            // by creating and destroying spheres, this is important
            external_storage->spheres.setpointers();
            std::array<double, 3> refpos = agg->myspheres[0]->get_position();
            for (const auto &sph : agg->myspheres) {
                // change the relative position of the new aggregate
                sph->set_relative_position(sph->get_position() - refpos);
            }
            // we have to recompute all the caracteristic of this new aggregate
            //agg->update_distances_and_overlapping();
            agg->update();
            // we have to move all the spheres (periodicity)
            refpos = agg->get_position() - agg->get_relative_position();
            for (const auto &sph : agg->myspheres) {
                sph->set_position(refpos + sph->get_relative_position());
            }
            // electric charge
            if (physicalmodel->with_dynamic_random_charges) {
                agg->electric_charge = physicalmodel->get_random_charge(*(agg->d_m));
            } else { // electric charges preservation
                int total_charge=0;
                for (const auto &sph : agg->myspheres) {
                    total_charge += sph->electric_charge;
                }
                agg->electric_charge = total_charge;
            }
        }
    }
    return have_splitted;
}
void Aggregate::remove_sphere(const size_t &id) noexcept {
    for (size_t local_id = 0; local_id < n_spheres; local_id++) {
        if (myspheres[local_id]->index_in_storage == id) {
            myspheres.list.erase(myspheres.list.begin() + long(local_id));
            break;
        }
    }
    external_storage->spheres.remove(id);
    n_spheres = myspheres.size();
    alpha_vs_extreme = 1.0/static_cast<double>(n_spheres);
    if (n_spheres > 0) {
        std::array<double, 3> refpos = myspheres[0]->get_position();
        for (const auto &sph : myspheres) {
            // change the relative position of the new aggregate
            sph->set_relative_position(sph->get_position() - refpos);
        }
        // we have to recompute all the caracteristic of this new aggregate
        //update_distances_and_overlapping();
        update();
        // we have to move all the spheres (periodicity)
        refpos = get_position() - get_relative_position();
        for (const auto &sph : myspheres) {
            sph->set_position(refpos + sph->get_relative_position());
        }
    } // else it should be deleted ASAP
}
// This function below convert one aggregate to a mass equivalent sphere
//void Aggregate::agg_to_sphere() noexcept {
//    // 1. Select one sphere to keep and move it to the agg. center of mass
//    std::array<double, 3> mypos{{*x, *y, *z}};
//    myspheres[0].set_position(mypos);
//    // 2. Increase the radius of the selected sphere - mass conservation
//    double sph_0_radius = 0.5 * pow(6 * (*agregat_volume) / _pi, 1.0 / 3.0);
//    myspheres[0].set_radius(sph_0_radius);
//    // 3. Update myspheres, delete additional ones
//    // TODO: this part is pending
//    // 4. Update aggregate properties according to selected sphere
//    n_spheres = myspheres.size();
//    update_distances_and_overlapping();
//    update();
//}
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
    std::cout << "    C/H ratio         : " << *CH_ratio << std::endl;
    std::cout << "    bulk density      : " << bulk_density << std::endl;
    std::cout << "    Position          : " << *x << " " << *y << " " << *z << std::endl;
    std::cout << "    Proper time       : " << *proper_time << std::endl;
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
// In the function below we do both, update distances and calc. n_c_avg and c_ov_avg
void Aggregate::update_distances_and_overlapping() noexcept {
    *overlapping = *coordination_number = 0.0; // average overlapping and coord. numbers calculated if needed
    double c_ij(0);
    int intersections(0);
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
            double dist = relative_distance(*myspheres[i], *myspheres[j]);
#ifdef FULL_INTERNAL_DISTANCES
            distances[i][j] = dist;
            // distances are symetric !
            distances[j][i] = dist;
#endif
            // if in contact (with some tolerence)
            if (dist <= (1. + _COORDINATION_EPSILON) * (myspheres[i]->get_radius() + myspheres[j]->get_radius())) {
#ifndef FULL_INTERNAL_DISTANCES
                distances[i][j] = dist;
                // distances are symetric !
                distances[j][i] = dist;
#endif
                // calculate overlapping coefficient
                c_ij = (myspheres[i]->get_radius() + myspheres[j]->get_radius() - dist)
                       / (myspheres[i]->get_radius() + myspheres[j]->get_radius());
                *overlapping += 2.0*c_ij;
                intersections += 2;
            }
        }
    }
    if (intersections > 0) {
        *overlapping /= static_cast<double>(intersections);
        *coordination_number = static_cast<double>(intersections) / static_cast<double>(n_spheres);
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
