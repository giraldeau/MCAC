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
#ifndef INCLUDE_AGGREGATS_AGGREGAT_HPP
#define INCLUDE_AGGREGATS_AGGREGAT_HPP 1
#include "constants.hpp"
#include "elem_storage/elem_storage.hpp"
#include "spheres/sphere_list.hpp"
#include "sbl/volume_surface.hpp"
#include "tools/contact_info.hpp"
#include <list>
#include "physical_model/physical_model_interpotential.hpp"


namespace mcac {
class AggregatList;

class Verlet;

class Aggregate :
    public ElemStorage<AggregatesFields::AGGREGAT_NFIELDS, AggregatList> {
    friend class AggregatList;

private:
    double *rg;                     // Gyration Radius
    double *f_agg;                  // Friction coeff
    double *lpm;                    // Apparent Mean Free Path
    double *time_step;              // Time to move along lpm
    double *rmax;                   // Radius of the sphere containing Agg
    double *agregat_volume;         // Etimation of the aggregate's volume
    double *agregat_surface;        // Estimation of the aggregate's surface
    double *x, *y, *z;              // position of the gravity center
    double *rx, *ry, *rz;           // position of the gravity center
    double *proper_time;            // Proper time of the aggregate
    double *dp;                     // Mean monomer diameters
    double *dg_over_dp;
    double *overlapping;            // Average overlapping coefficient
    double *coordination_number;    // Average coordination number
    double *d_m;                     // Mobility diameter
    int electric_charge;            // Total agg. electric charge
    size_t n_spheres;               // Number of spheres
    size_t label;                   // Uniq label of the aggregat
    double bulk_density;

#ifdef FULL_INTERNAL_DISTANCES
    std::vector<std::vector<double> > distances;
#else
    std::vector<std::unordered_map<size_t, double> > distances;
#endif
    std::vector<double> distances_center;
    std::vector<double> volumes;
    std::vector<double> surfaces;
    /*
    bool padding2,padding3,padding4;
    int padding1;
    */
    std::array<size_t, 3> compute_index_verlet() noexcept;
    void update_verlet() noexcept;
    void update_distances_and_overlapping() noexcept;
    double internal_sphere_distance(size_t i, size_t j) const noexcept;
    double sphere_distance_center(size_t i) const noexcept;
public:
    const PhysicalModel *physicalmodel;
    Verlet *verlet;
    std::array<size_t, 3> index_verlet;
    SphereList myspheres;
    /* getters */
    const double &get_rg() const noexcept;
    const double &get_dp() const noexcept;
    const double &get_f_agg() const noexcept;
    const double &get_lpm() const noexcept;
    const double &get_time_step() const noexcept;
    const double &get_rmax() const noexcept;
    const double &get_agregat_volume() const noexcept;
    const double &get_agregat_surface() const noexcept;
    std::array<double, 3> get_position() const noexcept;
    std::array<double, 3> get_relative_position() const noexcept;
    std::array<size_t, 3> get_verlet_index() const noexcept;
    const double &get_proper_time() const noexcept;
    const size_t &size() const noexcept;
    const size_t &get_label() const noexcept;
    const int &get_electric_charge() const noexcept;
    /* modifiers */
    void decrease_label() noexcept;
    void set_verlet(Verlet *) noexcept;
    void unset_verlet() noexcept;
    void set_proper_time(double newtime) noexcept;
    void set_bulk_density() noexcept;
    void time_forward(double deltatemps) noexcept;
    void set_position(const std::array<double, 3> &position) noexcept;
    void translate(std::array<double, 3> vector) noexcept;
    bool croissance_surface(double dt);
    //    void init();
    void init(size_t new_label,
              size_t sphere_index);
    void init(const PhysicalModel &,
              SphereList *,
              Verlet *,
              double new_time,
              size_t new_label,
              size_t sphere_index,
              const std::array<double, 3> &position,
              double sphere_diameter) noexcept;
    void update_partial() noexcept;
    void update() noexcept;
    void compute_volume_surface();
    void compute_mass_center() noexcept;
    void compute_max_radius() noexcept;
    void compute_giration_radius() noexcept;
    /* other */
    bool merge(std::shared_ptr<Aggregate> other, AggregateContactInfo contact_info) noexcept;
    bool split();
    void remove_sphere(const size_t &id) noexcept;
    // void agg_to_sphere() noexcept;
    void print() const noexcept;
    /* Storage specific */
private:
    void setpointers() noexcept;
public:
    /** Default constructor in local storage */
    Aggregate() noexcept;
    //    explicit Aggregate(PhysicalModel &);
    /** Constructor with external storage */
    Aggregate(AggregatList *, size_t newlabel) noexcept;
    /** Destructor */
    ~Aggregate() noexcept;
    /** Copy constructor */
    Aggregate(const Aggregate &, AggregatList *, size_t newlabel) noexcept;
    Aggregate(const Aggregate &) noexcept = delete;
    /** Move constructor */
    Aggregate(Aggregate &&) noexcept = delete;
    /** Copy assignment operator */
    Aggregate &operator=(const Aggregate &other) noexcept = delete;
    /** Move assignment operator */
    Aggregate &operator=(Aggregate &&other) noexcept = delete;
};
}  // namespace mcac


#endif //INCLUDE_AGGREGATS_AGGREGAT_HPP
