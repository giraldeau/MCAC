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
#include <list>


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
    double *time;                   // Proper time of the aggregate
    double *dp;
    double *dg_over_dp;
    size_t n_spheres;               // Number of spheres
    size_t label;                   // Uniq label of the aggregat

    std::vector<std::list<std::pair<size_t, double> > > distances;
    std::vector<double> distances_center;
    std::vector<double> volumes;
    std::vector<double> surfaces;
    /*
    bool padding2,padding3,padding4;
    int padding1;
    */
    void update_verlet_index() noexcept;
    void update_distances() noexcept;
    double internal_sphere_distance(size_t i, size_t j) const noexcept;
    double sphere_distance_center(size_t i) const noexcept;
public:
    const PhysicalModel *physicalmodel;
    Verlet *verlet;
    std::array<size_t, 3> index_verlet;
    SphereList myspheres;
    /* getters */
    double get_rg() const noexcept;
    double get_f_agg() const noexcept;
    double get_lpm() const noexcept;
    double get_time_step() const noexcept;
    double get_rmax() const noexcept;
    double get_agregat_volume() const noexcept;
    double get_agregat_surface() const noexcept;
    std::array<double, 3> get_position() const noexcept;
    std::array<double, 3> get_relative_position() const noexcept;
    std::array<size_t, 3> get_verlet_index() const noexcept;
    double get_time() const noexcept;
    size_t size() const noexcept;
    size_t get_label() const noexcept;
    /* modifiers */
    void decrease_label() noexcept;
    void set_verlet(Verlet *) noexcept;
    void unset_verlet() noexcept;
    void set_time(double newtime) noexcept;
    void time_forward(double deltatemps) noexcept;
    void set_position(const std::array<double, 3> &position) noexcept;
    void translate(std::array<double, 3> vector) noexcept;
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
    void update() noexcept;
    void compute_volume() noexcept;
    void compute_mass_center() noexcept;
    void compute_max_radius() noexcept;
    void compute_giration_radius() noexcept;
    /* other */
    void merge(Aggregate *) noexcept;
    bool split();
    void remove(const size_t &id) noexcept ;
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
    Aggregate(const Aggregate &, AggregatList *) noexcept;
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
