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
#ifndef INCLUDE_AGGREGATS_AGGREGAT_LIST_HPP
#define INCLUDE_AGGREGATS_AGGREGAT_LIST_HPP 1
#include "aggregat.hpp"
#include "list_storage/list_storage.hpp"
#include "physical_model/physical_model.hpp"
#include "spheres/sphere_list.hpp"
#include "verlet/verlet.hpp"
#include <gsl/gsl>


namespace mcac {
class AggregatList : public ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate> {
    friend class Aggregate;

private:
    PhysicalModel *physicalmodel;
    double maxradius;
    double avg_npp;
    double max_time_step;
    std::vector<size_t> index_sorted_time_steps;
    std::vector<double> cumulative_time_steps;
    std::vector<double>::iterator ptr_deb;
    std::vector<double>::iterator ptr_fin;
    std::unique_ptr<ThreadedIO> writer;
    size_t last_saved;
public:
    SphereList spheres;
    Verlet verlet;
    /* getters */
    double get_total_volume() const;
    double get_avg_npp() const;
    double get_max_time_step() const;
    double get_time_step(double max) const;
    size_t pick_random() const;
    size_t pick_last() const;
    /* modifiers */
    using ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate>::add;
    void add(size_t);
    void refresh();
    void sort_time_steps(double factor);
    void duplication();
    bool split();
    bool split_individual(const size_t numagg);
    bool merge(AggregateContactInfo contact_info);
    InterPotentialRegime check_InterPotentialRegime(AggregateContactInfo contact_info);
    /* other */
    bool croissance_surface(double dt);
    bool croissance_surface_individual(const double dt, const size_t index);
    std::tuple<bool, double, double, double> get_instantaneous_fractal_law() const;
    AggregateContactInfo distance_to_next_contact(const size_t source,
                                                  const std::array<double, 3> &direction,
                                                  const double distance) const;
    std::vector<size_t> get_neighborhood(const size_t source,
                                         const std::array<double, 3> &direction,
                                         const double distance) const;
    std::multimap<double, size_t> filter_neighborhood(const size_t moving_aggregate,
                                                      const std::array<double, 3> &direction,
                                                      const std::vector<size_t> &neighborhood,
                                                      const double distance) const;
    bool test_free_space(std::array<double, 3> pos, double diameter) const;
    /* I/O */
    void save();
    auto get_data() const;
    std::vector<double> format_position() const;
    std::vector<double> format_rg() const;
    std::vector<int> format_n_spheres() const;
    std::vector<double> format_f_agg() const;
    std::vector<double> format_lpm() const;
    std::vector<double> format_time_step() const;
    std::vector<double> format_rmax() const;
    std::vector<double> format_agregat_volume() const;
    std::vector<double> format_agregat_surface() const;
    std::vector<double> format_coordination_number() const;
    std::vector<double> format_proper_time() const;
    std::vector<double> format_overlapping() const;
    std::vector<double> format_d_m() const;
    std::vector<int> format_electric_charge() const;
    std::vector<int> format_label() const;
    void remove(const size_t &id) noexcept;
    void remove_sphere(const size_t &id) noexcept;
    /* Storage specific */
private:
    void setpointers();
public:
    /** Default constructor in local storage */
    explicit AggregatList(PhysicalModel *);
    /** Destructor */
    ~AggregatList() noexcept;
    /** Copy constructor */
    AggregatList(const AggregatList &other) noexcept = delete;
    /** Move constructor */
    AggregatList(AggregatList &&) noexcept = delete;
    /** Copy assignment operator */
    AggregatList &operator=(const AggregatList &other) noexcept = delete;
    /** Move assignment operator */
    AggregatList &operator=(AggregatList &&other) noexcept = delete;
};
}// namespace mcac


#endif //INCLUDE_AGGREGATS_AGGREGAT_LIST_HPP
