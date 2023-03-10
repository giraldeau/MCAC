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
#ifndef INCLUDE_PHYSICAL_MODEL_PHYSICAL_MODEL_HPP
#define INCLUDE_PHYSICAL_MODEL_PHYSICAL_MODEL_HPP 1
#include "constants.hpp"
#include "io/xmf_includes.hpp"
#include "physical_model_flame_coupling.hpp"
#include "physical_model_interpotential.hpp"
#include <array>
#include <ctime>
#include <experimental/filesystem>

namespace mcac {
class PhysicalModel {
  public:
    double fractal_dimension, fractal_prefactor;
    double flux_surfgrowth, u_sg;             // surface growth molecular flux [kg/m2/s], SG velocity u_sg [m/s]
    double flux_nucleation, nucleation_accum; // Nucleation flux [#/m^3-flame/s], accumulated nucleation [#]
    double pressure, temperature, gaz_mean_free_path, viscosity, density;
    double mean_diameter, dispersion_diameter;
    double mass_nuclei, mean_diameter_nucleation, dispersion_diameter_nucleation;
    double mean_massic_radius, friction_exponnant;
    double time;
    double volume_fraction, box_length, box_volume, aggregate_concentration, monomer_concentration;
    double total_surface_concent, total_volume_concent;
    double rp_min_oxid; // Minimum d_pp below which the particle disappear
    size_t n_verlet_divisions;
    PickMethods pick_method;
    VolSurfMethods volsurf_method;
    size_t n_monomeres;
    size_t n_time_per_file;
    MonomeresInitialisationMode monomeres_initialisation_type;
    size_t n_iter_without_event;
    size_t last_timestep_written;
    clock_t cpu_start;
    clock_t cpu_last_event;
    double cpu_limit, cpu_event_limit, physical_time_limit;
    double write_Delta_t;
    int mean_monomere_per_aggregate_limit;
    size_t number_of_aggregates_limit;
    int n_iter_without_event_limit;
    int random_seed;
    size_t write_events_frequency;
    size_t write_between_event_frequency;
    size_t full_aggregate_update_frequency;
    bool finished_by_flame;
    std::experimental::filesystem::path output_dir;
    std::string flame_file, interpotential_file;
    bool with_domain_duplication;
    bool with_domain_reduction;
    bool with_nucleation;
    bool with_collisions;
    bool with_surface_reactions;
    bool with_flame_coupling;
    bool enforce_volume_fraction;
    bool individual_surf_reactions;
    bool with_potentials, with_external_potentials, with_dynamic_random_charges;
    bool with_electric_charges;
    bool with_maturity;
    FlameCoupling flame;
    Interpotential intpotential_info;

    explicit PhysicalModel(const std::string &fichier_param);
    [[gnu::pure]] double cunningham(double r) const;
    [[gnu::pure]] double random_diameter(double mean_diameter, double dispersion_diameter) const;
    [[gnu::pure]] double grow(double r, double dt) const;
    [[gnu::pure]] double friction_exponent(double sphere_radius) const;
    [[gnu::pure]] double friction_coeff(double aggregate_volume, double sphere_volume, double sphere_radius) const;
    [[gnu::pure]] double diffusivity(double) const;
    [[gnu::pure]] static double relax_time(double masse, double);
    [[gnu::pure]] double mobility_diameter(double aggregate_volume, double sphere_volume, double sphere_radius) const;
    [[gnu::pure]] int get_random_charge(double mobility_diameter) const;
    void print() const;
    void update(size_t n_aggregates, size_t n_monomers, double new_total_volume, double new_total_surface) noexcept;
    void nucleation(double dt) noexcept;
    void update_from_flame();
    void update_temperature(double new_temperature) noexcept;
    [[gnu::pure]] bool finished(size_t number_of_aggregates, double mean_monomere_per_aggregate) const;
    [[gnu::pure]] bool time_to_write(size_t total_events);
    XMF_OUTPUT xmf_write() const;
};

[[gnu::const]] std::experimental::filesystem::path extract_path(const std::string &file);
[[gnu::const]] inline double periodic_distance(double dist, double dim) {
    double periodic_distance(dist);
    double half_dim(0.5 * dim);
    while (periodic_distance < -half_dim) {
        periodic_distance += dim;
    }
    while (periodic_distance >= half_dim) {
        periodic_distance -= dim;
    }
    return periodic_distance;
}
[[gnu::const]] inline std::array<double, 3> periodic_distance(std::array<double, 3> dist, double dim) {
    return {periodic_distance(dist[0], dim), periodic_distance(dist[1], dim), periodic_distance(dist[2], dim)};
}
[[gnu::const]] inline double periodic_position(double position, double dim) {
    double periodic_position(position);
    while (periodic_position < 0) {
        periodic_position += dim;
    }
    while (periodic_position >= dim) {
        periodic_position -= dim;
    }
    return periodic_position;
}
[[gnu::const]] inline std::array<double, 3> periodic_position(std::array<double, 3> position, double dim) {
    return {
        periodic_position(position[0], dim), periodic_position(position[1], dim), periodic_position(position[2], dim)};
}
} // namespace mcac

#endif // INCLUDE_PHYSICAL_MODEL_PHYSICAL_MODEL_HPP
