#ifndef INCLUDE_PHYSICAL_MODEL_PHYSICAL_MODEL_HPP_
#define INCLUDE_PHYSICAL_MODEL_PHYSICAL_MODEL_HPP_ 1
#include "constants.hpp"
#include "io/xmf_includes.hpp"
#include <ctime>
#include <array>
#include <experimental/filesystem>


namespace MCAC {
class PhysicalModel {
public:
    double fractal_dimension, fractal_prefactor;
    double x_surfgrowth, coeff_b;
    double a_surfgrowth;
    double pressure, temperature, gaz_mean_free_path, viscosity, density;
    double mean_diameter, dispersion_diameter;
    double mean_massic_radius, friction_exponnant;
    double time;
    double volume_fraction, box_lenght;
    size_t n_verlet_divisions;
    size_t n_monomeres;
    size_t n_time_per_file;
    int monomeres_initialisation_type, n_iter_without_contact;
    clock_t cpu_start;
    double cpu_limit, physical_time_limit;
    size_t mean_monomere_per_aggregate_limit;
    size_t number_of_aggregates_limit;
    int n_iter_without_contact_limit;
    std::experimental::filesystem::path output_dir;
public:
    explicit PhysicalModel(const std::string &fichier_param) noexcept;
    [[gnu::pure]] double cunningham(double r) const;
    [[gnu::pure]] double grow(double r, double dt) const;
    [[gnu::pure]] double friction_coeff(size_t npp) const;
    [[gnu::pure]] double friction_coeff2(double rgg) const;
    [[gnu::pure]] double diffusivity(double) const;
    [[gnu::pure]] double relax_time(double masse, double) const;
    [[gnu::pure]] void print() const;
    [[gnu::pure]] bool finished(size_t number_of_aggregates, double mean_monomere_per_aggregate) const;
    XMF_OUTPUT xmf_write() const;
};

[[gnu::const]] double inverfc(double p);
[[gnu::const]] double inverf(double p);
[[gnu::const]] std::experimental::filesystem::path extract_path(const std::string &file);
[[gnu::const]] std::tuple<bool, double, double, double> linreg(const std::vector<double> &x,
                                                               const std::vector<double> &y);
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
    return {periodic_distance(dist[0], dim),
            periodic_distance(dist[1], dim),
            periodic_distance(dist[2], dim)};
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
[[gnu::const]] inline std::array<double, 3> periodic_position(std::array<double, 3> position,
                                                              double dim) {
    return {periodic_position(position[0], dim),
            periodic_position(position[1], dim),
            periodic_position(position[2], dim)};
}
}  // namespace MCAC


#endif //INCLUDE_PHYSICAL_MODEL_PHYSICAL_MODEL_HPP_

