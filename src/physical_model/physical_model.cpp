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
#include "constants.hpp"
#include "physical_model/physical_model.hpp"
#include "tools/tools.hpp"
#include "exceptions.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <inipp.h>
#include <iostream>


namespace fs = std::experimental::filesystem;
namespace mcac {
PhysicalModel::PhysicalModel(const std::string &fichier_param) :
    fractal_dimension(1.4),
    fractal_prefactor(1.8),
    flux_surfgrowth(0.),
    u_sg(0.),
    pressure(_pressure_ref),
    temperature(_temperature_ref),
    gaz_mean_free_path(_fluid_mean_free_path_ref),
    viscosity(_viscosity_ref),
    density(1800.),
    mean_diameter(30.),
    dispersion_diameter(1.25),
    mean_massic_radius(0.),
    friction_exponnant(0.),
    time(0.),
    volume_fraction(1e-3),
    box_lenght(0.),
    aggregate_concentration(0.0),
    n_verlet_divisions(10),
    pick_method(PickMethods::PICK_RANDOM),
    n_monomeres(2500),
    n_time_per_file(10),
    monomeres_initialisation_type(MonomeresInitialisationMode::NORMAL_INITIALISATION),
    n_iter_without_event(0),
    cpu_start(),
    cpu_limit(-1),
    physical_time_limit(-1),
    mean_monomere_per_aggregate_limit(-1),
    number_of_aggregates_limit(1),
    n_iter_without_event_limit(-1),
    write_between_event_frequency(100),
    full_aggregate_update_frequency(1),
    output_dir("MCAC_output"),
    with_collisions(true),
    with_surface_reactions(false) {
    std::string default_str;
    // read the config file
    inipp::Ini<char> ini;
    std::ifstream is(fichier_param);
    ini.parse(is);
    ini.interpolate();
    std::cout << "Directory, input parameters: " << extract_path(fichier_param) << std::endl;
    // monomeres
    inipp::extract(ini.sections["monomers"]["number"], n_monomeres);
    inipp::extract(ini.sections["monomers"]["density"], density);
    inipp::extract(ini.sections["monomers"]["dispersion_diameter"], dispersion_diameter);
    inipp::extract(ini.sections["monomers"]["mean_diameter"], mean_diameter);
    default_str = resolve_monomeres_initialisation_mode(monomeres_initialisation_type);
    inipp::extract(ini.sections["monomers"]["initialisation_mode"], default_str);
    monomeres_initialisation_type = resolve_monomeres_initialisation_mode(default_str);
    if (monomeres_initialisation_type == MonomeresInitialisationMode::INVALID_INITIALISATION) {
        throw InputError("Monomere initialisation mode unknown: " + default_str);
    }
    // environment
    inipp::extract(ini.sections["environment"]["initial_time"], time);
    inipp::extract(ini.sections["environment"]["volume_fraction"], volume_fraction);
    inipp::extract(ini.sections["environment"]["temperature"], temperature);
    inipp::extract(ini.sections["environment"]["pressure"], pressure);
    inipp::extract(ini.sections["environment"]["fractal_prefactor"], fractal_prefactor);
    inipp::extract(ini.sections["environment"]["fractal_dimension"], fractal_dimension);
    // surface_growth
    inipp::extract(ini.sections["surface_growth"]["flux_surfgrowth"], flux_surfgrowth);
    inipp::extract(ini.sections["surface_growth"]["full_aggregate_update_frequency"], full_aggregate_update_frequency);
    // limits
    inipp::extract(ini.sections["limits"]["number_of_aggregates"], number_of_aggregates_limit);
    inipp::extract(ini.sections["limits"]["n_iter_without_event"], n_iter_without_event_limit);
    inipp::extract(ini.sections["limits"]["cpu"], cpu_limit);
    inipp::extract(ini.sections["limits"]["physical_time"], physical_time_limit);
    inipp::extract(ini.sections["limits"]["mean_monomere_per_aggregate"], mean_monomere_per_aggregate_limit);
    // numerics
    inipp::extract(ini.sections["numerics"]["with_collisions"], with_collisions);
    inipp::extract(ini.sections["numerics"]["with_surface_reactions"], with_surface_reactions);
    inipp::extract(ini.sections["numerics"]["n_verlet_divisions"], n_verlet_divisions);
    default_str = resolve_pick_method(pick_method);
    inipp::extract(ini.sections["numerics"]["pick_method"], default_str);
    pick_method = resolve_pick_method(default_str);
    // output
    inipp::extract(ini.sections["output"]["output_dir"], output_dir);
    inipp::extract(ini.sections["output"]["n_time_per_file"], n_time_per_file);
    inipp::extract(ini.sections["output"]["write_between_event_frequency"], write_between_event_frequency);
    output_dir = extract_path(fichier_param) / output_dir;
    if (fs::exists(output_dir) && fs::is_directory(output_dir)) {
        std::string answer;
        std::string valid_answer = "rRIiAa";
        do {
            std::cout << "folder " << output_dir << " already exists." << std::endl;
            std::cout << " (R)emove, (I)gnore, (A)bandon " << std::endl;
            std::cin >> answer;
        } while (std::find(std::begin(valid_answer), std::end(valid_answer), answer[0]) == std::end(valid_answer));
        if (answer == "r" || answer == "R") {
            fs::remove_all(output_dir);
        }
        if (answer == "a" || answer == "A") {
            throw AbandonError();
        }
    }
    if (!fs::exists(output_dir)) {
        if (!fs::create_directory(output_dir)) {
            throw IOError("Error creating directory " + output_dir.string());
        }
    } else {
        if (!fs::is_directory(output_dir)) {
            throw InputError("Error output_dir is not a directory " + output_dir.string());
        }
    }
    // secondary variables
    if (monomeres_initialisation_type == MonomeresInitialisationMode::NORMAL_INITIALISATION) {
        box_lenght = mean_diameter * 1E-9 *
                     std::pow(static_cast<double>(n_monomeres) * _pi / 6. / volume_fraction
                              * (1. + 3. * std::pow(dispersion_diameter / mean_diameter, 2)),
                              1. / 3.);
        mean_massic_radius = 0.5 * 1E-9 * (std::pow(mean_diameter, 4)
                                           + 6 * std::pow(mean_diameter, 2) * std::pow(dispersion_diameter, 2)
                                           + 3 * std::pow(dispersion_diameter, 4))
                             / (std::pow(mean_diameter, 3) + 3 * mean_diameter * std::pow(dispersion_diameter, 2));
    } else if (monomeres_initialisation_type == MonomeresInitialisationMode::LOG_NORMAL_INITIALISATION) {
        box_lenght = mean_diameter * 1E-9 *
                     std::pow(static_cast<double>(n_monomeres) * _pi / 6. / volume_fraction
                              * std::exp(9. / 2. * std::pow(std::log(dispersion_diameter), 2)),
                              1. / 3.);
        mean_massic_radius = 0.5 * mean_diameter * 1E-9
                             * std::exp(1.5 * std::pow(std::log(dispersion_diameter), 2));  // Hatch-Choate
    } else {
        throw InputError("Monomere initialisation mode unknown");
    }
    update_temperature(temperature);
    u_sg = flux_surfgrowth / density;     // Surface growth velocity [m/s], u_sg=dr_p/dt
    // Particle number concentration
    aggregate_concentration = static_cast<double>(n_monomeres) / std::pow(box_lenght, 3);
    print();
    std::ofstream os(output_dir / "params.ini");
    ini.generate(os);
    cpu_start = clock();
}
[[gnu::pure]] bool PhysicalModel::finished(size_t number_of_aggregates, double mean_monomere_per_aggregate) const {
    if (number_of_aggregates <= number_of_aggregates_limit) {
        std::cout << "We reach the AggMin condition" << std::endl << std::endl;
        return true;
    }
    if (n_iter_without_event_limit > 0
        && n_iter_without_event >= static_cast<size_t>(n_iter_without_event_limit)) {
        std::cout << "We reach the WaitLimit condition" << std::endl << std::endl;
        return true;
    }
    if (cpu_limit > 0) {
        clock_t current_cpu = clock();
        double elapse = double(current_cpu - cpu_start) / CLOCKS_PER_SEC;
        if (elapse >= cpu_limit) {
            std::cout << "We reach the CPULimit condition" << std::endl << std::endl;
            return true;
        }
    }
    if (physical_time_limit > 0
        && time >= physical_time_limit) {
        std::cout << "We reach the Maximum physical time condition " << time << "/" << physical_time_limit << std::endl;
        return true;
    }
    if (mean_monomere_per_aggregate_limit > 0
        && mean_monomere_per_aggregate >= mean_monomere_per_aggregate_limit) {
        std::cout << "We reach the NPP_avg_limit condition " << mean_monomere_per_aggregate << "/"
                  << mean_monomere_per_aggregate_limit << std::endl;
        return true;
    }
    return false;
}
void PhysicalModel::print() const {
    std::cout << "PARTICLES PROPERTIES:" << std::endl
              << " density  : " << density << " (kg/m^3)" << std::endl
              << " Dpm      : " << mean_diameter << " (nm)" << std::endl
              << " sigmaDpm : " << dispersion_diameter << " (- Lognormal, nm Normal)" << std::endl
              << " dfe      : " << fractal_dimension << " (-)" << std::endl
              << " kfe      : " << fractal_prefactor << " (-)" << std::endl
              << " rpeqmass : " << mean_massic_radius << " (m)" << std::endl
              << " gamma_   : " << friction_exponnant << " (-)" << std::endl
              << std::endl
              << "FLUID PROPERTIES:" << std::endl
              << " Pressure    : " << pressure << " (Pa)" << std::endl
              << " Temperature : " << temperature << " (K)" << std::endl
              << " viscosity   : " << viscosity << " (kg/m*s)" << std::endl
              << " lambda      : " << gaz_mean_free_path << " (m)" << std::endl
              << std::endl
              << "SIMULATION OPTIONS:" << std::endl
              << " Initial aggregate concentration : " << aggregate_concentration << " (#/m^3)" << std::endl
              << " Initial Nagg                    : " << n_monomeres << " (-)" << std::endl
              << " Box size                        : " << box_lenght << " (m)" << std::endl
              << " FV                              : " << volume_fraction << " (-)" << std::endl;
    if (with_collisions) {
        std::cout << " With collisions" << std::endl;
    } else{
        std::cout << " Without collision" << std::endl;
    }
    if (with_surface_reactions) {
        std::cout << " With surface reations" << std::endl
                  << "  flux_surfgrowth                 : " << flux_surfgrowth << " (kg/m^2/s)" << std::endl
                  << "  u_sg                            : " << u_sg << " (m/s)" << std::endl;
    } else {
        std::cout << " Without surface reations" << std::endl;
    }
    std::cout << std::endl
              << "Options for Pysical model: " << std::endl
              << " Initialisation mode : " << resolve_monomeres_initialisation_mode(monomeres_initialisation_type)
              << std::endl
              << " Pick method         : " << resolve_pick_method(pick_method) << std::endl
              << std::endl
              << "Ending simulation when:" << std::endl
              << " - There are " << number_of_aggregates_limit << " aggregats left or less" << std::endl;
    if (n_iter_without_event_limit > 0) {
        std::cout << " - It has been " << n_iter_without_event_limit << " iterations without collision" << std::endl;
    }
    if (cpu_limit > 0) {
        std::cout << " - The simulations is running for more than " << cpu_limit << " seconds" << std::endl;
    }
    if (physical_time_limit > 0) {
        std::cout << " - The residence time is larger than " << physical_time_limit << " seconds" << std::endl;
    }
    if (mean_monomere_per_aggregate_limit > 0) {
        std::cout << " - The average Npp per aggregate reach " << mean_monomere_per_aggregate_limit << " monomers"
                  << std::endl;
    }
    std::cout << std::endl;
}
void PhysicalModel::update(size_t n_aggregates, double total_volume) noexcept {
    aggregate_concentration = static_cast<double>(n_aggregates) / std::pow(box_lenght, 3);
    volume_fraction = total_volume / std::pow(box_lenght, 3);
}
//#####################################################################################################################
void PhysicalModel::update_temperature(double new_temperature) noexcept {
    temperature = new_temperature;
    viscosity = _viscosity_ref
                * (_sutterland_interpolation_constant + _temperature_ref)
                / (_sutterland_interpolation_constant + temperature)
                * std::pow(temperature / _temperature_ref, 1.5); // Schlichting 1979
    gaz_mean_free_path = _fluid_mean_free_path_ref
                         * (_pressure_ref / pressure)
                         * (temperature / _temperature_ref)
                         * (1. + _sutterland_interpolation_constant / _temperature_ref)
                         / (1. + _sutterland_interpolation_constant / temperature); // Willeke 1976
    // Yon et al. Journal of Aerosol Science vol.87 2015 p.2837 eq.7
    friction_exponnant = 0.689 * (1. + std::erf(((gaz_mean_free_path / mean_massic_radius) + 4.454) / 10.628));
}

//#####################################################################################################################

[[gnu::pure]]  double PhysicalModel::cunningham(double r) const {
    double a = 1.142;
    double b = 0.558;
    double c = 0.999;
    return 1.0 + a * gaz_mean_free_path / r
           + b * gaz_mean_free_path / r * std::exp(-c * r / gaz_mean_free_path);
}


//############################# Fonctions pour le calcul du diamètre de mobilité ################################

/*
 Fonction permettant de retrouver le rayon de mobilité en régime transitoire
 On obtient le bon rayon de mobilité lorsque la fonction retourne 0
*/

[[gnu::pure]]  double PhysicalModel::grow(double r, double dt) const {
    // Based on the molecular flux
    return r + u_sg * dt;
}
[[gnu::pure]] double PhysicalModel::friction_coeff(double aggregate_volume,
                                                   double sphere_volume,
                                                   double sphere_radius) const {
    double friction_exp = 0.689 * (1. + std::erf(((gaz_mean_free_path / sphere_radius) + 4.454) / 10.628));
    double cc = cunningham(sphere_radius);
    return (6. * _pi * viscosity * sphere_radius / cc)
           * std::pow(aggregate_volume / sphere_volume, friction_exp / fractal_dimension);
}
[[gnu::pure]] double PhysicalModel::diffusivity(double f_agg) const {
    return _boltzmann * temperature / f_agg;
}
[[gnu::pure]] double PhysicalModel::relax_time(double masse, double f_agg) {
    return masse / f_agg;
}
[[gnu::const]] double inverfc(double p) {
    double x;
    double t;
    double pp;
    if (p >= 2.) {
        return -100.;
    }
    if (p <= 0.0) {
        return 100.;
    }
    pp = (p < 1.0) ? p : 2. - p;
    t = std::sqrt(-2. * std::log(pp / 2.));
    x = -0.70711 * ((2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t);
    for (int j = 0; j < 2; j++) {
        double err = std::erfc(x) - pp;
        x += err / (1.12837916709551257 * std::exp(-(x * x)) - x * err);
    }
    return (p < 1.0 ? x : -x);
}
[[gnu::const]] double inverf(double p) {
    return inverfc(1. - p);
}
[[gnu::const]] fs::path extract_path(const std::string &filename) {
    fs::path path = filename;
    fs::path fullpath = fs::absolute(path);
    if (!fs::exists(fullpath)) {
        throw InputError("File does not exist");
    }
    fs::path parentpath = fullpath.parent_path();
    return parentpath;
}
}  // namespace mcac

