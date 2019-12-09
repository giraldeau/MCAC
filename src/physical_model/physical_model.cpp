#include "physical_model/physical_model.hpp"
#include "constants.hpp"
#include "tools/tools.hpp"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <inipp.h>


namespace fs = std::experimental::filesystem;
using namespace std;
namespace MCAC {
PhysicalModel::PhysicalModel(const string &fichier_param) noexcept:
    dfe(1.4),
    kfe(1.8),
    x_surfgrowth(2),
    pressure(_pressure_ref),
    temperature(_temperature_ref),
    density(1800),
    mean_diameter(30),
    dispersion_diameter(1.25),
    volume_fraction(1e-3),
    n_verlet_divisions(10),
    n_monomeres(2500),
    n_time_per_file(10),
    monomeres_initialisation_type(MonomeresInitialisationMode::NORMAL_INITIALISATION),
    n_iter_without_contact(0),
    cpu_limit(-1),
    physical_time_limit(-1),
    mean_monomere_per_aggregate_limit(-1),
    number_of_aggregates_limit(1),
    n_iter_without_contact_limit(-1) {
    // read the config file
    inipp::Ini<char> ini;
    std::ifstream is(fichier_param);
    ini.parse(is);
    ini.interpolate();
    // monomeres
    inipp::extract(ini.sections["monomeres"]["number"], n_monomeres);
    inipp::extract(ini.sections["monomeres"]["density"], density);
    inipp::extract(ini.sections["monomeres"]["dispersion_diameter"], dispersion_diameter);
    inipp::extract(ini.sections["monomeres"]["mean_diameter"], mean_diameter);
    string default_mode = "normal";
    inipp::extract(ini.sections["monomeres"]["initialisation_mode"], default_mode);
    monomeres_initialisation_type = resolveMonomeresInitialisationMode(default_mode);
    // environment
    inipp::extract(ini.sections["environment"]["volume_fraction"], volume_fraction);
    inipp::extract(ini.sections["environment"]["temperature"], temperature);
    inipp::extract(ini.sections["environment"]["pressure"], pressure);
    // surface_growth
    inipp::extract(ini.sections["surface_growth"]["kfe"], kfe);
    inipp::extract(ini.sections["surface_growth"]["dfe"], dfe);
    inipp::extract(ini.sections["surface_growth"]["coeff_b"], coeff_b);
    inipp::extract(ini.sections["surface_growth"]["x_surfgrowth"], x_surfgrowth);
    // limits
    inipp::extract(ini.sections["limits"]["number_of_aggregates"], number_of_aggregates_limit);
    inipp::extract(ini.sections["limits"]["n_iter_without_contact"], n_iter_without_contact_limit);
    inipp::extract(ini.sections["limits"]["cpu"], cpu_limit);
    inipp::extract(ini.sections["limits"]["physical_time"], physical_time_limit);
    inipp::extract(ini.sections["limits"]["mean_monomere_per_aggregate"], mean_monomere_per_aggregate_limit);
    // numerics
    inipp::extract(ini.sections["numerics"]["n_verlet_divisions"], n_verlet_divisions);
    // output
    inipp::extract(ini.sections["output"]["output_dir"], output_dir);
    inipp::extract(ini.sections["output"]["n_time_per_file"], n_time_per_file);
    // checks
    number_of_aggregates_limit = MAX(number_of_aggregates_limit, 1);
    if (monomeres_initialisation_type == MonomeresInitialisationMode::INVALID_INITIALISATION) {
        exit(ErrorCodes::INPUT_ERROR);
    }
    output_dir = extract_path(fichier_param) / output_dir;
    if (!fs::exists(output_dir)) {
        if (!fs::create_directory(output_dir)) {
            cout << "Error creating directory " << output_dir << endl;
            exit(ErrorCodes::IO_ERROR);
        }
    } else {
        if (!fs::is_directory(output_dir)) {
            cout << "Error not a directory " << output_dir << endl;
            exit(ErrorCodes::INPUT_ERROR);
        }
    }
    // secondary variables
    viscosity = _viscosity_ref
                * (_sutterland_interpolation_constant + _temperature_ref)
                / (_sutterland_interpolation_constant + temperature)
                * pow(temperature / _temperature_ref, 1.5); // Schlichting 1979
    gaz_mean_free_path = _mean_free_path_ref
                         * (_pressure_ref / pressure)
                         * (temperature / _temperature_ref)
                         * (1. + _sutterland_interpolation_constant / _temperature_ref)
                         / (1. + _sutterland_interpolation_constant / temperature); // Willeke 1976
    mean_massic_radius = 0.5 * mean_diameter * 1E-9
                         * exp(1.5 * POW_2(log(dispersion_diameter)));  // Hatch-Choate
    friction_exponnant = 0.689 * (1.
                                  + erf(((gaz_mean_free_path / mean_massic_radius) + 4.454)
                                        / 10.628)); // Yon et al. Journal of Aerosol Science vol.87 2015 p.2837 eq.7
    if (monomeres_initialisation_type == MonomeresInitialisationMode::NORMAL_INITIALISATION) {
        box_lenght = mean_diameter * 1E-9 *
                     pow(static_cast<double>(n_monomeres) * _pi / 6. / volume_fraction
                         * (1. + 3. * POW_2(dispersion_diameter / mean_diameter)),
                         1. / 3.);
    } else if (monomeres_initialisation_type == MonomeresInitialisationMode::LOG_NORMAL_INITIALISATION) {
        box_lenght = mean_diameter * 1E-9 *
                     pow(static_cast<double>(n_monomeres) * _pi / 6. / volume_fraction
                         * exp(9. / 2. * POW_2(log(dispersion_diameter))),
                         1. / 3.);
    } else {
        exit(ErrorCodes::INPUT_ERROR);
    }
    a_surfgrowth = coeff_b * 1E-3;
    print();
    std::ofstream os(output_dir / "params.ini");
    ini.generate(os);
    cpu_start = clock();
}
[[gnu::pure]] bool PhysicalModel::finished(size_t number_of_aggregates, double mean_monomere_per_aggregate) const {
    if (number_of_aggregates <= number_of_aggregates_limit) {
        cout << "We reach the AggMin condition" << endl << endl;
        return true;
    }
    if (n_iter_without_contact_limit > 0
        && n_iter_without_contact >= n_iter_without_contact_limit) {
        cout << "We reach the WaitLimit condition" << endl << endl;
        return true;
    }
    if (cpu_limit > 0) {
        clock_t current_cpu = clock();
        double elapse = double(current_cpu - cpu_start) / CLOCKS_PER_SEC;
        if (elapse >= cpu_limit) {
            cout << "We reach the CPULimit condition" << endl << endl;
            return true;
        }
    }
    if (physical_time_limit > 0
        && time >= physical_time_limit) {
        cout << "We reach the Maximum physical time condition " << time << "/" << physical_time_limit << endl;
        return true;
    }
    if (mean_monomere_per_aggregate_limit > 0
        && mean_monomere_per_aggregate >= mean_monomere_per_aggregate_limit) {
        cout << "We reach the NPP_avg_limit condition " << mean_monomere_per_aggregate << "/"
             << mean_monomere_per_aggregate_limit << endl;
        return true;
    }
    return false;
}
[[gnu::pure]] void PhysicalModel::print() const {
    cout << "Physical parameters:" << endl
         << " Initial Nagg : " << n_monomeres << endl
         << " Box size     : " << box_lenght << endl
         << " Pressure     : " << pressure << endl
         << " Temperature  : " << temperature << endl
         << " diffusivity  : " << viscosity << endl
         << " FV           : " << volume_fraction << endl
         << " density      : " << density << endl
         << " Dpm          : " << mean_diameter << endl
         << " sigmaDpm     : " << dispersion_diameter << endl
         << " Asurfgrowth  : " << a_surfgrowth << endl
         << " xsurfgrowth  : " << x_surfgrowth << endl
         << " coeffB       : " << coeff_b << endl
         << " dfe          : " << dfe << endl
         << " kfe          : " << kfe << endl
         << " lambda       : " << gaz_mean_free_path << endl
         << " rpeqmass     : " << mean_massic_radius << endl
         << " gamma_       : " << friction_exponnant << endl
         << endl
         << "Options for Pysical model: " << endl
         << " Mode : " << monomeres_initialisation_type << endl
         << endl
         << "Ending calcul when:" << endl
         << " - There is " << number_of_aggregates_limit << " aggregats left or less" << endl;
    if (n_iter_without_contact_limit > 0) {
        cout << " - It has been " << n_iter_without_contact_limit << " iterations without collision" << endl;
    }
    if (cpu_limit > 0) {
        cout << " - The simulations is running for more than " << cpu_limit << " seconds" << endl;
    }
    if (physical_time_limit > 0) {
        cout << " - The residence time is larger than " << physical_time_limit << " seconds" << endl;
    }
    if (mean_monomere_per_aggregate_limit > 0) {
        cout << " - The average Npp per aggregate reach " << mean_monomere_per_aggregate_limit << " monomers" << endl;
    }
    cout << endl;
}


//#####################################################################################################################

[[gnu::pure]]  double PhysicalModel::cunningham(double _r) const {
    double A = 1.142;
    double B = 0.558;
    double C = 0.999;
    return 1.0 + A * gaz_mean_free_path / _r
           + B * gaz_mean_free_path / _r * exp(-C * _r / gaz_mean_free_path);
}


//############################# Fonctions pour le calcul du diamètre de mobilité ################################

/*
 Fonction permettant de retrouver le rayon de mobilité en régime transitoire
 On obtient le bon rayon de mobilité lorsque la fonction retourne 0
*/

[[gnu::pure]]  double PhysicalModel::grow(double r, double dt) const {
    return r + a_surfgrowth * pow(r, x_surfgrowth - 2) * dt;
}
[[gnu::pure]] double PhysicalModel::friction_coeff(size_t npp) const {
    //to be expressed in terms of V_agg/V_pp
    double cc = cunningham(mean_massic_radius);
    return (6.0 * _pi * viscosity * mean_massic_radius / cc) * pow(npp, friction_exponnant / dfe);
}
[[gnu::pure]] double PhysicalModel::friction_coeff2(double rgg) const {
    double rmm = rgg * 1.3;
    double cc = cunningham(rmm);
    return (6.0 * _pi * viscosity * rmm / cc);
}
[[gnu::pure]] double PhysicalModel::diffusivity(double f_agg) const {
    return _boltzmann * temperature / f_agg;
}
[[gnu::pure]] double PhysicalModel::relax_time(double masse, double f_agg) const {
    return masse / f_agg;
}
[[gnu::const]] double inverfc(double p) {
    double x;
    double t;
    double pp;
    if (p >= 2.0) {
        return -100.;
    }
    if (p <= 0.0) {
        return 100.;
    }
    pp = (p < 1.0) ? p : 2. - p;
    t = sqrt(-2. * log(pp / 2.));
    x = -0.70711 * ((2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t);
    for (int j = 0; j < 2; j++) {
        double err = erfc(x) - pp;
        x += err / (1.12837916709551257 * exp(-(x * x)) - x * err);
    }
    return (p < 1.0 ? x : -x);
}
[[gnu::const]] double inverf(const double p) {
    return inverfc(1. - p);
}
[[gnu::const]] fs::path extract_path(const string &filename) {
    fs::path path = filename;
    fs::path fullpath = fs::absolute(path);
    if (!fs::exists(fullpath)) {
        cout << "File does not exist\n" << endl;
        exit(1);
    }
    fs::path parentpath = fullpath.parent_path();
    return parentpath; //Cette variable ne retient que le chemin du fichier param
}
[[gnu::const]] tuple<bool, double, double, double> linreg(const vector<double> &x, const vector<double> &y) {
    double sumx = 0.0;                        /* sum of x                      */
    double sumx2 = 0.0;                       /* sum of x**2                   */
    double sumxy = 0.0;                       /* sum of x * y                  */
    double sumy = 0.0;                        /* sum of y                      */
    double sumy2 = 0.0;                       /* sum of y**2                   */

    size_t n = x.size();
    auto N = double(n);
    for (size_t i = 0; i < n; i++) {
        double X(log(x[i]));
        double Y(log(y[i]));
        sumx += X;
        sumx2 += POW_2(X);
        sumxy += X * Y;
        sumy += Y;
        sumy2 += POW_2(Y);
    }
    double denom = (N * sumx2 - POW_2(sumx));
    if (n == 0 || abs(denom) < 1e-9) {
        // singular matrix. can't solve the problem.
        tuple<bool, double, double, double> res{false, 0., 0., 0.};
        return res;
    }
    double a = (N * sumxy - sumx * sumy) / denom;
    double b = (sumy * sumx2 - sumx * sumxy) / denom;

    /* compute correlation coeff     */
    double r = (sumxy - sumx * sumy / N) /
               POW_2((sumx2 - POW_2(sumx) / N) *
                     (sumy2 - POW_2(sumy) / N));
    tuple<bool, double, double, double> res{true, a, b, r};
    return res;
}
}  // namespace MCAC

