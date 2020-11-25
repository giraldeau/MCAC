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
#ifndef INCLUDE_CONSTANTS_HPP
#define INCLUDE_CONSTANTS_HPP 1
#include <cmath>
#include <map>
#include <string>


#ifdef WITH_HDF5
#define XMF_OUTPUT shared_ptr<XdmfInformation>
#else
#define XMF_OUTPUT void
#endif
namespace mcac {
enum class WriterStatus {
    IDLE, APPENDING, WRITING
};
enum SpheresFields {
    SPHERE_X,
    SPHERE_Y,
    SPHERE_Z,
    SPHERE_R,
    SPHERE_VOLUME,
    SPHERE_SURFACE,
    SPHERE_RX,
    SPHERE_RY,
    SPHERE_RZ,
    SPHERE_NFIELDS
};
enum AggregatesFields {
    AGGREGAT_RG,
    AGGREGAT_F_AGG,
    AGGREGAT_LPM,
    AGGREGAT_TIME_STEP,
    AGGREGAT_RMAX,
    AGGREGAT_VOLUME,
    AGGREGAT_SURFACE,
    AGGREGAT_X,
    AGGREGAT_Y,
    AGGREGAT_Z,
    AGGREGAT_RX,
    AGGREGAT_RY,
    AGGREGAT_RZ,
    AGGREGAT_TIME,
    AGGREGAT_DP,
    AGGREGAT_DG_OVER_DP,
    AGGREGAT_OVERLAPPING,
    AGGREGAT_COORDINATION_NUMBER,
    AGGREGAT_D_M,
    AGGREGAT_NFIELDS
};
enum ErrorCodes {
    NO_ERROR,
    UNKNOWN_ERROR,
    IO_ERROR,
    VERLET_ERROR,
    INPUT_ERROR,
    ABANDON_ERROR,
    TOO_DENSE_ERROR,
    SBL_ERROR,
    VOL_SURF_ERROR,
    MERGE_ERROR,
};
enum PickMethods {
    PICK_RANDOM,
    PICK_LAST,
    INVALID_PICK_METHOD
};
static const std::array<const std::string, 2> PICK_METHODS_STR{{"random", "last"}};
enum MonomeresInitialisationMode {
    LOG_NORMAL_INITIALISATION,
    NORMAL_INITIALISATION,
    INVALID_INITIALISATION
};
static const std::array<const std::string, 2> MONOMERES_INITIALISATION_MODE_STR{{"lognormal", "normal"}};
enum VolSurfMethods {
    SPHERICAL_CAPS,
    EXACT_SBL,
    EXACT_ARVO,
    ALPHAS,
    NONE,
    INVALID_VOLSURF_METHOD
};
static const std::array<const std::string, 5> VOLSURF_METHODS_STR{{"caps", "sbl", "arvo", "alphas", "none"}};

const double _pi = atan(1.0) * 4;
const double _volume_factor = 4 * _pi / 3;
const double _surface_factor = 4 * _pi;
const double _boltzmann = 1.38066E-23;
const double _fluid_mean_free_path_ref = 66.5E-9;
const double _temperature_ref = 293.15;
const double _sutterland_interpolation_constant = 110;
const double _pressure_ref = 101300;
const double _viscosity_ref = 18.203E-6;
}  //namespace mcac
#endif //INCLUDE_CONSTANTS_HPP
