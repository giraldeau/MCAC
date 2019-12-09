#ifndef INCLUDE_CONSTANTS_HPP
#define INCLUDE_CONSTANTS_HPP 1
#include <cmath>
#include <map>

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
    AGGREGAT_NFIELDS
};
enum ErrorCodes {
    NO_ERROR,
    UNKNOWN_ERROR,
    IO_ERROR,
    VERLET_ERROR,
    INPUT_ERROR
};
enum MonomeresInitialisationMode {
    LOG_NORMAL_INITIALISATION,
    NORMAL_INITIALISATION,
    INVALID_INITIALISATION,
};

const double _pi = atan(1.0) * 4;
const double _volume_factor = 4 * _pi / 3;
const double _surface_factor = 4 * _pi;
const double _boltzmann = 1.38066E-23;
const double _mean_free_path_ref = 66.5E-9;
const double _temperature_ref = 293.15;
const double _sutterland_interpolation_constant = 110;
const double _pressure_ref = 101300;
const double _viscosity_ref = 18.203E-6;
}  //namespace mcac
#endif //INCLUDE_CONSTANTS_HPP
