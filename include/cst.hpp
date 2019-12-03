#ifndef INCLUDE_CST_HPP_
#define INCLUDE_CST_HPP_ 1
#include <cmath>


#ifdef WITH_HDF5
#define XMF_OUTPUT shared_ptr<XdmfInformation>
#else
#define XMF_OUTPUT void
#endif
namespace MCAC {
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
    AGGREGAT_NFIELDS
};
enum ErrorCodes {
    NO_ERROR,
    IO_ERROR,
    VERLET_ERROR,
};
const double _pi = atan(1.0) * 4;
const double _facvol = 4 * _pi / 3;
const double _facsurf = 4 * _pi;
}  //namespace MCAC
#endif //INCLUDE_CST_HPP_
