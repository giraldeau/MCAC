#ifndef INCLUDE_CST_HPP_
#define INCLUDE_CST_HPP_ 1
#include <cmath>

#define AGG_N_FIELD 15

#ifdef WITH_HDF5
#define XMF_OUTPUT shared_ptr<XdmfInformation>

#else
#define XMF_OUTPUT void
#endif



namespace MCAC {

enum class WriterStatus { IDLE, APPENDING, WRITING };
enum SpheresFields { X, Y, Z, R, VOLUME, SURFACE, RX, RY, RZ, NFIELD };

const double _pi = atan(1.0) * 4;
const double _facvol = 4 * _pi / 3;
const double _facsurf = 4 * _pi;

}  //namespace MCAC
#endif //INCLUDE_CST_HPP_
