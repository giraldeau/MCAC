#ifndef SBLVOLUMESURFACE_H
#define SBLVOLUMESURFACE_H

#ifdef WITH_SBL

#include <Spherelist.hpp>

#include <vector>
#include <array>
#include <utility>

using namespace std;
namespace DLCA{

pair<vector<double>, vector<double>> compute_volume_surface(const ListSphere& spherelist);

}

#endif
#endif // SBLVOLUMESURFACE_H
