#ifndef INCLUDE_SBL_VOLUME_SURFACE_HPP
#define INCLUDE_SBL_VOLUME_SURFACE_HPP 1

#ifdef WITH_SBL

#include <Spherelist.hpp>

#include <vector>

using namespace std;
namespace mcac{

pair<vector<double>, vector<double>> compute_volume_surface(const SphereList& spherelist);

}

#endif // WITH_SBL
#endif //INCLUDE_SBL_VOLUME_SURFACE_HPP
