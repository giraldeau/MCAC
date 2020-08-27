//********************************************************************************
// Wrapper for ARVO library
// Jose MORAN, CORIA Laboratory, France, 2020
// Based on: Computer Physics Communications 183 (2012) 2494–2497
//********************************************************************************
#include "arvo/arvo_mcac_call.hpp"
#include "spheres/sphere.hpp"
#include "tools/tools.hpp"
#include <iostream>

/********************************************************************************
* Main program
********************************************************************************/
namespace mcac {
std::pair<std::vector<double>, std::vector<double>>
arvo_call(const SphereList &spherelist,
          const std::vector<std::unordered_map<size_t, double>> &distances) {
    std::vector<double> one_column_data(spherelist.size() * 4);
    std::vector<double> volumes(spherelist.size());
    std::vector<double> surfaces(spherelist.size());

    // convert data for ARVO
    // one column array (x,y,z,r,...)
    size_t i(0);
    for (const auto &sphere :spherelist) {
        std::array<double, 3> pos(sphere->get_position());
        one_column_data[i * 4] = pos[0] * 1e+09; // convert to nm for better accuracy
        one_column_data[i * 4 + 1] = pos[1] * 1e+09;
        one_column_data[i * 4 + 2] = pos[2] * 1e+09;
        one_column_data[i * 4 + 3] = (sphere->get_radius()) * 1e+09;
        i += 1;
    }
    // a parallelized version is invoked only for large aggregates
    main_arvo(one_column_data.data(), volumes.data(), surfaces.data(), spherelist.size(), distances);
    return std::make_pair(volumes, surfaces);
}
}  // namespace mcac
