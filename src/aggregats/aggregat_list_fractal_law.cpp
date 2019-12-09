#include "aggregats/aggregat_list.hpp"
#include "tools/tools.hpp"


namespace mcac {
std::tuple<bool, double, double, double> AggregatList::get_instantaneous_fractal_law() const {
    std::vector<double> sizes;
    std::vector<double> dg_over_dps;
    sizes.reserve(size());
    dg_over_dps.reserve(size());
    for (const Aggregate *agg : list) {
        sizes.push_back(double(agg->size()));
        dg_over_dps.push_back(*agg->dg_over_dp);
    }
    return linreg(dg_over_dps, sizes);
}
}// namespace mcac

