#include "aggregats/aggregat_list.hpp"
#include "tools/tools.hpp"


using namespace std;
namespace MCAC {
tuple<bool, double, double, double> AggregatList::get_instantaneous_fractal_law() const {
    vector<double> sizes;
    vector<double> _dg_over_dps;
    sizes.reserve(size());
    _dg_over_dps.reserve(size());
    for (const Aggregate *Agg : list) {
        sizes.push_back(double(Agg->size()));
        _dg_over_dps.push_back(*Agg->dg_over_dp);
    }
    return linreg(_dg_over_dps, sizes);
}
}// namespace MCAC

