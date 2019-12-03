#ifndef INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP_
#define INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP_ 1
#include "aggregat.hpp"
#include <array>


namespace MCAC {
double distance(const Aggregate &agregate_1, const Aggregate &agregate_2, std::array<double, 3> direction) noexcept;
bool contact(const Aggregate &aggregate_1, const Aggregate &aggregate_2) noexcept;
}  // namespace MCAC


#endif //INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP_
