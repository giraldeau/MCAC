#ifndef INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP
#define INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP 1
#include "aggregat.hpp"
#include <array>


namespace mcac {
double distance(const Aggregate &aggregate_1, const Aggregate &aggregate_2, std::array<double, 3> direction) noexcept;
bool contact(const Aggregate &aggregate_1, const Aggregate &aggregate_2) noexcept;
}  // namespace mcac


#endif //INCLUDE_AGGREGATS_AGGREGAT_DISTANCE_HPP
