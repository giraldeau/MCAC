#ifndef INCLUDE_TOOLS_TOOLS_HPP
#define INCLUDE_TOOLS_TOOLS_HPP 1
#include "constants.hpp"
#include <array>
#include <vector>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
#define POW_3(a) ((a)*(a)*(a))
namespace mcac {
template<class T>
std::array<T, 3> operator+(const std::array<T, 3> &a_1, const std::array<T, 3> &a_2) {
    return {a_1[0] + a_2[0], a_1[1] + a_2[1], a_1[2] + a_2[2]};
}
template<class T>
std::array<T, 3> operator-(const std::array<T, 3> &a_1, const std::array<T, 3> &a_2) {
    return {a_1[0] - a_2[0], a_1[1] - a_2[1], a_1[2] - a_2[2]};
}
template<class T>
std::array<T, 3> operator*(const std::array<T, 3> &a_1, const std::array<T, 3> &a_2) {
    return {a_1[0] * a_2[0], a_1[1] * a_2[1], a_1[2] * a_2[2]};
}
template<class T>
std::array<T, 3> operator*(double f, const std::array<T, 3> &a) {
    return {f * a[0], f * a[1], f * a[2]};
}
template<class T>
std::array<T, 3> operator*(const std::array<T, 3> &a, double f) {
    return {a[0] * f, a[1] * f, a[2] * f};
}
template<class T>
T sum(std::array<T, 3> a) {
    return a[0] + a[1] + a[2];
}
template<class T>
std::array<T, 3> &operator+=(std::array<T, 3> &a_1, const std::array<T, 3> &a_2) {
    a_1[0] += a_2[0];
    a_1[1] += a_2[1];
    a_1[2] += a_2[2];
    return a_1;
}
template<class T>
std::array<T, 3> &operator/=(std::array<T, 3> &a, const T &rhs) {
    a[0] /= rhs;
    a[1] /= rhs;
    a[2] /= rhs;
    return a;
}
void init_random();
double random();
MonomeresInitialisationMode resolve_monomeres_initialisation_mode(const std::string& input);
[[gnu::const]] std::tuple<bool, double, double, double> linreg(const std::vector<double> &x,
                                                               const std::vector<double> &y);
}  //namespace mcac
#endif //INCLUDE_TOOLS_TOOLS_HPP
