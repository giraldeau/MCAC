#ifndef INCLUDE_TOOLS_TOOLS_HPP_
#define INCLUDE_TOOLS_TOOLS_HPP_ 1
#include "constants.hpp"
#include <array>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW_2(a) ((a)*(a))
#define POW_3(a) ((a)*(a)*(a))
namespace MCAC {
template<class T>
std::array<T, 3> operator+(const std::array<T, 3> &a1, const std::array<T, 3> &a2) {
    return {a1[0] + a2[0], a1[1] + a2[1], a1[2] + a2[2]};
}
template<class T>
std::array<T, 3> operator-(const std::array<T, 3> &a1, const std::array<T, 3> &a2) {
    return {a1[0] - a2[0], a1[1] - a2[1], a1[2] - a2[2]};
}
template<class T>
std::array<T, 3> operator*(const std::array<T, 3> &a1, const std::array<T, 3> &a2) {
    return {a1[0] * a2[0], a1[1] * a2[1], a1[2] * a2[2]};
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
std::array<T, 3> &operator+=(std::array<T, 3> &a1, const std::array<T, 3> &a2) {
    a1[0] += a2[0];
    a1[1] += a2[1];
    a1[2] += a2[2];
    return a1;
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
MonomeresInitialisationMode resolve_monomeres_initialisation_mode(std::string input);
}  //namespace MCAC
#endif //INCLUDE_TOOLS_TOOLS_HPP_
