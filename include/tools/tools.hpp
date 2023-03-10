/*
 * MCAC
 * Copyright (C) 2020 CORIA
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef INCLUDE_TOOLS_TOOLS_HPP
#define INCLUDE_TOOLS_TOOLS_HPP 1
#include "constants.hpp"
#include <array>
#include <vector>


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
unsigned long mix(unsigned long a, unsigned long b, unsigned long c);
void init_random(int random_seed);
[[gnu::const]] double inverfc(double p);
[[gnu::const]] double inverf(double p);
double random();
double random_normal(const double mean, const double sigma);
std::array<double, 3> random_direction();
MonomeresInitialisationMode resolve_monomeres_initialisation_mode(const std::string &input);
std::string resolve_monomeres_initialisation_mode(MonomeresInitialisationMode mode);
PickMethods resolve_pick_method(const std::string &input);
std::string resolve_pick_method(PickMethods method);
VolSurfMethods resolve_surfvol_method(const std::string &input);
std::string resolve_surfvol_method(VolSurfMethods method);
[[gnu::const]] std::tuple<bool, double, double, double> linreg(const std::vector<double> &x,
                                                               const std::vector<double> &y);
[[gnu::pure]] double interpolate_2d(const double& f_x1_y1,
                                    const double& f_x1_y2,
                                    const double& f_x2_y1,
                                    const double& f_x2_y2,
                                    const double& dx_over_Dx,
                                    const double& dy_over_Dy);
}  //namespace mcac
#endif //INCLUDE_TOOLS_TOOLS_HPP
