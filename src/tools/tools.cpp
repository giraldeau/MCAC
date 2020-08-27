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
#include "tools/tools.hpp"
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>


namespace mcac {
void init_random() {
#ifdef DEBUG_MCAC
    srand(0);
#else
    time_t t;
    time(&t);
    srand(uint(t));
#endif
}
double random() {
    double v = rand();
    v = v / RAND_MAX;
    return v;
}
MonomeresInitialisationMode resolve_monomeres_initialisation_mode(const std::string &input) {
    auto itr = std::find(MONOMERES_INITIALISATION_MODE_STR.begin(),
                         MONOMERES_INITIALISATION_MODE_STR.end(),
                         input);
    if (itr != MONOMERES_INITIALISATION_MODE_STR.end()) {
        return MonomeresInitialisationMode(std::distance(MONOMERES_INITIALISATION_MODE_STR.begin(), itr));
    }
    return INVALID_INITIALISATION;
}
std::string resolve_monomeres_initialisation_mode(MonomeresInitialisationMode mode) {
    return MONOMERES_INITIALISATION_MODE_STR[mode];
}
PickMethods resolve_pick_method(const std::string &input) {
    auto itr = std::find(PICK_METHODS_STR.begin(),
                         PICK_METHODS_STR.end(),
                         input);
    if (itr != PICK_METHODS_STR.end()) {
        return PickMethods(std::distance(PICK_METHODS_STR.begin(), itr));
    }
    return INVALID_PICK_METHOD;
}
std::string resolve_pick_method(PickMethods method) {
    return PICK_METHODS_STR[method];
}
VolSurfMethods resolve_surfvol_method(const std::string &input) {
    auto itr = std::find(VOLSURF_METHODS_STR.begin(),
                         VOLSURF_METHODS_STR.end(),
                         input);
    if (itr != VOLSURF_METHODS_STR.end()) {
        return VolSurfMethods(std::distance(VOLSURF_METHODS_STR.begin(), itr));
    }
    return INVALID_VOLSURF_METHOD;
}
std::string resolve_surfvol_method(VolSurfMethods method) {
    return VOLSURF_METHODS_STR[method];
}
[[gnu::const]] std::tuple<bool, double, double, double> linreg(const std::vector<double> &x,
                                                               const std::vector<double> &y) {
    double sumx = 0.0;                        /* sum of x                      */
    double sumx_2 = 0.0;                       /* sum of x**2                   */
    double sumxy = 0.0;                       /* sum of x * y                  */
    double sumy = 0.0;                        /* sum of y                      */
    double sumy_2 = 0.0;                       /* sum of y**2                   */

    auto n = static_cast<double>(x.size());
    for (size_t i = 0; i < static_cast<size_t>(n); i++) {
        double log_x(std::log(x[i]));
        double log_y(std::log(y[i]));
        sumx += log_x;
        sumx_2 += std::pow(log_x, 2);
        sumxy += log_x * log_y;
        sumy += log_y;
        sumy_2 += std::pow(log_y, 2);
    }
    double denom = (n * sumx_2 - std::pow(sumx, 2));
    if (static_cast<size_t>(n) == 0 || std::abs(denom) < 1e-9) {
        // singular matrix. can't solve the problem.
        return {false, 0., 0., 0.};
    }
    double a = (n * sumxy - sumx * sumy) / denom;
    double b = (sumy * sumx_2 - sumx * sumxy) / denom;

    /* compute correlation coeff     */
    double r = (sumxy - sumx * sumy / n) /
               std::pow((sumx_2 - std::pow(sumx, 2) / n) *
                        (sumy_2 - std::pow(sumy, 2) / n), 2);
    return {true, a, b, r};
}
}  // namespace mcac
