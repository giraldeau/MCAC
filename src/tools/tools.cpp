#include "tools/tools.hpp"
#include <vector>


namespace mcac {
void init_random() {
    time_t t;
    time(&t);
    //srand(uint(t));
    srand(0);
}
double random() {
    double v = rand();
    v = v / RAND_MAX;
    return v;
}
MonomeresInitialisationMode resolve_monomeres_initialisation_mode(const std::string& input) {
    const std::map<std::string, MonomeresInitialisationMode> _initialisation_mode_strings{
        {"lognormal", LOG_NORMAL_INITIALISATION},
        {"normal", NORMAL_INITIALISATION},
    };
    auto itr = _initialisation_mode_strings.find(input);
    if (itr != _initialisation_mode_strings.end()) {
        return itr->second;
    }
    return INVALID_INITIALISATION;
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
        double log_x(log(x[i]));
        double log_y(log(y[i]));
        sumx += log_x;
        sumx_2 += POW_2(log_x);
        sumxy += log_x * log_y;
        sumy += log_y;
        sumy_2 += POW_2(log_y);
    }
    double denom = (n * sumx_2 - POW_2(sumx));
    if (static_cast<size_t>(n) == 0 || fabs(denom) < 1e-9) {
        // singular matrix. can't solve the problem.
        return {false, 0., 0., 0.};
    }
    double a = (n * sumxy - sumx * sumy) / denom;
    double b = (sumy * sumx_2 - sumx * sumxy) / denom;

    /* compute correlation coeff     */
    double r = (sumxy - sumx * sumy / n) /
               POW_2((sumx_2 - POW_2(sumx) / n) *
                     (sumy_2 - POW_2(sumy) / n));
    return {true, a, b, r};
}
}  // namespace mcac
