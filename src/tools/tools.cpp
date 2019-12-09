#include "tools/tools.hpp"
#include <vector>


namespace MCAC {
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
MonomeresInitialisationMode resolve_monomeres_initialisation_mode(std::string input) {
    const std::map<std::string, MonomeresInitialisationMode> InitialisationModeStrings{
        {"lognormal", LOG_NORMAL_INITIALISATION},
        {"normal", NORMAL_INITIALISATION},
    };
    auto itr = InitialisationModeStrings.find(input);
    if (itr != InitialisationModeStrings.end()) {
        return itr->second;
    }
    return INVALID_INITIALISATION;
}
[[gnu::const]] std::tuple<bool, double, double, double> linreg(const std::vector<double> &x,
                                                               const std::vector<double> &y) {
    double sumx = 0.0;                        /* sum of x                      */
    double sumx2 = 0.0;                       /* sum of x**2                   */
    double sumxy = 0.0;                       /* sum of x * y                  */
    double sumy = 0.0;                        /* sum of y                      */
    double sumy2 = 0.0;                       /* sum of y**2                   */

    size_t n = x.size();
    auto N = double(n);
    for (size_t i = 0; i < n; i++) {
        double X(log(x[i]));
        double Y(log(y[i]));
        sumx += X;
        sumx2 += POW_2(X);
        sumxy += X * Y;
        sumy += Y;
        sumy2 += POW_2(Y);
    }
    double denom = (N * sumx2 - POW_2(sumx));
    if (n == 0 || fabs(denom) < 1e-9) {
        // singular matrix. can't solve the problem.
        return {false, 0., 0., 0.};
    }
    double a = (N * sumxy - sumx * sumy) / denom;
    double b = (sumy * sumx2 - sumx * sumxy) / denom;

    /* compute correlation coeff     */
    double r = (sumxy - sumx * sumy / N) /
               POW_2((sumx2 - POW_2(sumx) / N) *
                     (sumy2 - POW_2(sumy) / N));
    return {true, a, b, r};
}
}  // namespace MCAC
