#include "tools/tools.hpp"


namespace MCAC {
void InitRandom() {
    time_t t;
    time(&t);
    //srand(uint(t));
    srand(0);
}
double Random() {
    double v = rand();
    v = v / RAND_MAX;
    return v;
}
}  // namespace MCAC
