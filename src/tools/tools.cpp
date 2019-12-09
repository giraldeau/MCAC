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
MonomeresInitialisationMode resolveMonomeresInitialisationMode(std::string input) {
    const std::map<std::string, MonomeresInitialisationMode> InitialisationModeStrings {
        { "lognormal", LOG_NORMAL_INITIALISATION },
        { "normal", NORMAL_INITIALISATION },
    };

    auto itr = InitialisationModeStrings.find(input);
    if (itr != InitialisationModeStrings.end()) {
        return itr->second;
    }
    return INVALID_INITIALISATION;
}
}  // namespace MCAC
