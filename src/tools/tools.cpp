#include "tools/tools.hpp"


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
