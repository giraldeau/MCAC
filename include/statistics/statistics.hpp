#ifndef STATISTICS_H
#define STATISTICS_H 1

namespace MCAC {
class StatisticStorage;

class AggregatList;

class Aggregate;

// This routine will be called in order to filter what we keep in the stats
struct StatcmpAgg {
    bool operator()(const Aggregate &lhs, const Aggregate &rhs) const;
};

struct Statcmpdouble {
    bool operator()(const double &lhs, const double &rhs) const;
};

class StatisicsData {
    friend class Aggregate;

    friend class StatisticStorage;

public:

    // Here you can declare all the variables you need for statistical analysis

    double Dp;
    double DgOverDp;

    //double Dp3;

    //double DmOverDp;
    //double DgeoOverDp;
    //double SurfaceOverVolume;

    //double cov;
    //double Npe;
    //double Tv;
    //double Ts;

    //std::size_t Nc;

public:
    // This routine will be called for each update of the aggregate
    void partialStatistics();
    // This routine will be called after each aggregation
    void fullStatistics();
};
}  // namespace MCAC

#endif // STATISTICS_H
