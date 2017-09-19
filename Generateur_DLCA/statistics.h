#ifndef STATISTICS_H
#define STATISTICS_H

#include "aggregat.h"
#include "aggregatList.h"

#include <vector>

namespace DLCA{

class AnalizeAggregate;
class Statistics;

class AnalizeAggregate: public Aggregate
{
    friend class Statistics;

private:
    int Nc;

    double Dp;
    double Dp3;

    double DgOverDp;
    double DmOverDp;
    double DgeoOverDp;
    double SurfaceOverVolume;

    double cov;
    double Npe;
    double Tv;
    double Ts;

public:
    AnalizeAggregate(const Aggregate& source);
};

class Statistics
{
    std::vector<bool> target;

private:
    std::vector<AnalizeAggregate> SavedAggregates;
    PhysicalModel* physicalmodel;

public:

    Statistics(PhysicalModel& _physicalmodel);
    void Init(void);

    void Analyze(const ListAggregat& current);
    void append(const AnalizeAggregate& Agg);
    void print(void) const;
    bool keep(const Aggregate& Agg) const;


};
}
#endif // STATISTICS_H
