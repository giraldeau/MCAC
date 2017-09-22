#ifndef STATISTICS_H
#define STATISTICS_H

#include "aggregat.hpp"
#include "aggregatList.hpp"

#include <vector>

namespace DLCA{

class AnalizedAggregate;
class StatisticStorage;

class AnalizedAggregate: public Aggregate
{
    friend class StatisticStorage;

private:
    size_t Nc;

    double Dp3;

    double DmOverDp;
    double DgeoOverDp;
    double SurfaceOverVolume;

    double cov;
    double Npe;
    double Tv;
    double Ts;

public:
    explicit AnalizedAggregate(const Aggregate& source);
    void Analyze();

};

class StatisticStorage
{
private:

    std::vector<bool> targetNp;

    std::vector<AnalizedAggregate> SavedAggregates;
    PhysicalModel* physicalmodel;

public:

    void Init();

    void Analyze(const ListAggregat& current);
    void append(const AnalizedAggregate& Agg);
    void print() const;
    bool InsertIfNew(const Aggregate& Agg);

public:
    /** Default constructor */
    explicit StatisticStorage(PhysicalModel& _physicalmodel);

    /** Copy constructor */
    StatisticStorage(const StatisticStorage&);

    /** Move constructor */
    StatisticStorage (StatisticStorage&&) noexcept; /* noexcept needed to enable optimizations in containers */

    /** Destructor */
    ~StatisticStorage() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
    StatisticStorage& operator= (const StatisticStorage& other);

    /** Move assignment operator */
    StatisticStorage& operator= (StatisticStorage&& other) noexcept;
};
}  // namespace DLCA

#endif // STATISTICS_H
