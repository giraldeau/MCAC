#ifndef STATISTICS_H
#define STATISTICS_H 1

#include "io/threaded_io.hpp"
#include "physical_model.hpp"
#include <set>
#include <vector>
#include <tuple>


namespace MCAC{

class StatisticStorage;
class ListAggregat;
class Aggregate;

// This routine will be called in order to filter what we keep in the stats
struct StatcmpAgg{
  bool operator() (const Aggregate& lhs, const Aggregate& rhs) const;
};

struct Statcmpdouble{
  bool operator() (const double& lhs, const double& rhs) const;
};

class StatisicsData
{
    friend class Aggregate;
    friend class StatisticStorage;
    friend struct StatcmpAgg;

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



class StatisticStorage
{
private:
    PhysicalModel* physicalmodel;
    ListAggregat* SavedAggregates;
    std::vector< std::set<double,Statcmpdouble> > FractalLaw;
    std::vector< double > times;

    ThreadedIO* WriterAgg;
    ThreadedIO* WriterSph;


public:

    void Init();

    void Analyze(const ListAggregat& current);
    void print() const;
    bool InsertIfNew(const Aggregate& Agg);

    void Save();
    void Save(bool finish);

    std::vector<double> FormatTimeData() const;
    std::tuple<bool,double,double,double> getCompleteFractalLaw() const;

public:
    /** Default constructor */
    explicit StatisticStorage(PhysicalModel& _physicalmodel);

    /** Destructor */
    ~StatisticStorage() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

    /** Copy assignment operator */
     StatisticStorage& operator= (const StatisticStorage& other);

    /** Move assignment operator */
     StatisticStorage& operator= (StatisticStorage&& other) noexcept;
};

}  // namespace MCAC

#endif // STATISTICS_H