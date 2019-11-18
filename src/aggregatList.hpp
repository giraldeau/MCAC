#ifndef AGGREGATLIST_H
#define AGGREGATLIST_H 1

#include "Sphere.hpp"
#include "IO.hpp"
#include "Spherelist.hpp"
#include "aggregat.hpp"
#include "physical_model.hpp"
#include "storage.hpp"
#include "storagelist.hpp"
#include "verlet.hpp"

namespace DLCA{

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;
class StatisticStorage;


class ListAggregat : public storage_list<15,Aggregate>
{
    friend class Aggregate;
    friend class StatisticStorage;

    private:
        PhysicalModel* physicalmodel;

        double maxradius;

        std::vector<size_t> indexSortedTimeSteps;
        std::vector<double> CumulativeTimeSteps;

        std::vector<double>::iterator ptr_deb;
        std::vector<double>::iterator ptr_fin;

        ThreadedIO* Writer;
        size_t lastSaved;


    public:
        ListSphere spheres;
        Verlet verlet;
        void Init(PhysicalModel&, size_t _size);
        void Duplication();
        size_t GetAvg_npp();

        std::vector< std::pair<size_t,double> > SortSearchSpace(size_t MovingAgg,std::array<double,3> Vectdir, std::vector<size_t> SearchSpace) const;
        std::vector<size_t> GetSearchSpace(size_t source, std::array<double,3> Vectdir) const;
        int DistanceToNextContact(size_t source, std::array<double,3> Vectdir, double &distmin) const;
        bool TestFreeSpace(const std::array<double,3> pos, const double diameter) const;

        double GetMaxTimeStep() const;
        double GetTimeStep(const double max) const;

        size_t Merge(size_t first, size_t second);
        void SortTimeSteps(double factor);
        size_t RandomPick(double random);
        size_t PickLast();

        std::tuple<bool,double,double,double> getInstantaneousFractalLaw() const;

        void save();
        void save(bool);

        auto GetData() const;

        std::vector<double> FormatPositionData() const;
        std::vector<double> FormatRgData() const;
        std::vector<int>    FormatNpData() const;
        std::vector<double> FormatFaggData() const;
        std::vector<double> FormatlpmData() const;
        std::vector<double> FormatdeltatData() const;
        std::vector<double> FormatRmaxData() const;
        std::vector<double> FormatVolumeData() const;
        std::vector<double> FormatSurfaceData() const;
        std::vector<int>    FormatLabelData() const;
        std::vector<int>    FormatIndexData() const;


        /* Storage specific */
    private:
        void setpointers();

    public:

        Aggregate* add(const Aggregate& newAgg);

        /** Default constructor in local storage */
        ListAggregat();
        explicit ListAggregat(PhysicalModel& _physicalmodel, size_t _size);

        /** Copy constructor */
        ListAggregat(const ListAggregat& other);

        /** Move constructor */
        ListAggregat (ListAggregat&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~ListAggregat() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        ListAggregat& operator= (const ListAggregat& other);

        /** Move assignment operator */
        ListAggregat& operator= (ListAggregat&& other) noexcept;
};

}// namespace DLCA


#endif // AGGREGATLIST_H
