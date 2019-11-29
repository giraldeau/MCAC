#ifndef AGGREGATLIST_H
#define AGGREGATLIST_H 1

#include "spheres/sphere.hpp"
#include "spheres/Spherelist.hpp"
#include "aggregat.hpp"
#include "physical_model.hpp"
#include "storage.hpp"
#include "storagelist.hpp"
#include "verlet.hpp"

namespace MCAC{

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

        std::vector< std::pair<size_t,double> > SortSearchSpace(size_t MovingAgg,
            std::array<double,3> Vectdir,
            const std::vector<size_t>& SearchSpace) const;
        std::vector<size_t> GetSearchSpace(size_t source, std::array<double,3> Vectdir) const;
        int DistanceToNextContact(size_t source, std::array<double,3> Vectdir, double &distmin) const;
        bool TestFreeSpace(std::array<double,3> pos, double diameter) const;

        double GetMaxTimeStep() const;
        double GetTimeStep(double max) const;

        size_t Merge(size_t first, size_t second);
        void SortTimeSteps(double factor);
        size_t RandomPick(double random);
        size_t PickLast();

        std::tuple<bool,double,double,double> getInstantaneousFractalLaw() const;

        void save(){save(false);};
        void save(bool);

        auto get_data() const;

        std::vector<double> Format_Position() const;
        std::vector<double> Format_rg() const;
        std::vector<int>    Format_Np() const;
        std::vector<double> Format_f_agg() const;
        std::vector<double> Format_lpm() const;
        std::vector<double> Format_time_step() const;
        std::vector<double> Format_rmax() const;
        std::vector<double> Format_volAgregat() const;
        std::vector<double> Format_surfAgregat() const;
        std::vector<int>    Format_Label() const;


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

}// namespace MCAC


#endif // AGGREGATLIST_H
