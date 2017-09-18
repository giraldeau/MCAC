#ifndef AGGREGATLIST_H
#define AGGREGATLIST_H

#include "Sphere.h"
#include "Spherelist.h"
#include "aggregat.h"
#include "aggregatList.h"
#include "physical_model.h"
#include "storage.h"
#include "storagelist.h"
#include "verlet.h"
#include "IO.h"


class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;

class ListAggregat : public storage_list<13,Aggregate>
{
    friend class Aggregate;

    private:
        PhysicalModel* physicalmodel;

        double maxradius;
        std::vector<int> indexSortedTimeSteps;
        std::vector<double> CumulativeTimeSteps;

        std::vector<double>::iterator ptr_deb;
        std::vector<double>::iterator ptr_fin;

        ThreadedIO* Writer;

    public:
        ListSphere spheres;
        Verlet verlet;
        void Init(PhysicalModel&, const int _size);

        std::vector<int> PotentialCollision(int id,std::array<double,3> Vectdir, std::vector<int> SearchSpace) const;
        std::vector<int> GetSearchSpace(int source, std::array<double,3> Vectdir) const;
        int DistanceToNextContact(const int source, const std::array<double,3> Vectdir, double &distmin) const;

        double GetMaxTimeStep() const;

        int Merge(const int first, const int second);
        void SortTimeSteps(double factor);
        int RandomPick(double &deltatemps, const double random);



        void save(void) const;
        void save(const bool) const;

        std::vector<double> FormatPositionData() const;
        std::vector<double> FormatRgData() const;
        std::vector<int>    FormatNpData() const;
        std::vector<double> FormatDmData() const;
        std::vector<double> FormatlpmData() const;
        std::vector<double> FormatdeltatData() const;
        std::vector<double> FormatRmaxData() const;
        std::vector<double> FormatVoumeData() const;
        std::vector<double> FormatSurfaceData() const;
        std::vector<int>    FormatLabelData() const;


        /* Storage specific */
    private:
        void setpointers();

    public:
        /** Default constructor in local storage */
        ListAggregat(void);
        ListAggregat(PhysicalModel& _physicalmodel, const int _size);

        /** Copy constructor */
        ListAggregat(const ListAggregat& other);

        /** Move constructor */
        ListAggregat (ListAggregat&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~ListAggregat(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        ListAggregat& operator= (const ListAggregat& other);

        /** Move assignment operator */
        ListAggregat& operator= (ListAggregat&& other) noexcept;
};



#endif // AGGREGATLIST_H
