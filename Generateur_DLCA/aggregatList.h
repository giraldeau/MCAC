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


using namespace std;

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;

class ListAggregat : public storage_list<16,Aggregate>
{
    friend class Aggregate;

    private:
        PhysicalModel* physicalmodel;
        double maxradius;
        double maxtime_step;

    public:
        ListSphere spheres;
        Verlet verlet;
        void Init(PhysicalModel&, const int size);

        vector<int> PotentialContacts(int id,array<double,4> Vectdir, vector<int> SearchSpace);
        vector<int> GetSearchSpace(int source, array<double,4> Vectdir);
        int DistanceToNextContact(const int source, const array<double,4> Vectdir, double &distmin);

        double GetMaxTimeStep();

        /* Storage specific */
    private:
        void setpointers();

    public:
        /** Default constructor in local storage */
        ListAggregat(void);
        ListAggregat(PhysicalModel& _physicalmodel, const int size);

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
