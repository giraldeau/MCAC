#ifndef AGGREGAT_H
#define AGGREGAT_H

#include "Sphere.h"
#include "Spherelist.h"
#include "aggregat.h"
#include "aggregatList.h"
#include "storage.h"
#include "storagelist.h"
#include <list>

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;
class Verlet;

class Aggregate : public storage_elem<13,ListAggregat>
{

    friend class ListAggregat;
    friend class Sphere;

    private:
        PhysicalModel* physicalmodel;

        ListSphere myspheres;

        Verlet* verlet;
        std::array<int, 3> IndexVerlet;

        std::vector<std::vector <double > > distances;
        std::vector<double> volumes;
        std::vector<double> surfaces;

        double *rg;  //Gyration Radius
        double *dm;  //Mobility Diameter
        double *lpm; //Mean Free Path
        double *time_step;
        double *rmax;                   //Radius of the sphere containing Agg
        double *volAgregat;             //Etimation of the Aggregate's volume
        double *surfAgregat;            //Estimation of the sufrace of the aggregate

        double *x,*y,*z; // position of the gravity center
        double *rx,*ry,*rz; // position of the gravity center

        int Np;     //Number of spheres

        bool InVerlet;

    public:

        void Init(void);
        void Init(PhysicalModel& ,Verlet&,const std::array<double, 3> position ,const int _label, ListSphere&,double r);

        double GetLpm() const noexcept;
        double GetVolAgregat() const noexcept;
        int GetLabel() const noexcept;


        const std::array<double, 3> GetPosition(void) const noexcept;
        std::array<int, 3> GetVerletIndex() noexcept;


        void SetPosition(const std::array<double, 3> position) noexcept;
        void SetPosition(const double x,const double y,const double z) noexcept;
        void Translate(const std::array<double, 3> vector) noexcept;


        void Update();
        void Volume();
        void MassCenter();
        void CalcRadius();
        void RayonGiration(void);
        bool Contact(Aggregate&) const noexcept;
        double Distance(Aggregate&, std::array<double,3> Vectdir) const;
        void Merge(Aggregate&);
        void DecreaseLabel(void) noexcept;

        void UpdateDistances(void) noexcept;

        void check(void);

        /* Storage specific */
    private:
        void setpointers(void);

    public:
        /** Default constructor in local storage */
        Aggregate(void);
        Aggregate(PhysicalModel&);

        /** Constructor in local storage with initialization */
        Aggregate(PhysicalModel&, const double x, const double y, const double z, const double r);
        Aggregate(PhysicalModel&, const std::array<double, 3> position, const double r);
        Aggregate(const Sphere&);


        /** Constructor with external storage */
        Aggregate(ListSphere& Storage, const int id);
        Aggregate(Sphere);
        Aggregate(Aggregate Agg1, Aggregate Agg2);
        Aggregate(ListAggregat&, const int label);

        /** Copy constructor */
        Aggregate(const Aggregate&);

        /** Move constructor */
        Aggregate (Aggregate&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~Aggregate(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        Aggregate& operator= (const Aggregate& other);

        /** Move assignment operator */
        Aggregate& operator= (Aggregate&& other) noexcept;
};
#endif // AGGREGAT_H
