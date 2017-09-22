#ifndef AGGREGAT_H
#define AGGREGAT_H

#include "Sphere.hpp"
#include "Spherelist.hpp"
#include "aggregatList.hpp"
#include "storage.hpp"
#include "storagelist.hpp"

namespace DLCA{

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;
class Verlet;
class AnalizedAggregate;

class Aggregate : public storage_elem<15,ListAggregat>
{

    friend class ListAggregat;
    friend class AnalizedAggregate;
    friend class StatisticStorage;
    friend class Sphere;

    private:
        PhysicalModel* physicalmodel;

        ListSphere myspheres;

        Verlet* verlet;
        std::array<size_t, 3> IndexVerlet;

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

        // for statistics filtering
        double *Dp,*DgOverDp;

        size_t Np;     //Number of spheres

        bool InVerlet;

        /*
        bool padding2,padding3,padding4;
        int padding1;
        */

    public:

        void Init();
        void Init(PhysicalModel& ,Verlet&, std::array<double, 3> position ,size_t _label, ListSphere&, double D);

        double GetLpm() const noexcept;
        double GetVolAgregat() const noexcept;
        size_t GetLabel() const noexcept;


        const std::array<double, 3> GetPosition() const noexcept;
        std::array<size_t, 3> GetVerletIndex() noexcept;


        void SetPosition(std::array<double, 3> position) noexcept;
        void SetPosition(double newx, double newy, double newz) noexcept;
        void Translate(std::array<double, 3> vector) noexcept;


        void Update();
        void Volume();
        void MassCenter();
        void CalcRadius();
        void RayonGiration();
        bool Contact(Aggregate&) const noexcept;
        double Distance(Aggregate&, std::array<double,3> Vectdir) const;
        void Merge(Aggregate&);
        void DecreaseLabel() noexcept;

        void UpdateDistances() noexcept;

        void check();

        // for statistics filtering
        void Statistics();
        bool operator <(const Aggregate&) const;

        /* Storage specific */
    private:
        void setpointers();

    public:
        /** Default constructor in local storage */
        Aggregate();
        explicit Aggregate(PhysicalModel&);

        /** Constructor in local storage with initialization */
        explicit Aggregate(PhysicalModel&, double x, double y, double z, double r);
        explicit Aggregate(PhysicalModel&, std::array<double, 3> position, double r);
        explicit Aggregate(Sphere&);

        /** Constructor with external storage */
        explicit Aggregate(ListAggregat&, size_t label);

        /** Copy constructor */
        Aggregate(const Aggregate&);

        /** Move constructor */
        Aggregate (Aggregate&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~Aggregate() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        Aggregate& operator= (const Aggregate& other);

        /** Move assignment operator */
        Aggregate& operator= (Aggregate&& other) noexcept;
};
}  // namespace DLCA


#endif // AGGREGAT_H
