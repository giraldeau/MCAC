#ifndef AGGREGAT_H
#define AGGREGAT_H 1

#include "spheres/sphere.hpp"
#include "spheres/sphere_list.hpp"
#include "aggregatList.hpp"
#include "statistics.hpp"
#include "storage.hpp"
#include "storagelist.hpp"

#include <list>

namespace MCAC{

class Sphere;
class SphereList;
class Aggregate;
class ListAggregat;
class Verlet;
class AnalizedAggregate;
class StatisicsData;

class Aggregate :
        public storage_elem<15,ListAggregat>,
        public StatisicsData
{

    friend class ListAggregat;
    friend class AnalizedAggregate;
    friend class StatisticStorage;
    friend class Sphere;
    friend class storage_list<15,Aggregate>;


    using storage_elem<15,ListAggregat>::external_storage;
    using storage_elem<15,ListAggregat>::Storage;
    using storage_elem<15,ListAggregat>::indexInStorage;


    private:
        PhysicalModel* physicalmodel;

        SphereList myspheres;

        Verlet* verlet;
        std::array<size_t, 3> IndexVerlet;

        std::vector<std::list < std::pair<size_t, double > > > _distances;
        std::vector<double> distances_center;
        std::vector<double> volumes;
        std::vector<double> surfaces;

        double *rg;                     //Gyration Radius
        double *f_agg;                  //Friction coeff
        double *lpm;                    //Apparent Mean Free Path
        double *time_step;              //Time to move along lpm
        double *rmax;                   //Radius of the sphere containing Agg
        double *volAgregat;             //Etimation of the Aggregate's volume
        double *surfAgregat;            //Estimation of the sufrace of the aggregate

        double *x,*y,*z;                // position of the gravity center
        double *rx,*ry,*rz;             // position of the gravity center
        double *time;                   // Proper time of the aggregate

        size_t Np;                      //Number of spheres
        size_t Label;


        bool InVerlet;

        /*
        bool padding2,padding3,padding4;
        int padding1;
        */

    public:

        void Init();
        void Init(PhysicalModel&, Verlet&, std::array<double, 3> position , size_t _label, SphereList&, double D);

        double SphereDistance(size_t i, size_t j) const;
        double SphereDistance(size_t i) const;
        double GetLpm() const noexcept;
        double GetTimeStep() const noexcept;
        double GetVolAgregat() const noexcept;
        size_t GetLabel() const noexcept;


        std::array<double, 3> GetPosition() const noexcept;
        std::array<size_t, 3> GetVerletIndex() noexcept;


        void SetPosition(std::array<double, 3> position) noexcept;
        void SetPosition(double newx, double newy, double newz) noexcept;
        void Translate(std::array<double, 3> vector) noexcept;
        void TimeForward(double deltatemps) noexcept;


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

        void print() const ;

        void partialStatistics();
        void fullStatistics();

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
        Aggregate(const Aggregate&, ListAggregat&);

        /** Move constructor */
        Aggregate (Aggregate&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~Aggregate() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        Aggregate& operator= (const Aggregate& other);

        /** Move assignment operator */
        Aggregate& operator= (Aggregate&& other) noexcept;
};
}  // namespace MCAC


#endif // AGGREGAT_H