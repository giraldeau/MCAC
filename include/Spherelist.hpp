#ifndef SPHERELIST_H
#define SPHERELIST_H 1

#include "io/threaded_io.hpp"
#include "Sphere.hpp"
#include "physical_model.hpp"
#include "storagelist.hpp"




/*

Sphere.h and Sphere.cpp defines the data storage.

 * Sphere *
 This is an object representing a sphere (!) with its volume and surface
 Data can be stored in an external Aggregat for vectorization purposes

 Beyond managing its properties, it can compute
 - its distance to a point or an another sphere
 - its intersection with an another sphere (volume and surface) (TODO : NOT CORRECTLY COMPUTED)
 - detect a collision with an another sphere

 * Aggregat *
 This is container for an aggregat which is an enhanced list of spheres
 Data can be shared between multiple Aggregat

*/
namespace MCAC{

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;
class Verlet;

class ListSphere : public storage_list<9,Sphere>
{
    friend class Sphere;
    friend class Aggregate;
    friend class StatisticStorage;

    /* Generic */

    private:
        PhysicalModel* physicalmodel;

        std::vector<double>::iterator ptr_deb;
        std::vector<double>::iterator ptr_fin;

        ThreadedIO* Writer;
        size_t lastSaved;


    public:

        void Init(PhysicalModel& _physicalmodel, size_t _size);
        void DecreaseLabel() noexcept;

        void CroissanceSurface(double dt);

        void print() const;

        void save(){save(false);};
        void save(bool _finish);

        auto get_data() const;

        std::vector<double> Format_Position() const;
        std::vector<double> Format_Radius() const;
        std::vector<long>    Format_AggLabel() const;

    /* Storage specific */
    private:
        void setpointers();

    public:
        /** Default constructor in local storage */
        ListSphere();
        explicit ListSphere(PhysicalModel& _physicalmodel, size_t _size);

        /** Constructor with external storage */
        explicit ListSphere(ListSphere& parent,std::vector<size_t>& _index);
        explicit ListSphere(ListSphere& parent,std::vector<size_t> _index);

        /** Copy constructor */
        ListSphere(const ListSphere& other);
        ListSphere(const ListSphere& other, ListSphere& _Storage);


        /** Move constructor */
        ListSphere (ListSphere&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~ListSphere() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        ListSphere& operator= (const ListSphere& other);

        /** Move assignment operator */
        ListSphere& operator= (ListSphere&& other) noexcept;

        friend bool operator==(const ListSphere&, const ListSphere&);
        friend bool operator!=(const ListSphere&, const ListSphere&);

};

std::string filename(int, int);
}  // namespace MCAC

#endif // SPHERELIST_H
