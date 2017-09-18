#ifndef SPHERELIST_H
#define SPHERELIST_H

#include "Sphere.h"
#include "Spherelist.h"
#include "storagelist.h"
#include "physical_model.h"
#include "IO.h"


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

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;
class Verlet;

class ListSphere : public storage_list<9,Sphere>
{
    friend class Sphere;

    /* Generic */

    private:
        PhysicalModel* physicalmodel;

        std::vector<double>::iterator ptr_deb;
        std::vector<double>::iterator ptr_fin;

        ThreadedIO* Writer;

    public:

        void Init(PhysicalModel& _physicalmodel, const int _size);
        void DecreaseLabel(void) noexcept;

        void CroissanceSurface(const double dt);

        void save(void) const;
        void save(const bool finish) const;
        std::vector<double> FormatPositionData() const;
        std::vector<double> FormatRadiusData() const;
        std::vector<int>    FormatLabelData() const;

    /* Storage specific */
    private:
        void setpointers(void);

    public:
        /** Default constructor in local storage */
        ListSphere(void);
        ListSphere(PhysicalModel& _physicalmodel, const int _size);

        /** Constructor with external storage */
        ListSphere(ListSphere& parent,std::vector<int> ListIndex);

        /** Copy constructor */
        ListSphere(const ListSphere& other);

        /** Move constructor */
        ListSphere (ListSphere&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~ListSphere(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        ListSphere& operator= (const ListSphere& other);

        /** Move assignment operator */
        ListSphere& operator= (ListSphere&& other) noexcept;

        friend bool operator==(const ListSphere&, const ListSphere&);
        friend bool operator!=(const ListSphere&, const ListSphere&);

};

std::string filename(const int, const int);

#endif // SPHERELIST_H
