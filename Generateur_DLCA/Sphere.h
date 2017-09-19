#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include "Sphere.h"
#include "Spherelist.h"
#include "storage.h"
#include "physical_model.h"
#include <vector>
#include <array>

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


namespace DLCA{

class Sphere;
class ListSphere;
class Aggregate;
class ListAggregat;
class Verlet;

class Sphere : public storage_elem<9,ListSphere>
{
    friend class ListSphere;
    friend class Aggregate;

    /* Generic */

    private:
        double* x;
        double* y;
        double* z;
        double* r;
        double* rx;
        double* ry;
        double* rz;
        double* volume;
        double* surface;
        PhysicalModel* physicalmodel;

        int AggLabel;

        void UpdateVolAndSurf(void) noexcept;

    public:
        void InitVal(void);
        void InitVal(const double x, const double y, const double z, const double r);
        void InitVal(const std::array<double, 3> position,const double r);


        void SetPosition(const double x, const double y, const double z) noexcept;
        void SetPosition(const std::array<double, 3> position) noexcept;

        void Translate(const double x,const double y,const double z) noexcept;
        void Translate(const std::array<double, 3>  vector) noexcept;
        void RelativeTranslate(const double x,const double y,const double z) noexcept;

        void SetLabel(const int) noexcept;
        void DecreaseLabel(void) noexcept;

        void CroissanceSurface(const double dt);

        double Volume(void) const noexcept;
        double Surface(void) const noexcept;
        double Radius(void) const noexcept;
        const std::array<double, 3> Position(void) const noexcept;

        std::string str(const double coef) const ;
        void Aff(const double coef) const ;

        double Distance(const Sphere&) const noexcept;
        double Distance(const std::array<double, 3> point) const noexcept;
        double Distance(const double x, const double y, const double z) const noexcept;
        double RelativeDistance(const double x, const double y, const double z) const noexcept;
        double RelativeDistance(const Sphere&) const noexcept;


        double Distance2(const Sphere&) const noexcept;
        double Distance2(const std::array<double, 3> point) const noexcept;
        double Distance2(const double x, const double y, const double z) const noexcept;
        double RelativeDistance2(const Sphere&) const noexcept;
        double RelativeDistance2(const double x, const double y, const double z) const noexcept;



        bool Contact(const Sphere&) const noexcept;

        double Collision(const Sphere&, const std::array<double,3> vector) const;
        double Intersection(const Sphere& c, const double dist,double& vol1, double& vol2, double& surf1, double& surf2 ) const;

    /* Storage specific */
    private:
        void setpointers(void);

    public:
        /** Default constructor in local storage */
        Sphere(void);
        Sphere(PhysicalModel&);
        Sphere(const Aggregate&);

        /** Constructor in local storage with initialization */
        Sphere(PhysicalModel&, const double x, const double y, const double z, const double r);
        Sphere(PhysicalModel&, const std::array<double, 3> position, const double r);

        /** Constructor with external storage */
        Sphere(ListSphere& Storage, const int id);

        /** Copy constructor */
        Sphere(const Sphere&);

        /** Move constructor */
        Sphere (Sphere&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~Sphere(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        Sphere& operator= (const Sphere& other);

        /** Move assignment operator */
        Sphere& operator= (Sphere&& other) noexcept;

};

}
#endif // SPHERE_H
