#ifndef SPHERE_H
#define SPHERE_H

#include "physical_model.hpp"
#include "Spherelist.hpp"
#include "storage.hpp"

#include <array>
#include <iostream>
#include <vector>
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

class Sphere : public storage_elem<9,ListSphere>
{
    friend class ListSphere;
    friend class Aggregate;
    friend class ListAggregat;

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

        long AggLabel;

        void UpdateVolAndSurf() noexcept;

    public:
        void InitVal();
        void InitVal(double newx, double newy, double newz, double newr);
        void InitVal(std::array<double, 3> newposition, double newr);


        void SetPosition(double newx, double newy, double newz) noexcept;
        void SetPosition(std::array<double, 3> newposition) noexcept;

        void Translate(double transx,double transy, double transz) noexcept;
        void Translate(std::array<double, 3>  trans) noexcept;
        void RelativeTranslate(double transx, double transy, double transz) noexcept;

        void SetLabel( long) noexcept;
        void DecreaseLabel() noexcept;

        void CroissanceSurface( double dt);

        double Volume() const noexcept;
        double Surface() const noexcept;
        double Radius() const noexcept;
        const std::array<double, 3> Position() const noexcept;

        std::string str(double coef) const ;
        void Aff(double coef) const ;
        void print() const;

        double Distance(const Sphere&) const noexcept;
        double Distance(std::array<double, 3> point) const noexcept;
        double Distance(double otherx, double othery, double otherz) const noexcept;
        double RelativeDistance(double otherx, double othery, double otherz) const noexcept;
        double RelativeDistance(const Sphere&) const noexcept;


        double Distance2(const Sphere&) const noexcept;
        double Distance2(std::array<double, 3> point) const noexcept;
        double Distance2(double otherx, double othery, double otherz) const noexcept;
        double RelativeDistance2(const Sphere&) const noexcept;
        double RelativeDistance2(double otherx, double othery, double otherz) const noexcept;



        bool Contact(const Sphere&) const noexcept;

        std::pair<bool,double> CollisionR(const Sphere& c, std::array < std::array < double, 3>, 3> RotMat) const;
        std::pair<bool,double> Collision(const Sphere&, std::array<double,3> vectordir) const;
        std::vector<double> Collisions(const ListSphere& list, std::array<double,3> vectordir) const;

        double Intersection(const Sphere& c, double dist,double& vol1, double& vol2, double& surf1, double& surf2 ) const;

    /* Storage specific */
    private:
        void setpointers();

    public:
        /** Default constructor in local storage */
        Sphere();
        explicit Sphere(PhysicalModel&);
        explicit Sphere(const Aggregate&);

        /** Constructor in local storage with initialization */
        explicit Sphere(PhysicalModel&, double newx, double newy, double newz, double newr);
        explicit Sphere(PhysicalModel&, std::array<double, 3> newposition, double newr);

        /** Constructor with external storage */
        explicit Sphere(ListSphere& aggregat, size_t id);

        /** Copy constructor */
        Sphere(const Sphere&);
        Sphere(const Sphere&, ListSphere& aggregat, size_t id);


        /** Move constructor */
        Sphere (Sphere&&) noexcept; /* noexcept needed to enable optimizations in containers */

        /** Destructor */
        ~Sphere() noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

        /** Copy assignment operator */
        Sphere& operator= (const Sphere& other);

        /** Move assignment operator */
        Sphere& operator= (Sphere&& other) noexcept;

};

std::array < std::array < double, 3>, 3> GetRotMat(std::array<double,3> vectordir);


}  // namespace MCAC

#endif // SPHERE_H
