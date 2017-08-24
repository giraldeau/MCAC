#ifndef SPHERE
#define SPHERE

#include <iostream>
#include <physical_model.h>
#include <storage.h>
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

using namespace std;

class Aggregate;
class ListSphere;
class Sphere;

class Sphere : public storage_elem<7,ListSphere>
{
    friend class ListSphere;
    friend class Aggregate;

    /* Generic */

    private:
        double* x;
        double* y;
        double* z;
        double* r;
        double* volume;
        double* surface;
        PhysicalModel* physicalmodel;

        int AggLabel;

        void UpdateVolAndSurf(void);

    public:

        void Set(const double x, const double y, const double z, const double r);
        void Set(const double position[], const double r);
        void Set(const array<double, 4> position,const double r);

        void Copy(Sphere&);
        void Init(void);

        void SetPosition(const double x, const double y, const double z);
        void SetPosition(const double position[]);
        void SetPosition(const array<double, 4> position);

        void Translate(const double x,const double y,const double z);
        void Translate(const double vector[]);
        void Translate(const array<double, 4>  vector);

        void SetLabel(const int);
        void DecreaseLabel(void);

        double Volume(void) ;
        double Surface(void) ;
        double Radius(void) ;
        const array<double, 4> Position(void) ;

        string str(const double coef) ;
        void Aff(const double coef) ;

        double Distance(Sphere&) ;
        double Distance(const double point[]) ;
        double Distance(const array<double, 4> point) ;
        double Distance(const double x, const double y, const double z) ;

        double Distance2(Sphere&) ;
        double Distance2(const double point[]) ;
        double Distance2(const array<double, 4> point) ;
        double Distance2(const double x, const double y, const double z) ;

        double Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) ;
        double Collision(Sphere& c, const array<double,4> vector,const  double  distmax) ;
        void CroissanceSurface(const double dt);

    /* Storage specific */
    private:
        void setpointers(void);

    public:
        /** Default constructor in local storage */
        Sphere(void);
        Sphere(PhysicalModel&);

        /** Constructor in local storage with initialization */
        Sphere(PhysicalModel&, const double x, const double y, const double z, const double r);
        Sphere(PhysicalModel&, const array<double, 4> position, const double r);
        Sphere(PhysicalModel&, const double position[], const double r);

        /** Constructor with external storage */
        Sphere(ListSphere& Storage, const int id);

//        /** Copy constructor */
//        Sphere(const Sphere&);

//        /** Move constructor */
//        Sphere (Sphere&&) noexcept; /* noexcept needed to enable optimizations in containers */

//        /** Destructor */
//        ~Sphere(void) noexcept; /* explicitly specified destructors should be annotated noexcept as best-practice */

//        /** Copy assignment operator */
//        Sphere& operator= (const Sphere& other);

//        /** Move assignment operator */
//        Sphere& operator= (Sphere&& other) noexcept;

    };


class ListSphere : public storage_list<7,Sphere>
{
    friend class Sphere;

    /* Generic */

    private:
        PhysicalModel* physicalmodel;

    public:

        void Init(PhysicalModel& _physicalmodel, const int size);
        void DecreaseLabel(void);

        void CroissanceSurface(const double dt);


    /* Storage specific */
    private:
        void setpointers(void);

    public:
        /** Default constructor in local storage */
        ListSphere(void);
//        ListSphere(PhysicalModel& _physicalmodel, const int size);

        /** Constructor with external storage */
        ListSphere(ListSphere& parent,int index[]);
        ListSphere(ListSphere& parent,int* index[],const int start,const int end);

//        /** Copy constructor */
//        ListSphere(const ListSphere& other);

//        /** Move constructor */
//        ListSphere (ListSphere&&) noexcept; /* noexcept needed to enable optimizations in containers */

//        /** Copy assignment operator */
//        ListSphere& operator= (const ListSphere& other);

//        /** Move assignment operator */
//        ListSphere& operator= (ListSphere&& other) noexcept;

//        friend void swap(ListSphere& first, ListSphere& second);

};

double periodicPosition(const double x, const double dim);
double periodicDistance(const double x, const double dim);

#endif // SPHERE
