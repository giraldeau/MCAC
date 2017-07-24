#ifndef SPHERE
#define SPHERE

#include <iostream>
#include <physical_model.h>
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

class Sphere
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

        int AggLabel, SphereLabel;

        void UpdateVolAndSurf(void);

    public:

        void Init(const double x, const double y, const double z, const double r);
        void Init(const double* position, const double r);
        void Init(const array<double, 4> position,const double r);

        void Init(Sphere&);
        void Init(void);


        void SetLabel(const int);
        void DecreaseLabel(void);
        void Translate(const double* vector);
        void Translate(const array<double, 4>  vector);
        string str(const double coef) ;
        void Aff(const double coef) ;
        double Volume(void) ;
        double Surface(void) ;
        double Radius(void) ;
        const array<double, 4> Position(void) ;

        double Distance(Sphere&) ;
        double Distance(const double* point) ;
        double Distance(const array<double, 4> point) ;
        double Distance(const double x, const double y, const double z) ;

        double Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) ;
        double Collision(Sphere& c, const double* vector,const  double  distmax) ;
        void CroissanceSurface(const double dt);

    /* Storage specific */

    public:
        Sphere(void);
        Sphere(ListSphere& Storage, const int id);
        Sphere(PhysicalModel&);
        ~Sphere(void);

        Sphere(PhysicalModel&, const double x, const double y, const double z, const double r);
        Sphere(PhysicalModel&, const array<double, 4> position, const double r);
        Sphere(PhysicalModel&, const double* position, const double r);
        Sphere(Sphere&);

    private:
        array< vector<double>, 7>* Storage;
        ListSphere* external_storage;

        void setpointers(void);
        void add(void);
        double& operator[](const int);

    };


class ListSphere
{
    friend class Sphere;

    /* Generic */

    private:
        int N;
        vector < Sphere* > spheres;
        PhysicalModel* physicalmodel;

    public:
        void CroissanceSurface(const double dt);

    /* Storage specific */
    private:
        array< vector<double>, 7>* Storage;
        const ListSphere* external_storage;
        vector < int > index;

        void setpointers();

    public:
        ListSphere(void);
        ListSphere(PhysicalModel& _physicalmodel, const int N);
        ListSphere(ListSphere& parent,int* index);
        ListSphere(ListSphere& parent,int** index,const int start,const int end);

        ~ListSphere(void);
        Sphere& operator[](const int);

        void Init(PhysicalModel& _physicalmodel, const int N);
        void Destroy();

        int size() const;


};
#endif // SPHERE
