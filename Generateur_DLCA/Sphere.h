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

        array< vector<double>, 7>* Storage;
        ListSphere* external_storage;

        void UpdateVolAndSurf(void);

    public:

        void Init(const double x, const double y, const double z, const double r);
        void Init(const double position[], const double r);
        void Init(const array<double, 4> position,const double r);

        void Init(Sphere&);
        void Init(void);

        void SetPosition(const double x, const double y, const double z);
        void SetPosition(const double position[]);
        void SetPosition(const array<double, 4> position);

        void SetLabel(const int);
        void DecreaseLabel(void);
        void Translate(const double x,const double y,const double z);
        void Translate(const double vector[]);
        void Translate(const array<double, 4>  vector);
        string str(const double coef) ;
        void Aff(const double coef) ;
        double Volume(void) ;
        double Surface(void) ;
        double Radius(void) ;
        const array<double, 4> Position(void) ;

        double Distance(Sphere&) ;
        double Distance(const double point[]) ;
        double Distance(const array<double, 4> point) ;
        double Distance(const double x, const double y, const double z) ;

        double Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) ;
        double Collision(Sphere& c, const array<double,4> vector,const  double  distmax) ;
        void CroissanceSurface(const double dt);

    /* Storage specific */

    public:
        Sphere(void);
        Sphere(ListSphere& Storage, const int id);
        Sphere(PhysicalModel&);
        ~Sphere(void);

        Sphere(PhysicalModel&, const double x, const double y, const double z, const double r);
        Sphere(PhysicalModel&, const array<double, 4> position, const double r);
        Sphere(PhysicalModel&, const double position[], const double r);
        Sphere(const Sphere&);

    private:

        void setpointers(void);
        void add(void);
        double& operator[](const int);
        Sphere& operator=(const Sphere& other);
        Sphere& operator=(const Sphere&& other) noexcept;

    };


class ListSphere
{
    friend class Sphere;

    /* Generic */

    private:
        PhysicalModel* physicalmodel;
        vector < Sphere* > spheres;
        vector < int > index;

        array< vector<double>, 7>* Storage;
        const ListSphere* external_storage;

        int N;

    public:
        void CroissanceSurface(const double dt);
        void DecreaseLabel(void);



    /* Storage specific */
    private:
        void setpointers();

    public:
        ListSphere(void);
        ListSphere(PhysicalModel& _physicalmodel, const int N);
        ListSphere(ListSphere& parent,int index[]);
        ListSphere(ListSphere& parent,int* index[],const int start,const int end);

        ~ListSphere(void);
        Sphere& operator[](const int);

        void Init(PhysicalModel& _physicalmodel, const int N);
        void Destroy();

        int size() const;
        ListSphere(const ListSphere& other);

        ListSphere& operator=(ListSphere other);
        friend void swap(ListSphere& first, ListSphere& second);

};

double periodicPosition(const double x, const double dim);
double periodicDistance(const double x, const double dim);

#endif // SPHERE
