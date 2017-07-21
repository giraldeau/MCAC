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
class Aggregat;
class ListSphere;
class Sphere;


class Sphere
{
    friend class ListSphere;
    friend class Aggregat;

    /* Generic */

    private:
        double* x;
        double* y;
        double* z;
        double* r;
        double* volume;
        double* surface;
        int AggLabel, SphereLabel;
        PhysicalModel* physicalmodel;

        void UpdateVolAndSurf(void);

    public:
        Sphere(void);
        Sphere(ListSphere& Storage);
        Sphere(PhysicalModel&);
        ~Sphere(void);

        Sphere(PhysicalModel&, const double x, const double y, const double z, const double r);
        Sphere(PhysicalModel&, const double* position, const double newr);
        Sphere(Sphere&);


        void Init(const double x, const double y, const double z, const double r);
        void Init(const double* position, const double r);
        void Init(Sphere&);

        void SetLabel(const int);
        void DecreaseLabel(void);
        void Translate(const double* vector);
        std::string str(const double coef) ;
        void Aff(const double coef) ;
        double Volume(void) ;
        double Surface(void) ;
        double Radius(void) ;
        const double* Position(void) ;

        double Distance(Sphere&) ;
        double Distance(const double* point) ;
        double Distance(const double x, const double y, const double z) ;

        double Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) ;
        double Collision(Sphere& c, const double* vector,const  double  distmax) ;
        void CroissanceSurface(const double dt);

    /* Storage specific */

    private:
        std::array< std::vector<double>, 7>* Storage;
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
        std::vector < Sphere* > spheres;
        PhysicalModel* physicalmodel;

    public:
        void CroissanceSurface(const double dt);

    /* Storage specific */
    private:
        std::array< std::vector<double>, 7>* Storage;
        const ListSphere* external_storage;
        std::vector < int > index;

        void setpointers();

    public:
        ListSphere(void);
        ListSphere(PhysicalModel& _physicalmodel, const int N);
        ListSphere(ListSphere& parent,int* index);
        ListSphere(ListSphere& parent,int** index,const int start,const int end);

        ~ListSphere(void);
        Sphere& operator[](const int);

        int size() const;


};
#endif // SPHERE
