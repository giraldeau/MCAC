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
class Sphere;


class Sphere
{
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
        double Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) ;
        double Collision(Sphere& c, const double* vd,const  double  distmax,double& distance_contact) ;
        void CroissanceSurface(const double dt);

    /* Storage specific */

    private:
        std::array< std::vector<double>, 7>* Storage;
        Aggregat* external_storage;

        void setpointers(void);
        void add(void);
        double& operator[](const int i);

    public:
        Sphere(void);
        Sphere(Aggregat* Storage);
        Sphere(PhysicalModel& _physicalmodel);
        ~Sphere(void);

        void Update(const double newx, const double newy, const double newz, const double newr);
        void Update(const double* newp, const double newr);
        void Update(Sphere& c);
        void SetLabel(const int value);
        void DecreaseLabel(void);
        void Translate(const double* trans);
        std::string str(const double coef) ;
        void Aff(const double coef) ;
        double Volume(void) ;
        double Surface(void) ;
        double Radius(void) ;
        const double* Position(void) ;

        double Distance(Sphere& c) ;
        double Distance(const double* point) ;
        double Distance(const double otherx, const double othery, const double otherz) ;

    };


class Aggregat
{
    friend class Sphere;

    /* Generic */

    private:
        int N;
        Sphere** spheres;
        PhysicalModel* physicalmodel;

    public:
        void CroissanceSurface(const double dt);

    /* Storage specific */
    private:
        std::array< std::vector<double>, 7>* Storage;
        const Aggregat* external_storage;
        int* index;

        void setpointers();

    public:
        Aggregat(void);
        Aggregat(const Aggregat* parent,int* AggLabels);

        ~Aggregat(void);
        Sphere& operator[](const int i);

        int size() const;

        void Init(const int N,PhysicalModel& _physicalmodel);

        Aggregat extract(const int, int** AggLabels) const;
        Aggregat extractplus(const int, int** AggLabels,const int NAgg) const;

};
#endif // SPHERE
