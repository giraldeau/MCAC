#include <iostream>
#include <physical_model.h>
#include <vector>
#include <array>
#ifndef SPHERE
#define SPHERE


class SphereList;
class Sphere;


class Sphere
{
    friend class SphereList;

    private:
        double* x;
        double* y;
        double* z;
        double* r;
        double* volume;
        double* surface;

        int AggLabel, SphereLabel;

        std::array< std::vector<double>, 7>* Storage;
        SphereList* external_storage;

        void UpdateVolAndSurf(void);
        void setpointers(void);
        void add(void);

    public:
        Sphere(void);
        Sphere(SphereList* Storage);
        ~Sphere(void);
        double& operator[](const int i);

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

        double Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) ;
        double Collision(Sphere& c, const double* vd,const  double  distmax,double& distance_contact) ;
    };


class SphereList
{
    friend class Sphere;

    private:
        int N;
        Sphere** spheres;
        PhysicalModel* physicalmodel;

        std::array< std::vector<double>, 7>* Storage;
        const SphereList* external_storage;
        int* index;

        void setpointers();

    public:
        SphereList(void);
        SphereList(const SphereList* parent,int* AggLabels);

        ~SphereList(void);
        Sphere& operator[](const int i);

        int size() const;

        void Init(const int N,PhysicalModel& _physicalmodel);

        SphereList extract(const int, int** AggLabels) const;
        SphereList extractplus(const int, int** AggLabels,const int NAgg) const;

        void CroissanceSurface(const double dt);

};
#endif // SPHERE
