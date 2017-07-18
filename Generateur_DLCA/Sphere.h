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
        Sphere(SphereList* Storage,const int i);
        ~Sphere(void);
        double& operator[](const int i);

        void Update(const double newx, const double newy, const double newz, const double newr);
        void Update(const double* newp, const double newr);
        void Update(const Sphere& c);
        void SetLabel(const int value);
        void DecreaseLabel(void);
        void Translate(const double* trans);
        std::string str(const double coef) const;
        void Aff(const double coef) const;
        double Volume(void) const;
        double Surface(void) const;
        double Radius(void) const;
        const double* Position(void) const;

        double Distance(const Sphere& c) const;
        double Distance(const double* point) const;
        double Distance(const double otherx, const double othery, const double otherz) const;

        double Intersection(const Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) const;
        double Collision(const Sphere& c, const double* vd,const  double  distmax,double& distance_contact) const;
    };


class SphereList
{
    friend class Sphere;

    private:
        int N;
        Sphere** spheres;
        PhysicalModel* physicalmodel;

        std::array< std::vector<double>, 7>* Storage;
        SphereList* external_storage;

        void setpointers();

    public:
        ~SphereList(void);
        Sphere& operator[](const int i);

        int size() const;

        void Init(const int N,PhysicalModel& _physicalmodel);
        void extract(const int, int** AggLabels, SphereList& res) const;
        void extractplus(const int, int** AggLabels,const int NAgg, SphereList& res) const;

        void CroissanceSurface(const double dt);

};
#endif // SPHERE
