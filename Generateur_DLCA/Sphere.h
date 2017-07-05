#include <iostream>
#include <physical_model.h>
#ifndef SPHERE
#define SPHERE

class Sphere
{
    friend class SphereList;

    private:
        double* pos;
        double r;
        int Label;

        void Update(void);

        double volume, surface;

    public:
        Sphere(void);
        ~Sphere(void);

        void Update(const double newx, const double newy, const double newz, const double newr);
        void Update(const double* newp, const double newr);
        void Update(const Sphere& c);
        void SetLabel(const int value);
        void DecreaseLabel(void);
        void Translate(const double* trans);
        std::string str(const double coef) const;
        void Aff(const double coef) const;
        double Distance(const Sphere& c) const;
        double Distance(const double* point) const;
        double Volume(void) const;
        double Surface(void) const;
        double Radius(void) const;
        const double* Position(void) const;
        double Intersection(const Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) const;
        double Collision(const Sphere& c, const double* vd,const  double  distmax,double& distance_contact) const;
    };


class SphereList
{
    private:
        Sphere* spheres;
        int N;
        const PhysicalModel* physicalmodel;

    public:

        void Init(const int N,const PhysicalModel _physicalmodel);
        ~SphereList(void);

        Sphere& operator[](const int i);
        void CroissanceSurface(double dt);
};
#endif // SPHERE
