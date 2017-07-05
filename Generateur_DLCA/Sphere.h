#include <iostream>

#pragma once

class Sphere
{
private:


public:
    double* pos;
    double r;
    int Label;
    Sphere(void);
    ~Sphere(void);
	
    void Update(const double newx, const double newy, const double newz, const double newr);
    void Update(const double* newp, const double newr);
    void Update(const Sphere& c);
    void SetLabel(const int value);
    void DecreaseLabel();
    void Translate(const double* trans);
    void Translate(const double trans);
    std::string str(const double coef) const;
    void Aff(const double coef) const;
    double Distance(const Sphere& c) const;
    double Distance(const double* point) const;
    double Volume() const;
    double Surface() const;
    double Radius() const;
    double* Position() const;
    double Intersection(const Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) const;
    double Collision(const Sphere& c, const double* vd,const  double  distmax,double& distance_contact) const;
};
