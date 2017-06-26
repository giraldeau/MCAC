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
    void SetLabel(const int value);
    void Translate(const double* trans);
    void Translate(const double trans);
    void Aff(const double coef) const;
    double Distance(const Sphere& c) const;
    double Intersection(const Sphere& c, const double* vd,const  double  distmax,double& distance_contact) const;
    double VolumeCalotteij(const Sphere& c) const;
    double SurfaceCalotteij(const Sphere& c) const;
};
