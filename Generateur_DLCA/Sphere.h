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
	
	void Update(double newx, double newy, double newz, double newr);
	void Update(double* newp, double newr);
	void SetLabel(int value);
	void Aff(double coef);
	double Intersection(Sphere c, double* vd, double  distmax,double& distance_contact);
	void Translate(double* trans);
	void Translate(double trans);
	double Distance(Sphere c);
    double VolumeCalotteij(Sphere c);
    double SurfaceCalotteij(Sphere c);
};
