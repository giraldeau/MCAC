#include "Sphere.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double PI = atan(1.0)*4;

Sphere::Sphere(void)
{
    pos = new double[4];
    Label = 1;
    Update(0, 0, 0, 1);
}

Sphere::~Sphere(void)
{
}

void Sphere::Update(double newx, double newy, double newz, double newr)
{
    pos[1] = newx;
    pos[2] = newy;
    pos[3] = newz;
    r = newr;
}

void Sphere::SetLabel(int value)
{
    Label = value;
}

void Sphere::Aff(double coef)
{
    printf("%8.3f %8.3f %8.3f %8.3f\n", pos[1]*coef, pos[2]*coef, pos[3]*coef, r*coef);
}

void Sphere::Update(double* newp, double newr)
{	
	// newp: position double[4] pour compatibilité....
    pos[1] = newp[1];
    pos[2] = newp[2];
    pos[3] = newp[3];
    r = newr;
}

double Sphere::Distance(Sphere c)
{
	return sqrt(pow(pos[1]-c.pos[1],2)+pow(pos[2]-c.pos[2],2)+pow(pos[3]-c.pos[3],2));
}

double Sphere::Intersection(Sphere c, double* vd,double distmax,double& distance_contact)
{
/*
     (vd): vecteur directeur double[4] : vd[1],vd[2],vd[3], vd[0] inutilisé
     résultat:
        dist=-1: les sphères se touchent
        dist=1 : la sphère courante ne peut pas toucher la sphère (c) en
        se déplaçant suivant le vecteur directeur (vd) et sur une distance maxi
        égale à (distmax).
        sinon  : il peut y avoir contact. La sphère courante doit être déplacée
        de (dist) en suivant (vd) pour toucher la sphère (c)
*/
  /*! 
   *  \image html cCSphereFIntersectionCSphereddd.png
   */


    double A, B, C;
    double dx, dy, dz;
    double DELTA;
    double K, K1, K2;
    double dist;

    dx = c.pos[1] - pos[1];
    dy = c.pos[2] - pos[2];
    dz = c.pos[3] - pos[3];

    //$ Compute signed distance for contact between two spheres
    distance_contact = sqrt(dx*dx+dy*dy+dz*dz)-r-c.r;

    if (distance_contact <= 0)
    {
       //$ There is contact between these two spheres
       dist = -1;
    }
    else
    {
       //$ Computing distance before contact

        dist = 1;
        K = -1;
        A = vd[1]*vd[1]+vd[2]*vd[2]+vd[3]*vd[3];
        B = -2*(dx*vd[1]+dy*vd[2]+dz*vd[3]);
        C = -pow(r+c.r,2)+(dx*dx+dy*dy+dz*dz);
        DELTA = B*B-4*A*C;
        if (DELTA >= 0)
        {
            DELTA = sqrt(DELTA);
            K1=(-B-DELTA)/2/A;
            K2=(-B+DELTA)/2/A;
            if (K1*K2 < 0) K = fmax(K1,K2);
            else if (K1 > 0) K = fmin(K1,K2);
            if (K > 0) dist = K*sqrt(A);
            if (dist > distmax) dist = 1;
        }
    }
    return dist;
}

void Sphere::Translate(double* trans)
{
    for (int i = 1; i <= 3; i++) pos[i] = pos[i] + trans[i];
}

void Sphere::Translate(double trans)
{
    for (int i = 1; i <= 3; i++) pos[i] = pos[i] + trans;
}

//Calcul du volume de la calotte sphérique de la sphère courante de rayon Ri due à la surestimation de la sphère c de rayon Rj
double Sphere::VolumeCalotteij(Sphere c)
{
    double Volcal, d, h;
    double Ri, Rj;

    Ri = r;
    Rj = c.r;

    d = sqrt(pow(pos[1]-c.pos[1],2.0) + pow(pos[2]-c.pos[2],2.0) + pow(pos[3]-c.pos[3],2.0));

    if (d >= Ri + Rj)
        Volcal = 0.0;
    else if ((0 <= d) && (d < Ri - Rj))
        Volcal = 4.0*PI*pow(Rj,3.0)/3.0;
    else if ((0 <= d) && (d < Rj - Ri))
        Volcal = 4.0*PI*pow(Ri,3.0)/3.0;
    else
    {
        h = (pow(Rj,2.0)-pow((Ri-d),2.0))/(2.0*d);
        Volcal = PI*pow(h,2.0)*(3*Ri-h)/3.0;
    }

    return Volcal;
}

//Calcul de la surface de la calotte sphérique de la sphère courante de rayon Ri due à la surestimation de la sphère c de rayon Rj
double Sphere::SurfaceCalotteij(Sphere c)
{
    double Surfcal, d, h;
    double Ri, Rj;

    Ri = r;
    Rj = c.r;

    d = sqrt(pow(pos[1]-c.pos[1], 2.0)+pow(pos[2]-c.pos[2], 2.0)+pow(pos[3]-c.pos[3], 2.0));

    if (d >= Ri + Rj)
        Surfcal = 0.0;
    else if ((0 <= d) && (d < Ri - Rj))
        Surfcal = 4.0*PI*pow(Rj, 2.0);
    else if ((0 <= d) && (d < Rj - Ri))
        Surfcal = 4.0*PI*pow(Ri, 2.0);
    else
    {
        h = (pow(Rj,2.0)-pow((Ri-d),2.0))/(2.0*d);
        Surfcal = 2*PI*Ri*h;
    }

    return Surfcal;
}
