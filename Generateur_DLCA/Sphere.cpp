#include "Sphere.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

const double PI = atan(1.0)*4;
const double facvol = 4*PI/3;
const double facsurf = 4*PI;

Sphere::Sphere(void)
{
    pos = new double[4];
    Label = 0;
    Update(0, 0, 0, 0);
}

Sphere::~Sphere(void)
{
    if (pos!=NULL) delete[] pos;
}

void Sphere::Update(const double newx,const double newy,const double newz,const double newr)
{
    pos[1] = newx;
    pos[2] = newy;
    pos[3] = newz;
    r = newr;
    Update();
}

void Sphere::Update(const double* newp,const double newr)
{
    Update(newp[1],newp[2],newp[3],newr);
}

void Sphere::Update(const Sphere& c)
{
    Update(c.pos,c.r);
}


void Sphere::SetLabel(const int value)
{
    Label = value;
}

void Sphere::DecreaseLabel()
{
    Label--;
}

void Sphere::Translate(const double* trans)
{
    for (int i = 1; i <= 3; i++) pos[i] = pos[i] + trans[i];
}

double Sphere::Volume() const
{
    return volume;
}

double Sphere::Surface() const
{
    return surface;
}

double Sphere::Radius() const
{
    return r;
}

const double* Sphere::Position() const
{
    return pos;
}

string Sphere::str(const double coef) const
{
    stringstream res;
    res
        << setw(5) << Label << "\t"
        << setprecision(5) << setw(11) << fixed << r*coef << "\t"
        << setprecision(5) << setw(11) << fixed << pos[1]*coef << "\t"
        << setprecision(5) << setw(11) << fixed << pos[2]*coef << "\t"
        << setprecision(5) << setw(11) << fixed << pos[3]*coef;
    return res.str();
}

void Sphere::Aff(const double coef) const
{  
    cout << str(coef) << endl;
}

void Sphere::Update()
{
    if (Label > 0)
    {
        volume = facvol*pow(r, 3);
        surface =  facsurf*pow(r, 2);
    }
}

double Sphere::Distance(const Sphere& c) const
{
    return Distance(c.pos);
}

double Sphere::Distance(const double* point) const
{
    return sqrt(pow(pos[1]-point[1],2)+pow(pos[2]-point[2],2)+pow(pos[3]-point[3],2));
}

double Sphere::Collision(const Sphere& c,const  double* vd,const double distmax,double& distance) const
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
    double B, C;
    double dx, dy, dz;
    double DELTA;
    double K, K1, K2;
    double dist,dist_contact;

    dx = c.pos[1] - pos[1];
    dy = c.pos[2] - pos[2];
    dz = c.pos[3] - pos[3];

    dist_contact = pow(r + c.r,2);

    //$ Compute signed distance for contact between two spheres
    distance = dx*dx+dy*dy+dz*dz;

    if (distance <= dist_contact)
    {
       //$ There is contact between these two spheres
       dist = -1;
    }
    else
    {
       //$ Computing distance before contact

        dist = 1;
        K = -1;
        B = -2*(dx*vd[1]+dy*vd[2]+dz*vd[3]);
        C = distance - dist_contact;
        DELTA = B*B-4*C;
        if (DELTA >= 0)
        {
            DELTA = sqrt(DELTA);
            K1=0.5*(-B-DELTA);
            K2=0.5*(-B+DELTA);
            K = MIN(K1,K2);
            if (K < 0) K = MAX(K1,K2);
            if (K > 0) dist = K;
            if (dist > distmax) dist = 1;
        }
    }
    return dist;
}

//Calcul du volume de la calotte sphérique de la sphère courante de rayon Ri due à la surestimation de la sphère c de rayon Rj
double Sphere::Intersection(const Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) const
{
    double d, h;
    double Ri, Rj;

    vol1 = vol2 = 0.;
    surf1 = surf2 = 0.;

    Ri = r;
    Rj = c.r;

    //$ Determination of the distance between the center of the 2 aggregates
    d = sqrt(pow(pos[1]-c.pos[1],2) + pow(pos[2]-c.pos[2],2) + pow(pos[3]-c.pos[3],2));

    //$ Check if they aren't in contact
    if (d >= Ri + Rj)
    {
        //$ Intersection is empty
        vol1 = vol2 = 0.;
        surf1 = surf2 = 0.;
    }
    else if(d >= fabs(Ri-Rj))
    {
        //$ Volume of the intersection is returned
        h = (pow(Rj,2)-pow((Ri-d),2))/(2.0*d);
        vol1= PI*pow(h,2)*(3*Ri-h)/3.0;
        surf1 = 2*PI*Ri*h;
        h = (pow(Ri,2)-pow((Rj-d),2))/(2.0*d);
        vol2= PI*pow(h,2)*(3*Rj-h)/3.0;
        surf2 = 2*PI*Rj*h;
    }
    //$ Check if j is completely absorbed by i
    else if (d < Ri - Rj)
    {
        //$ Volcal = VolJ
        vol1 = vol2 = c.Volume();
        surf1 = surf2 = c.Surface();
    }

    //$ Check if i is completely in j

    else // if (d < Rj - Ri)
    {
        //$ Volcal = Voli
        vol1 = vol2 = Volume();
        surf1 = surf2 = Surface();
    }

    return d;
}
