#include "Sphere.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace std;

const double PI = atan(1.0)*4;

Sphere::Sphere(void)
{
    pos = new double[4];
    Label = 1;
    Update(0, 0, 0, 1);
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
}

void Sphere::Update(const double* newp,const double newr)
{
    // newp: position double[4] pour compatibilité....
    pos[1] = newp[1];
    pos[2] = newp[2];
    pos[3] = newp[3];
    r = newr;
}


void Sphere::SetLabel(int value)
{
    Label = value;
}

void Sphere::Translate(const double* trans)
{
    for (int i = 1; i <= 3; i++) pos[i] = pos[i] + trans[i];
}

void Sphere::Translate(const double trans)
{
    for (int i = 1; i <= 3; i++) pos[i] = pos[i] + trans;
}

void Sphere::Aff(const double coef) const
{
    printf("%8.3f %8.3f %8.3f %8.3f\n", pos[1]*coef, pos[2]*coef, pos[3]*coef, r*coef);
}
double Sphere::Distance(const Sphere& c) const
{
    return sqrt(pow(pos[1]-c.pos[1],2)+pow(pos[2]-c.pos[2],2)+pow(pos[3]-c.pos[3],2));
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


double Sphere::Collision_opt(const Sphere& c,const  double* vd,const double distmax,double& distance) const
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
    double dist_proj,min_dist,dist_contact;
    double dx, dy, dz;
    double dist;

    dx = c.pos[1] - pos[1];
    dy = c.pos[2] - pos[2];
    dz = c.pos[3] - pos[3];

    dist_contact = pow(r + c.r,2);
    distance = dx*dx+dy*dy+dz*dz;

    if (distance <= dist_contact)
    {
        //$ There is already contact between these two spheres
        return -1;
    }
    dist_proj = dx*vd[1]+dy*vd[2]+dz*vd[3];
    min_dist = pow(dist_proj*vd[1] - dx,2) + pow(dist_proj*vd[2] - dy,2) + pow(dist_proj*vd[3] - dz,2);

    if (min_dist > dist_contact)
    {
        //$ There is not contact possible between these two spheres
        return  1;
    }
    dist = dist_proj-sqrt(dist_contact-min_dist);
    if (dist > distmax || dist < 0)
    {
        //$ The contact would be too late, or is behind
        return  1;
    }
    return dist;
}

//Calcul du volume de la calotte sphérique de la sphère courante de rayon Ri due à la surestimation de la sphère c de rayon Rj
double Sphere::VolumeCalotteij(const Sphere& c) const
{
    double Volcal, d, h;
    double Ri, Rj;

    Ri = r;
    Rj = c.r;
    //$ Determination of the distance between the center of the 2 aggregates
    d = sqrt(pow(pos[1]-c.pos[1],2) + pow(pos[2]-c.pos[2],2) + pow(pos[3]-c.pos[3],2));
    //$ Check if they aren't in contact
    if (d >= Ri + Rj)
    {
        //$ Volcal = 0

        return 0.0;
    }


    //$ Check if j is completely absorbed by i
    if ((0 <= d) && (d < Ri - Rj))
    {
        //$ Volcal = VolJ
        return 4.0*PI*pow(Rj,3)/3.0;

    }

    //$ Check if i is completely in j

    if ((0 <= d) && (d < Rj - Ri))
    {
        //$ Volcal = Voli
        return 4.0*PI*pow(Ri,3)/3.0;
    }

    //$ Volume of the intersection is returned
        h = (pow(Rj,2)-pow((Ri-d),2))/(2.0*d);

    return PI*pow(h,2)*(3*Ri-h)/3.0;

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
        vol1 = vol2 = 4.0*PI*pow(Rj,3)/3.0;
        surf1 = surf2 = 4.0*PI*pow(Rj,2);
    }

    //$ Check if i is completely in j

    else // if (d < Rj - Ri)
    {
        //$ Volcal = Voli
        vol1 = vol2 = 4.0*PI*pow(Ri,3)/3.0;
        surf1 = surf2 = 4.0*PI*pow(Ri,2);
    }

    return d;
}

//Calcul de la surface de la calotte sphérique de la sphère courante de rayon Ri due à la surestimation de la sphère c de rayon Rj
double Sphere::SurfaceCalotteij(const Sphere& c) const // Works exactly the same way as VolumeCalotte, the only things that changes is the formulas
{
    double Surfcal, d, h;
    double Ri, Rj;

    Ri = r;
    Rj = c.r;

    d = sqrt(pow(pos[1]-c.pos[1], 2)+pow(pos[2]-c.pos[2], 2)+pow(pos[3]-c.pos[3], 2));

    if (d >= Ri + Rj)
        Surfcal = 0.0;
    else if ((0 <= d) && (d < Ri - Rj))
        Surfcal = 4.0*PI*pow(Rj, 2);
    else if ((0 <= d) && (d < Rj - Ri))
        Surfcal = 4.0*PI*pow(Ri, 2);
    else
    {
        h = (pow(Rj,2)-pow((Ri-d),2))/(2.0*d);
        Surfcal = 2*PI*Ri*h;
    }

    return Surfcal;
}
