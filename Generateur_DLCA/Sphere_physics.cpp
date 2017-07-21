
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


/* #############################################################################################################
 * #################################                                       #####################################
 * #################################             SPHERE                    #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/


/* #############################################################################################################
 * ################################# Distance between a sphere and a point #####################################
 * #############################################################################################################*/
double Sphere::Distance(Sphere& c)
{
    return Distance(*c.x,*c.y,*c.z);
}

double Sphere::Distance(const double* point)
{
    return Distance(point[1],point[2],point[3]);
}
double Sphere::Distance(const std::array<double, 4> point)
{
    return Distance(point[1],point[2],point[3]);
}

double Sphere::Distance(const double otherx, const double othery, const double otherz)
{
    return sqrt(pow(*x-otherx,2)+pow(*y-othery,2)+pow(*z-otherz,2));
}

/* #############################################################################################################
 * ##################################### Volume and surface of a sphere ########################################
 * #############################################################################################################*/
void Sphere::UpdateVolAndSurf(void)
{
    if (AggLabel > 0)
    {
        *volume = facvol*pow(*r, 3);
        *surface =  facsurf*pow(*r, 2);
    }
}

/* #############################################################################################################
 * ########################################## Distance before collision ########################################
 * #############################################################################################################*/
double Sphere::Collision(Sphere& c,const  double* vd,const double distmax)
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
    double distance;

    dx = *c.x - *x;
    dy = *c.y - *y;
    dz = *c.z - *z;

    dist_contact = pow(*r + *c.r,2);

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


/* #############################################################################################################
 * ######################## Surface and volume of the intersection of two sphere ###############################
 * #############################################################################################################*/
double Sphere::Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 )
{
    double d, h;
    double Ri, Rj;

    vol1 = vol2 = 0.;
    surf1 = surf2 = 0.;

    Ri = *r;
    Rj = *c.r;


    //$ Determination of the distance between the center of the 2 aggregates
    d = Distance(c);

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
        vol1 = vol2 = *c.volume;
        surf1 = surf2 = *c.surface;
    }

    //$ Check if i is completely in j

    else // if (d < Rj - Ri)
    {
        //$ Volcal = Voli
        vol1 = vol2 = *volume;
        surf1 = surf2 = *surface;
    }

    return d;
}



/* #############################################################################################################
 * ############################################# Sphere growing ################################################
 * #############################################################################################################*/
void Sphere::CroissanceSurface(const double dt)
{
    double newR = physicalmodel->Grow(*r, dt);
    double newR2=newR*newR;
    double newR3=newR2*newR;
    *r = newR;
    *volume = facvol*newR3;
    *surface = facsurf*newR2;
}



/* #############################################################################################################
 * #################################                                       #####################################
 * #################################              AGREGATE                 #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/



/* #############################################################################################################
 * ########################################### Grow all spheres ################################################
 * #############################################################################################################*/
void ListSphere::CroissanceSurface(const double dt)
{
    const int listSize = N;
    #pragma omp for simd
    for (int i = 0; i < listSize; i++)
    {
        spheres[i]->CroissanceSurface(dt);
    }
}

