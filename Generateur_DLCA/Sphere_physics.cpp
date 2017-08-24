
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
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

const double PI = atan(1.0)*4;
const double facvol = 4*PI/3;
const double facsurf = 4*PI;

__attribute((const)) double periodicDistance(const double x, const double dim)
{
        double hdim(0.5*dim);
        if (x<-hdim)
            return x+dim;
        if (x>hdim)
            return x-dim;
    return x;
}


__attribute((const)) double periodicPosition(const double x, const double dim)
{
        if (x<0)
            return x+dim;
        if (x>dim)
            return x-dim;

    return x;

}


/* #############################################################################################################
 * #################################                                       #####################################
 * #################################             SPHERE                    #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/


/* #############################################################################################################
 * ################################# Distance between a sphere and a point #####################################
 * #############################################################################################################*/
__attribute((pure)) double Sphere::Distance(Sphere& c)
{
    return Distance(*c.x,*c.y,*c.z);
}

__attribute((pure)) double Sphere::Distance(const double point[])
{
    return Distance(point[1],point[2],point[3]);
}
__attribute((pure)) double Sphere::Distance(const array<double, 4> point)
{
    return Distance(point[1],point[2],point[3]);
}
__attribute((pure)) double Sphere::Distance2(Sphere& c)
{
    return Distance2(*c.x,*c.y,*c.z);
}

__attribute((pure)) double Sphere::Distance2(const double point[])
{
    return Distance2(point[1],point[2],point[3]);
}
__attribute((pure)) double Sphere::Distance2(const array<double, 4> point)
{
    return Distance2(point[1],point[2],point[3]);
}

__attribute__((pure)) double Sphere::Distance2(const double otherx, const double othery, const double otherz)
{
    if (physicalmodel != nullptr)
    {
       double dx(periodicDistance((*x-otherx),physicalmodel->L));
       double dy(periodicDistance((*y-othery),physicalmodel->L));
       double dz(periodicDistance((*z-otherz),physicalmodel->L));
       return POW2(dx)+POW2(dy)+POW2(dz);
    }
    else
    {
         return POW2(*x-otherx)+POW2(*y-othery)+POW2(*z-otherz);
    }
}

 __attribute__((pure)) double Sphere::Distance(const double otherx, const double othery, const double otherz)
{
    return sqrt(Distance2(otherx,othery,otherz));
}

/* #############################################################################################################
 * ##################################### Volume and surface of a sphere ########################################
 * #############################################################################################################*/
void Sphere::UpdateVolAndSurf(void)
{
    if (AggLabel > 0)
    {
        *volume = facvol*POW3(*r);
        *surface =  facsurf*POW2(*r);
    }
}

/* #############################################################################################################
 * ########################################## Distance before collision ########################################
 * #############################################################################################################*/
 __attribute__((pure)) double Sphere::Collision(Sphere& c,const array<double,4> vd,const double distmax)
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
    double DELTA;
    double K, K1, K2;
    double dist,dist_contact;
    double distance;

    double dx = periodicDistance((*c.x-*x),physicalmodel->L);
    double dy = periodicDistance((*c.y-*y),physicalmodel->L);
    double dz = periodicDistance((*c.z-*z),physicalmodel->L);

    dist_contact = POW2(*r + *c.r);

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
            //if (dist > distmax) dist = 1;
        }
    }
    return dist;
}


/* #############################################################################################################
 * ######################## Surface and volume of the intersection of two sphere ###############################
 * #############################################################################################################*/
double Sphere::Intersection(Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 )
{

    vol1 = vol2 = 0.;
    surf1 = surf2 = 0.;

    double Ri = *r;
    double Rj = *c.r;

    //$ Determination of the distance between the center of the 2 spheres
    double d = Distance(c);

    //$ Check if they are in contact
    if (d < Ri + Rj)
    {
        if(d >= fabs(Ri - Rj))
        {
            //$ Volume of the intersection is returned
            double h1 = (POW2(Rj)-POW2((Ri-d)))/(2.0*d);
            double h2 = (POW2(Ri)-POW2((Rj-d)))/(2.0*d);
            vol1= PI*POW2(h1)*(3*Ri-h1)/3.0;
            vol2= PI*POW2(h2)*(3*Rj-h2)/3.0;
            surf1 = 2*PI*Ri*h1;
            surf2 = 2*PI*Rj*h2;
        }
        //$ Check if one is completely absorbed by the other
        else if (Ri < Rj)
        {
            //$ Volcal = VolJ
            vol1 = *volume;
            surf1 = *surface;
        }
        else // if (Rj < Ri)
        {
            //$ Volcal = Voli
            vol2 = *c.volume;
            surf2 = *c.surface;
        }
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
    const int listSize = size;
    #pragma omp for simd
    for (int i = 0; i < listSize; i++)
    {
        list[i]->CroissanceSurface(dt);
    }
}

