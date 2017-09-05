
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




/* #############################################################################################################
 * #################################                                       #####################################
 * #################################             SPHERE                    #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/
void Sphere::SetPosition(const double newx, const double newy, const double newz) noexcept
{
    *x = periodicPosition(newx,physicalmodel->L);
    *y = periodicPosition(newy,physicalmodel->L);
    *z = periodicPosition(newz,physicalmodel->L);
}

/* #############################################################################################################
 * ################################# Distance between a sphere and a point #####################################
 * #############################################################################################################*/
__attribute((pure)) double Sphere::Distance(const Sphere& c) const noexcept
{
    return Distance(*c.x,*c.y,*c.z);
}

__attribute((pure)) double Sphere::Distance(const array<double, 3> point) const noexcept
{
    return Distance(point[0],point[1],point[2]);
}
__attribute((pure)) double Sphere::Distance2(const Sphere& c) const noexcept
{
    return Distance2(*c.x,*c.y,*c.z);
}

__attribute((pure)) double Sphere::Distance2(const array<double, 3> point) const noexcept
{
    return Distance2(point[0],point[1],point[2]);
}

__attribute__((pure)) double Sphere::Distance2(const double otherx, const double othery, const double otherz) const noexcept
{
   double dx(periodicDistance((*x-otherx),physicalmodel->L));
   double dy(periodicDistance((*y-othery),physicalmodel->L));
   double dz(periodicDistance((*z-otherz),physicalmodel->L));
   return POW2(dx)+POW2(dy)+POW2(dz);
}

 __attribute__((pure)) double Sphere::Distance(const double otherx, const double othery, const double otherz) const noexcept
{
    return sqrt(Distance2(otherx,othery,otherz));
}

/* #############################################################################################################
 * ##################################### Volume and surface of a sphere ########################################
 * #############################################################################################################*/

const double PI = atan(1.0)*4;
const double facvol = 4*PI/3;
const double facsurf = 4*PI;
void Sphere::UpdateVolAndSurf(void) noexcept
{
    if (AggLabel > -1)
    {
        *volume = facvol*POW3(*r);
        *surface =  facsurf*POW2(*r);
    }
}

 /* #############################################################################################################
  * ########################################## Distance before collision ########################################
  * #############################################################################################################*/
  __attribute__((pure)) bool Sphere::Contact(const Sphere& c) const noexcept
 {
     //$ Compute signed distance for contact between two spheres
     double distance = Distance2(c);

     //$ Compute minimum distance for contact
     double dist_contact = POW2(*r + *c.r);

     return (distance <= dist_contact);
 }

  __attribute__((pure)) double Sphere::Collision(const Sphere& c,const array<double,3> vd) const
  {
      /*
       * Denoting
       * V the unitary displacement vector
       * D the vector from the moving sphere to the other at start
       * C the vector from the moving sphere to the other at collision
       * x the distance we are looking for
       *
       * We have a triangle so :
       *  --> --> -->
       * x V + C - D = 0
       *
       * i.e
       * -->   -->  -->
       *  D - x V  = C
       *
       * taking the norm
       *              --> -->
       * d² + x² - 2 x V . D = c²
       *
       * i.e x is solution of
       *         --> -->
       * x² - 2 x V . D  + d² - c² = 0
       *
      */

      //$ Compute signed distance for contact between two spheres
      double distance = Distance2(c);

      //$ Compute minimum distance for contact
      double dist_contact = POW2(*r + *c.r);

      //$ Computing distance before contact
      double dx = periodicDistance((*c.x-*x),physicalmodel->L);
      double dy = periodicDistance((*c.y-*y),physicalmodel->L);
      double dz = periodicDistance((*c.z-*z),physicalmodel->L);
      double VD = -2*(dx*vd[0] + dy*vd[1] + dz*vd[2]);
      double DC = distance - dist_contact;
      double DELTA = VD*VD - 4*DC;

      if (DELTA >= 0)
      {
          DELTA = sqrt(DELTA);
          if (DELTA <= -VD)
              return 0.5*(-VD-DELTA);
          else
              return 0.5*(-VD+DELTA);
      }
      else
          return -1.;
}

/* #############################################################################################################
 * ######################## Surface and volume of the intersection of two sphere ###############################
 * #############################################################################################################*/
double Sphere::Intersection(const Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) const
{
    double dist = Distance(c);
    return Intersection(c,dist,vol1,vol2,surf1,surf2);
}

double Sphere::Intersection(const Sphere& c, const double dist,double& vol1, double& vol2, double& surf1, double& surf2 ) const
{

    vol1 = vol2 = 0.;
    surf1 = surf2 = 0.;

    double Ri = *r;
    double Rj = *c.r;

    //$ Check if they are in contact
    if (dist < Ri + Rj)
    {
        if(dist >= fabs(Ri - Rj))
        {
            //$ Volume of the intersection is returned
            double h1 = (POW2(Rj)-POW2((Ri-dist)))/(2.0*dist);
            double h2 = (POW2(Ri)-POW2((Rj-dist)))/(2.0*dist);
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
    return dist;
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
