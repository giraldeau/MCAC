
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



#include "Sphere.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))

using namespace std;
namespace DLCA{


/* #############################################################################################################
 * #################################                                       #####################################
 * #################################             SPHERE                    #####################################
 * #################################                                       #####################################
 * #############################################################################################################*/
void Sphere::SetPosition(const double newx, const double newy, const double newz) noexcept
{
    /*
    *x = periodicPosition(newx,physicalmodel->L);
    *y = periodicPosition(newy,physicalmodel->L);
    *z = periodicPosition(newz,physicalmodel->L);
    */
    *x = newx;
    *y = newy;
    *z = newz;
}

/* #############################################################################################################
 * ################################# Distance between a sphere and a point #####################################
 * #############################################################################################################*/
__attribute((pure)) double Sphere::RelativeDistance(const Sphere& c) const noexcept
{
    return RelativeDistance(*c.rx,*c.ry,*c.rz);
}

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

__attribute((pure)) double Sphere::RelativeDistance2(const Sphere& c) const noexcept
{
    return Distance2(*c.rx,*c.ry,*c.rz);
}
__attribute__((pure)) double Sphere::RelativeDistance2(const double otherx, const double othery, const double otherz) const noexcept
{
   return POW2(*rx-otherx)+POW2(*ry-othery)+POW2(*rz-otherz);
}

__attribute__((pure)) double Sphere::RelativeDistance(const double otherx, const double othery, const double otherz) const noexcept
{
   return sqrt(RelativeDistance2(otherx,othery,otherz));
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
void Sphere::UpdateVolAndSurf() noexcept
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

     // 1e-28 is for rounding error (1e-14 ^ 2)
     return (distance - dist_contact <= 1e-28);
 }

  __attribute__((pure)) vector<double> Sphere::Collisions(const ListSphere& list,const array<double,3> vectordir) const
{
      array < array < double, 3>, 3> RotMat = GetRotMat(vectordir);

      vector<double> distToCollision;
      for (const Sphere* c : list)
      {
          pair<bool,double> suspect = CollisionR(*c,RotMat);
          if (suspect.first)
          {
              distToCollision.push_back(suspect.second);
          }
      }
      return distToCollision;
}

__attribute__((pure)) pair<bool,double> Sphere::Collision(const Sphere& c,const array<double,3> vectordir) const
{
      array < array < double, 3>, 3> RotMat = GetRotMat(vectordir);
      return CollisionR(c,RotMat);
}


  __attribute__((pure)) pair<bool,double> Sphere::CollisionR(const Sphere& c,const array < array < double, 3>, 3> RotMat) const
  {
      /*
       * We use a change of axis system
       * the center is placed on the mobil sphere
       * the x axis is chosen to be the movement vector
       *
       * Thus we have two rotation to perform : one around z and one around y
       *
       * Then a collision is easy to detect :
       * The distance to the x axis must be less than the sum of the radius
       *
       * Moreover we are only interested in collision that can happend in the future,
       * wich is also easy to detect (x coordinate must be positive)
       *
       * Finally we consider 27 possible version of the other sphere in order to take into account
       * any periodicity effect.
      */

      double dist_contact = POW2(*r + *c.r);
      double dist = Distance2(c);
      double minval=0.;
      bool collision = true;

      if(dist>dist_contact)
      {
          double L = physicalmodel->L;

          collision = false;
          minval=10*L;

          array < double, 3> pos{};
          for(size_t l = 0; l < 3; ++l)
          {
              pos[l] = RotMat[l][0] * (*c.x - *x)
                     + RotMat[l][1] * (*c.y - *y)
                     + RotMat[l][2] * (*c.z - *z);
          }
          array < double, 3> perx{{L*RotMat[0][0],
                                  L*RotMat[1][0],
                                  L*RotMat[2][0],}};
          array < double, 3> pery{{L*RotMat[0][1],
                                  L*RotMat[1][1],
                                  L*RotMat[2][1],}};
          array < double, 3> perz{{L*RotMat[0][2],
                                  L*RotMat[1][2],
                                  L*RotMat[2][2],}};

          for (int i=-1;i<=1;i++)
          {
              for (int j=-1;j<=1;j++)
              {
                  for (int k=-1;k<=1;k++)
                  {
                      array < double, 3> tmp{};

                      tmp[0] = pos[0] + i*perx[0] + j*pery[0] + k*perz[0];

                      // in the future
                      if (tmp[0] <0.)
                      {
                        continue;
                      }

                      for(size_t l = 1; l < 3; ++l)
                      {
                          tmp[l] = pos[l]
                                 + i*perx[l]
                                 + j*pery[l]
                                 + k*perz[l];
                      }

                      double dist1 = POW2(tmp[1]) + POW2(tmp[2]);

                      // collision is possible
                      if (dist1 > dist_contact)
                      {
                        continue;
                      }

                      collision =true;

                      double sol = dist_contact - dist1;
                      sol = tmp[0] - sqrt(sol);

                      minval = MIN(minval,sol);
                  }
              }
          }
      }
      pair<bool,double> result = {collision,minval};
      return result;
}


/* #############################################################################################################
 * ######################## Surface and volume of the intersection of two sphere ###############################
 * #############################################################################################################*/
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


array < array < double, 3>, 3> GetRotMat(const array<double,3> vectordir)
{
    double anglez = -atan2(vectordir[1],vectordir[0]);
    array < array < double, 3>, 3> rotz{};
    rotz[0] = {{cos(anglez), -sin(anglez), 0}};
    rotz[1] = {{sin(anglez),  cos(anglez), 0}};
    rotz[2] = {{          0,            0, 1}};

    array < double, 3> tmp{};
    for(size_t i = 0; i < 3; ++i)
    {
        tmp[i] = 0;
        for(size_t j = 0; j < 3; ++j)
        {
            tmp[i] += rotz[i][j] * vectordir[j];
        }
    }

    double angley = atan2(tmp[2],tmp[0]);
    array < array < double, 3>, 3> roty{};
    roty[0] = {{ cos(angley), 0, sin(angley)}};
    roty[1] = {{           0, 1, 0          }};
    roty[2] = {{-sin(angley), 0, cos(angley)}};

    array < array < double, 3>, 3> matrot{};
    for(size_t i = 0; i < 3; ++i)
    {
        for(size_t j = 0; j < 3; ++j)
        {
            matrot[i][j]=0;
            for(size_t k = 0; k < 3; ++k)
            {
                matrot[i][j] += roty[i][k] * rotz[k][j];
            }
        }
    }
    return matrot;
}

}  // namespace DLCA

