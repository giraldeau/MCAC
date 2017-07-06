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

#define ROWMAJOR

using namespace std;

const double PI = atan(1.0)*4;
const double facvol = 4*PI/3;
const double facsurf = 4*PI;

Sphere::Sphere(void)
{
    arr[0] = new double[7];
    for (int i=1;i<=6;i++)
        arr[i]=&arr[0][i];

    x = arr[1];
    y = arr[2];
    z = arr[3];
    r = arr[4];
    volume = arr[5];
    surface = arr[6];
    AggLabel = 0;
    SphereLabel = 0;
    Update(0, 0, 0, 0);

    external_storage=false;
}

Sphere::Sphere(double** ext_arr,const int i)
{
    for (int j=0;j<=6;j++)
#ifdef ROWMAJOR
        arr[j]=&ext_arr[j][i];
#else
        arr[j]=&ext_arr[i][j];
#endif

    x = arr[1];
    y = arr[2];
    z = arr[3];
    r = arr[4];
    volume = arr[5];
    surface = arr[6];
    AggLabel = 0;
    SphereLabel = i;
    Update(0, 0, 0, 0);

    external_storage =true;
}

Sphere::~Sphere(void)
{
    if(!external_storage)
        if (arr[0]!=NULL) delete[] arr[0];

    for (int j=1;j<=6;j++)
        arr[j]=NULL;

    x = NULL;
    y = NULL;
    z = NULL;
    r = NULL;
    volume = NULL;
    surface = NULL;
}

double Sphere::Volume(void) const
{
    return *volume;
}

double Sphere::Surface(void) const
{
    return *surface;
}

double Sphere::Radius(void) const
{
    return *r;
}

const double* Sphere::Position(void) const
{

    double* mypos= new double[4];
    mypos[1]=*x;
    mypos[2]=*y;
    mypos[3]=*z;
    return mypos;

    //return arr[0];
}

void Sphere::Update(const double newx,const double newy,const double newz,const double newr)
{
    *x = newx;
    *y = newy;
    *z = newz;
    *r = newr;
    UpdateVolAndSurf();
}

void Sphere::Update(const double* newp,const double newr)
{
    Update(newp[1],newp[2],newp[3],newr);
}

void Sphere::Update(const Sphere& c)
{
    Update(*c.x,*c.y,*c.z,*c.r);
}


void Sphere::SetLabel(const int value)
{
    AggLabel = value;
}

void Sphere::DecreaseLabel(void)
{
    AggLabel--;
}

void Sphere::Translate(const double* trans)
{
    *x += trans[1];
    *y += trans[2];
    *z += trans[3];
}

string Sphere::str(const double coef) const
{
    stringstream res;
    res
        << setw(5) << AggLabel << "\t"
        << setprecision(5) << setw(11) << fixed << *r*coef << "\t"
        << setprecision(5) << setw(11) << fixed << *x*coef << "\t"
        << setprecision(5) << setw(11) << fixed << *y*coef << "\t"
        << setprecision(5) << setw(11) << fixed << *z*coef;
    return res.str();
}

void Sphere::Aff(const double coef) const
{  
    cout << str(coef) << endl;
}

void Sphere::UpdateVolAndSurf(void)
{
    if (AggLabel > 0)
    {
        *volume = facvol*pow(*r, 3);
        *surface =  facsurf*pow(*r, 2);
    }
}

double Sphere::Distance(const Sphere& c) const
{
    return Distance(*c.x,*c.y,*c.z);
}

double Sphere::Distance(const double* point) const
{
    return Distance(point[1],point[2],point[3]);
}

double Sphere::Distance(const double otherx, const double othery, const double otherz) const
{
    return sqrt(pow(*x-otherx,2)+pow(*y-othery,2)+pow(*z-otherz,2));
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

//Calcul du volume de la calotte sphérique de la sphère courante de rayon Ri due à la surestimation de la sphère c de rayon Rj
double Sphere::Intersection(const Sphere& c,double& vol1, double& vol2, double& surf1, double& surf2 ) const
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

void SphereList::Init(const int _N,PhysicalModel& _physicalmodel)
{
    N=_N;
    spheres = new Sphere*[N+1];
    physicalmodel=&_physicalmodel;

#ifdef ROWMAJOR
    array = new double*[7];
    for (int i = 1; i <= 6; i++)
        array[i] = new double[N+1];
#else
    array = new double*[N+1];
    for (int i = 1; i <= N; i++)
        array[i] = new double[7];
#endif

    for (int i = 1; i <= _N; i++)
        spheres[i] = new Sphere(array,i);

    external_storage = false;
}

SphereList::~SphereList(void)
{
    if(!external_storage)
        if (spheres!=NULL) delete[] spheres;
        if (array!=NULL) delete[] array;
    physicalmodel = NULL;
}

Sphere& SphereList::operator[](const int i)
{
    return *spheres[i];
}

//####################################### Croissance de surface des particules primaires ########################################
void SphereList::CroissanceSurface(double dt)
{
    #pragma omp for simd
    for (int i = 1; i <= N; i++)
    {
#ifdef ROWMAJOR
        array[4][i] = physicalmodel->Grow(array[4][i], dt);
        array[5][i] = facvol*pow(array[4][i], 3);
        array[6][i] = facsurf*pow(array[4][i], 2);
#else
        array[i][4] = physicalmodel->Grow(array[i][4], dt);
        array[i][5] = facvol*pow(array[i][4], 3);
        array[i][6] = facsurf*pow(array[i][4], 2);
#endif
    }
}
