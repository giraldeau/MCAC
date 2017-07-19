
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

Sphere::Sphere(void)
{
    Storage = new array< vector<double>, 7>;
    external_storage=NULL;
    physicalmodel = NULL;
    AggLabel = 0;

    add();
}
Sphere::Sphere(PhysicalModel& _physicalmodel)
{
    Storage = new array< vector<double>, 7>;
    external_storage=NULL;
    physicalmodel = &_physicalmodel;
    AggLabel = 0;

    add();
}

Sphere::Sphere(Aggregat* aggregat)
{
    external_storage =aggregat;

    Storage = (*external_storage).Storage;
    physicalmodel = (*external_storage).physicalmodel;
    AggLabel = 0;

    add();

    external_storage->setpointers();
}

void Sphere::add(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].push_back(0.);
    SphereLabel = (*Storage)[0].size()-1;

    setpointers();

    Update(0, 0, 0, 0);
}

Sphere::~Sphere(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].erase((*Storage)[j].begin() + SphereLabel);

    setpointers();
    if(external_storage!=NULL)
        external_storage->setpointers();

    if(external_storage==NULL)
        if (Storage!=NULL) delete Storage;

    external_storage=NULL;
    SphereLabel = 0;
    AggLabel = 0;

    x = NULL;
    y = NULL;
    z = NULL;
    r = NULL;
    volume = NULL;
    surface = NULL;
    physicalmodel = NULL;

}

double& Sphere::operator[](const int i)
{
    return (*Storage)[i][SphereLabel];
}

void Sphere::setpointers(void)
{
    x = &(*this)[1];
    y = &(*this)[2];
    z = &(*this)[3];
    r = &(*this)[4];
    volume = &(*this)[5];
    surface =&(*this)[6];
}

double Sphere::Volume(void)
{
    return *volume;
}

double Sphere::Surface(void)
{
    return *surface;
}

double Sphere::Radius(void)
{
    return *r;
}

const double* Sphere::Position(void)
{
    double* mypos= new double[4];
    mypos[1]=*x;
    mypos[2]=*y;
    mypos[3]=*z;
    return mypos;
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

void Sphere::Update(Sphere& c)
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

string Sphere::str(const double coef)
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

void Sphere::Aff(const double coef)
{  
    cout << str(coef) << endl;
}

double Sphere::Distance(Sphere& c)
{
    return Distance(*c.x,*c.y,*c.z);
}

double Sphere::Distance(const double* point)
{
    return Distance(point[1],point[2],point[3]);
}

double Sphere::Distance(const double otherx, const double othery, const double otherz)
{
    return sqrt(pow(*x-otherx,2)+pow(*y-othery,2)+pow(*z-otherz,2));
}


Aggregat::Aggregat(void)
{
        spheres = NULL;
        physicalmodel=NULL;

        external_storage=NULL;
        index = NULL;
        Storage = NULL;
}


Aggregat::Aggregat(const Aggregat* parent, int* _index)
{
        spheres = new Sphere*[_index[0]];
        physicalmodel=parent->physicalmodel;

        external_storage=parent;
        index = _index;
        Storage = external_storage->Storage;

        N = index[0];

        //#pragma omp for simd
        for (int i = 0; i < N; i++)
        {
            int iparent = external_storage->index[index[i+1]]-1;
            spheres[i] = external_storage->spheres[iparent];
        }
}


void Aggregat::Init(const int _N,PhysicalModel& _physicalmodel)
{
    spheres = new Sphere*[_N];
    physicalmodel=&_physicalmodel;

    Storage = new array< vector<double>, 7>;

    index = new int[_N+1];
    index[0]=_N;
    for (N = 0; N < _N; N++)
    {
        index[N+1] = N+1;
        spheres[N] = new Sphere(this);
    }

}


void Aggregat::setpointers()
{
    //#pragma omp for simd
    for (int i = 0; i < N-1; i++)
    {
        spheres[i]->setpointers();
    }
}

Aggregat::~Aggregat(void)
{
    if (external_storage==NULL)
    {
        if (spheres!=NULL)
        {
            delete[] spheres;
        }
        if (Storage!=NULL)
            delete Storage;
    }
    else
    {
        spheres = NULL;
        Storage = NULL;
    }
    physicalmodel = NULL;
}

Sphere& Aggregat::operator[](const int i)
{
    return *spheres[i-1];
}

int Aggregat::size() const
{
    return N;
}


//################################################## Recherche de sphères #############################################################################

Aggregat Aggregat::extract(const int id, int** AggLabels) const
{
    return Aggregat(this,AggLabels[id]);
}

Aggregat Aggregat::extractplus(const int id, int** AggLabels,const int NAgg) const
{

    int NSphereInAggrerat=0;

    for(int i=id;i<=NAgg;i++)
        NSphereInAggrerat += AggLabels[i][0];

    int* AggLabelPlus = new int[NSphereInAggrerat+1];

    AggLabelPlus[0] = NSphereInAggrerat;

    int m =0;
    for(int i=id;i<=NAgg;i++)
        for(int j=1;j<=AggLabels[i][0];j++)
        {
            m++;
            AggLabelPlus[m] = AggLabels[i][j];
        }

    Aggregat ret = Aggregat(this,AggLabelPlus);

    delete[] AggLabelPlus;

    return ret;

}
