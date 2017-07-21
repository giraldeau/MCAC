
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

    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);

    external_storage=NULL;
    physicalmodel = NULL;
    AggLabel = 0;
    SphereLabel = 0;

    Init();
}
Sphere::Sphere(PhysicalModel& _physicalmodel)
{
    Storage = new array< vector<double>, 7>;

    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);

    external_storage=NULL;
    physicalmodel = &_physicalmodel;
    AggLabel = 0;
    SphereLabel = 0;

    Init();
}

Sphere::Sphere(ListSphere& aggregat,const int id)
{
    external_storage =&aggregat;

    Storage = aggregat.Storage;
    physicalmodel = aggregat.physicalmodel;
    AggLabel = 0;
    SphereLabel = id;

    Init();

    external_storage->setpointers();
}


Sphere::Sphere(PhysicalModel& _physicalmodel, const double newx,const double newy,const double newz,const double newr)
{
    Storage = new array< vector<double>, 7>;

    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);

    external_storage=NULL;
    physicalmodel = &_physicalmodel;
    AggLabel = 0;
    SphereLabel = 0;

    Init();

    *x = newx;
    *y = newy;
    *z = newz;
    *r = newr;
    UpdateVolAndSurf();
}

Sphere::Sphere(PhysicalModel& _physicalmodel, const double* newp,const double newr) : Sphere(_physicalmodel,newp[1],newp[2],newp[3],newr){}

Sphere::Sphere(Sphere& c) : Sphere(*(c.physicalmodel), *c.x,*c.y,*c.z,*c.r){}

void Sphere::add(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].push_back(0.);
    SphereLabel = (*Storage)[0].size()-1;

    Init();
}

void Sphere::Init(void)
{
    setpointers();

    *x = 0.;
    *y = 0;
    *z = 0;
    *r = 0;
    *volume = 0.;
    *surface = 0.;
}


Sphere::~Sphere(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].erase((*Storage)[j].begin() + SphereLabel);

    //setpointers();
    if(external_storage!=NULL)
        external_storage->setpointers();
    else if (Storage!=NULL)
        delete Storage;

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

const array<double, 4> Sphere::Position(void)
{
    array<double, 4> mypos;
    mypos[1]=*x;
    mypos[2]=*y;
    mypos[3]=*z;
    return mypos;
}

void Sphere::Init(const double newx,const double newy,const double newz,const double newr)
{
    *x = newx;
    *y = newy;
    *z = newz;
    *r = newr;
    UpdateVolAndSurf();
}

void Sphere::Init(const double* newp,const double newr)
{
    Init(newp[1],newp[2],newp[3],newr);
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

ListSphere::ListSphere(void)
{
        physicalmodel=NULL;

        external_storage=NULL;
        Storage = NULL;
}


ListSphere::ListSphere(PhysicalModel& _physicalmodel, const int _N)
{
    Init(_physicalmodel,_N);
}

void ListSphere::Init(PhysicalModel& _physicalmodel, const int _N)
{
    Destroy();

    physicalmodel=&_physicalmodel;
    Storage = new array< vector<double>, 7>;
    external_storage=NULL;


    index.assign(_N+1, 0);
    spheres.assign(_N, NULL);
    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(_N, 0.);

    index[0]=_N;
    for (N = 0; N < _N; N++)
    {
        index[N+1] = N+1;
        spheres[N] = new Sphere(*this,N);
    }

}

ListSphere::ListSphere(ListSphere& parent, int* _index)
{
        physicalmodel=parent.physicalmodel;
        external_storage=&parent;
        Storage = external_storage->Storage;

        N = _index[0];

        index.assign(N+1, 0);
        spheres.assign(N, NULL);

        index[0]=N;
        //#pragma omp for simd
        for (int i = 0; i < N; i++)
        {
            index[i+1] = _index[i+1];
            int iparent = external_storage->index[index[i+1]]-1;
            spheres[i] = external_storage->spheres[iparent];
        }
}



ListSphere::ListSphere(ListSphere& parent, int** _index, const int start, const int end)
{
    physicalmodel=parent.physicalmodel;
    external_storage=&parent;
    Storage = external_storage->Storage;

    N = 0;
    for(int i=start;i<=end;i++)
        N += _index[i][0];

    index.assign(N+1, 0);
    spheres.assign(N, NULL);

    index[0]=N;
    int m=0;
    for(int i=start;i<=end;i++)
        for(int j=1;j<=_index[i][0];j++)
        {
            index[m+1] = _index[i][j];
            int iparent = external_storage->index[index[m+1]]-1;
            spheres[m] = external_storage->spheres[iparent];
            m++;
        }
}


void ListSphere::setpointers()
{
    //#pragma omp for simd
    for (int i = 0; i < N; i++)
    {
        spheres[i]->setpointers();
    }
}

ListSphere::~ListSphere(void)
{
    Destroy();
}

void ListSphere::Destroy(void)
{
    if (external_storage==NULL)
    {
        if (Storage!=NULL)
        {
            int _N = N;
            for (N = _N; N > 0; N--)
            {
                delete spheres[N-1];
            }
            delete Storage;
        }
    }
    else
    {
        Storage = NULL;
    }
    physicalmodel = NULL;
}

Sphere& ListSphere::operator[](const int i)
{
    return *spheres[i-1];
}

int ListSphere::size() const
{
    return N;
}
