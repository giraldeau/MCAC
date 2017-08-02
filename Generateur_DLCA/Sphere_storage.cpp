
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
#include <utility>

void Sphere::SetPosition(const double newx, const double newy, const double newz)
{
    if (physicalmodel != nullptr)
    {
        *x = periodicPosition(newx,physicalmodel->L);
        *y = periodicPosition(newy,physicalmodel->L);
        *z = periodicPosition(newz,physicalmodel->L);
    }
    else
    {
        *x = newx;
        *y = newy;
        *z = newz;
    }
}
void Sphere::SetPosition(const double position[])
{
    SetPosition(position[1],position[2],position[3]);
}
void Sphere::SetPosition(const array<double, 4> position)
{
    SetPosition(position[1],position[2],position[3]);
}

Sphere::Sphere(void):
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(nullptr),
    AggLabel(0),
    SphereLabel(0),
    Storage(new array< vector<double>, 7>),
    external_storage(nullptr)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(1, 0.);
    Init();
}

Sphere::Sphere(PhysicalModel& _physicalmodel) : Sphere()
{

    physicalmodel = &_physicalmodel;
}

Sphere::Sphere(PhysicalModel& _physicalmodel, const double newx,const double newy,const double newz,const double newr) : Sphere(_physicalmodel)
{
    SetPosition(newx, newy, newz);
    *r = newr;
    UpdateVolAndSurf();
}

Sphere::Sphere(ListSphere& aggregat,const int id):
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(aggregat.physicalmodel),
    AggLabel(0),
    SphereLabel(id),
    Storage(aggregat.Storage),
    external_storage(&aggregat)
{
    Init();
    external_storage->setpointers();
}



Sphere::Sphere(PhysicalModel& _physicalmodel, const double newp[],const double newr) : Sphere(_physicalmodel,newp[1],newp[2],newp[3],newr){}
Sphere::Sphere(PhysicalModel& _physicalmodel, const array<double, 4> newp,const double newr) : Sphere(_physicalmodel,newp[1],newp[2],newp[3],newr){}

Sphere::Sphere(const Sphere& c) : Sphere(*(c.physicalmodel), *c.x,*c.y,*c.z,*c.r){}

void Sphere::add(void)
{
    for (int j=0;j<=6;j++)
        (*Storage)[j].push_back(0.);
    SphereLabel = int((*Storage)[0].size()-1);

    Init();
}

void Sphere::Init(void)
{
    setpointers();

    SetPosition(0, 0, 0);
    *r = 0;
    *volume = 0.;
    *surface = 0.;
}


Sphere::~Sphere(void)
{
    //setpointers();
    if(external_storage!=nullptr)
        external_storage->setpointers();
    else if (Storage!=nullptr)
    {
        for (int j=0;j<=6;j++)
            (*Storage)[j].erase((*Storage)[j].begin() + SphereLabel);
        delete Storage;
    }
}

 __attribute__((pure)) double& Sphere::operator[](const int i)
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

 __attribute__((pure)) double Sphere::Volume(void)
{
    return *volume;
}

 __attribute__((pure)) double Sphere::Surface(void)
{
    return *surface;
}

 __attribute__((pure)) double Sphere::Radius(void)
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
    SetPosition(newx, newy, newz);
    *r = newr;
    UpdateVolAndSurf();
}

void Sphere::Init(const double newp[],const double newr)
{
    Init(newp[1],newp[2],newp[3],newr);
}

void Sphere::Init(const array<double, 4> newp,const double newr)
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

void Sphere::Translate(const double xnew,const double ynew,const double znew)
{
    SetPosition(*x + xnew, *y + ynew, *z + znew);
}

void Sphere::Translate(const double trans[])
{
    SetPosition(*x + trans[1], *y + trans[2], *z + trans[3]);
}

void Sphere::Translate(const array<double, 4> trans)
{
    SetPosition(*x + trans[1], *y + trans[2], *z + trans[3]);
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

ListSphere::ListSphere(void):
    physicalmodel(nullptr),
    spheres(),
    index(),
    Storage(nullptr),
    external_storage(nullptr),
    N(0)
{}

ListSphere::ListSphere(PhysicalModel& _physicalmodel, const int _N) : ListSphere()
{
    Init(_physicalmodel,_N);
}


ListSphere::ListSphere(ListSphere& parent, int _index[]):
    physicalmodel(parent.physicalmodel),
    spheres(),
    index(),
    Storage(parent.Storage),
    external_storage(&parent),
    N(_index[0])
{
        index.assign(N+1, 0);
        spheres.assign(N, nullptr);

        index[0]=N;
        const int listSize = N;
        //#pragma omp for simd
        for (int i = 0; i < listSize; i++)
        {
            index[i+1] = _index[i+1];
            int iparent = external_storage->index[index[i+1]]-1;
            spheres[i] = external_storage->spheres[iparent];
        }
}



ListSphere::ListSphere(ListSphere& parent, int* _index[], const int start, const int end):
    physicalmodel(parent.physicalmodel),
    spheres(),
    index(),
    Storage(parent.Storage),
    external_storage(&parent),
    N(0)
{

    for(int i=start;i<=end;i++)
        N += _index[i][0];

    index.assign(N+1, 0);
    spheres.assign(N, nullptr);

    index[0]=N;
    int m=0;
    for(int i=start;i<=end;i++)
    {
    #pragma omp for simd
        for(int j=1;j<=_index[i][0];j++)
        {
            index[m+1] = _index[i][j];
            int iparent = external_storage->index[index[m+1]]-1;
            spheres[m] = external_storage->spheres[iparent];
            m++;
        }
}
}

void ListSphere::Init(PhysicalModel& _physicalmodel, const int _N)
{
    Destroy();

    physicalmodel=&_physicalmodel;
    Storage = new array< vector<double>, 7>;
    external_storage=nullptr;


    index.assign(_N+1, 0);
    spheres.assign(_N, nullptr);
    for (int j=0;j<=6;j++)
        (*Storage)[j].assign(_N, 0.);

    index[0]=_N;
    const int listSize = _N;
    for (N = 0; N < listSize; N++)
    {
        index[N+1] = N+1;
        spheres[N] = new Sphere(*this,N);
    }

}

void swap(ListSphere& first, ListSphere& second)
{
    using std::swap;
    swap(first.physicalmodel, second.physicalmodel);
    swap(first.spheres, second.spheres);
    swap(first.index, second.index);
    swap(first.Storage, second.Storage);
    swap(first.external_storage, second.external_storage);
    swap(first.N, second.N);
}



ListSphere& ListSphere::operator=(ListSphere other)
{
    swap(*this, other);

    return *this;
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
    if (external_storage==nullptr)
    {
        if (Storage!=nullptr)
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
        Storage = nullptr;
    }
    physicalmodel = nullptr;
}

 __attribute__((pure)) Sphere& ListSphere::operator[](const int i)
{
    return *spheres[i-1];
}

 __attribute__((pure)) int ListSphere::size() const
{
    return N;
}
