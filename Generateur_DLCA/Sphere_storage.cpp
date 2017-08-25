
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

/* Getters */

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

/* Setters */

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

/* Alias for different type of arguments*/

void Sphere::SetPosition(const double position[])
{
    SetPosition(position[1],position[2],position[3]);
}
void Sphere::SetPosition(const array<double, 4> position)
{
    SetPosition(position[1],position[2],position[3]);
}

/* Basic initializer */

void Sphere::InitVal(void)
{
    InitVal(0, 0, 0, 0.);
}

void Sphere::InitVal(const double newx,const double newy,const double newz,const double newr)
{
    setpointers();
    SetPosition(newx, newy, newz);
    *r = newr;
    UpdateVolAndSurf();
}

void Sphere::InitVal(const double newp[],const double newr)
{
    InitVal(newp[1],newp[2],newp[3],newr);
}

void Sphere::InitVal(const array<double, 4> newp,const double newr)
{
    InitVal(newp[1],newp[2],newp[3],newr);
}

/* format */

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



/* Storage specific */

void Sphere::setpointers(void)
{
    x = &(*this)[1];
    y = &(*this)[2];
    z = &(*this)[3];
    r = &(*this)[4];
    volume = &(*this)[5];
    surface =&(*this)[6];
}




/** Default constructor in local storage */
Sphere::Sphere(void):
    storage_elem<7,ListSphere>(),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(nullptr),
    AggLabel(0)
{
    InitVal();
}

Sphere::Sphere(PhysicalModel& _physicalmodel) : Sphere()
{
    physicalmodel = &_physicalmodel;
}

/** Constructor in local storage with initialization */
Sphere::Sphere(PhysicalModel& _physicalmodel, const double newx,const double newy,const double newz,const double newr) : Sphere(_physicalmodel)
{
    InitVal(newx, newy, newz, newr);
}
Sphere::Sphere(PhysicalModel& _physicalmodel, const array<double, 4> newp,const double newr) : Sphere(_physicalmodel,newp[1],newp[2],newp[3],newr){}
Sphere::Sphere(PhysicalModel& _physicalmodel, const double newp[],const double newr) : Sphere(_physicalmodel,newp[1],newp[2],newp[3],newr){}

/** Constructor with external storage */
Sphere::Sphere(ListSphere& aggregat,const int id):
    storage_elem<7,ListSphere>(aggregat,id),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(aggregat.physicalmodel),
    AggLabel(0)
{
    InitVal();
    external_storage->setpointers();
}

/** Copy constructor */
Sphere::Sphere(const Sphere& other) :
    storage_elem<7,ListSphere>(other),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(other.physicalmodel),
    AggLabel(other.AggLabel)
{
    InitVal(*other.x,*other.y,*other.z,*other.r);
}

/** Move constructor */
Sphere::Sphere (Sphere&& other) noexcept : /* noexcept needed to enable optimizations in containers */
    storage_elem<7,ListSphere>(other),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(other.physicalmodel),
    AggLabel(other.AggLabel)
{
    setpointers();
    other.x=nullptr;
    other.y=nullptr;
    other.z=nullptr;
    other.r=nullptr;
    other.volume=nullptr;
    other.surface=nullptr;
    other.AggLabel=0;
}

/** Destructor */
Sphere::~Sphere(void) noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    // everything is already taken care of in the parent class storage_elem
}

/** Copy assignment operator */
Sphere& Sphere::operator= (const Sphere& other)
{
    Sphere tmp(other);      // re-use copy-constructor
    *this = std::move(tmp); // re-use move-assignment
    return *this;
}

/** Move assignment operator */
Sphere& Sphere::operator= (Sphere&& other) noexcept
{
    *this = static_cast<Sphere&>(storage_elem<7,ListSphere>::operator=(other));
    physicalmodel = other.physicalmodel;
    AggLabel = other.AggLabel;
    setpointers();
    other.setpointers();
    other.AggLabel=0;
    return *this;
}










void ListSphere::Init(PhysicalModel& _physicalmodel, const int _N)
{
    physicalmodel=&_physicalmodel;
    storage_list<7,Sphere>::Init(_N,*this);
}


void ListSphere::DecreaseLabel(void)
{
    #pragma omp for simd
    for (int i = 0; i < size; i++)
    {
        list[i]->DecreaseLabel();
    }
}


void ListSphere::setpointers()
{
    //#pragma omp for simd
    for (int i = 0; i < size; i++)
    {
        list[i]->setpointers();
    }
}



/** Default constructor in local storage */
ListSphere::ListSphere(void):
    storage_list<7,Sphere>(),
    physicalmodel(nullptr)
{}

ListSphere::ListSphere(PhysicalModel& _physicalmodel, const int _N) :
    ListSphere()
{
    Init(_physicalmodel,_N);
}

/** Constructor with external storage */
ListSphere::ListSphere(ListSphere& parent,int _index[]):
    storage_list<7,Sphere>(parent, _index),
    physicalmodel(parent.physicalmodel)
{}

ListSphere::ListSphere(ListSphere& parent,int* _index[],const int start,const int end):
    storage_list<7,Sphere>(parent, _index, start, end),
    physicalmodel(parent.physicalmodel)
{}

/** Copy constructor */
ListSphere::ListSphere(const ListSphere& other):
    storage_list<7,Sphere>(other),
    physicalmodel(other.physicalmodel)
{}

/** Move constructor */
ListSphere::ListSphere (ListSphere&& other) noexcept: /* noexcept needed to enable optimizations in containers */
    storage_list<7,Sphere>(other),
    physicalmodel(other.physicalmodel)
{}

/** Destructor */
ListSphere::~ListSphere(void) noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{}

/** Copy assignment operator */
ListSphere& ListSphere::operator= (const ListSphere& other)
{
    ListSphere tmp(other);      // re-use copy-constructor
    *this = std::move(tmp);     // re-use move-assignment
    return *this;
}

/** Move assignment operator */
ListSphere& ListSphere::operator= (ListSphere&& other) noexcept
{
    std::swap(static_cast<storage_list<7,Sphere>&>(*this),static_cast<storage_list<7,Sphere>&>(other));
    physicalmodel = other.physicalmodel;

    return *this;
}
