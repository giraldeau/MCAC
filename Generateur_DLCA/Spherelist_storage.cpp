
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

void ListSphere::Init(PhysicalModel& _physicalmodel, const int _N)
{
    physicalmodel=&_physicalmodel;
    storage_list<7,Sphere>::Init(_N,*this);
    setpointers();
}


void ListSphere::DecreaseLabel(void) noexcept
{
    for (Sphere* mysphere : list)
    {
        mysphere->DecreaseLabel();
    }
}


void ListSphere::setpointers()
{
    vector<double>::iterator newdeb((*Storage)[1].begin());
    vector<double>::iterator newfin((*Storage)[1].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin))
        return;
    for (Sphere* mysphere : list)
    {
        mysphere->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}



/** Default constructor in local storage */
ListSphere::ListSphere(void):
    storage_list<7,Sphere>(),
    physicalmodel(new PhysicalModel),
    ptr_deb(nullptr),
    ptr_fin(nullptr)
{}

ListSphere::ListSphere(PhysicalModel& _physicalmodel, const int _N) :
    storage_list<7,Sphere>(),
    physicalmodel(&_physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr)
{
    Init(_physicalmodel,_N);
}

/** Constructor with external storage */
ListSphere::ListSphere(ListSphere& parent,int _index[]):
    storage_list<7,Sphere>(parent, _index),
    physicalmodel(parent.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr)
{
    setpointers();
}

ListSphere::ListSphere(ListSphere& parent,int* _index[],const int start,const int end):
    storage_list<7,Sphere>(parent, _index, start, end),
    physicalmodel(parent.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr)
{
    setpointers();
}

/** Copy constructor */
ListSphere::ListSphere(const ListSphere& other):
    storage_list<7,Sphere>(other),
    physicalmodel(other.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr)
{
    setpointers();
}

/** Move constructor */
ListSphere::ListSphere (ListSphere&& other) noexcept: /* noexcept needed to enable optimizations in containers */
    storage_list<7,Sphere>(other),
    physicalmodel(other.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr)
{
    setpointers();
}

/** Destructor */
ListSphere::~ListSphere(void) noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;
}

/** Copy assignment operator */
ListSphere& ListSphere::operator= (const ListSphere& other)
{
    ListSphere tmp(other);      // re-use copy-constructor
    *this = std::move(tmp);     // re-use move-assignment
    setpointers();
    return *this;
}

/** Move assignment operator */
ListSphere& ListSphere::operator= (ListSphere&& other) noexcept
{
    if(physicalmodel->toBeDestroyed)
        delete physicalmodel;

    std::swap(static_cast<storage_list<7,Sphere>&>(*this),static_cast<storage_list<7,Sphere>&>(other));
    physicalmodel = other.physicalmodel;
    other.physicalmodel = new PhysicalModel;
    setpointers();
    return *this;
}
