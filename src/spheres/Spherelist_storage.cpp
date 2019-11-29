
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



#include "spheres/Spherelist.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

using namespace std;
namespace MCAC{

void ListSphere::Init(PhysicalModel& _physicalmodel, const size_t _size)
{
    if(physicalmodel && physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
    }
    if(Writer)
    {
        delete Writer;
    }
    physicalmodel=&_physicalmodel;
    Writer = new ThreadedIO(_physicalmodel, _size);
    storage_list<9,Sphere>::Init(_size,*this);
    setpointers();
}


void ListSphere::DecreaseLabel() noexcept
{
    for (Sphere* mysphere : list)
    {
        mysphere->decrease_label();
    }
}


void ListSphere::setpointers()
{
    auto newdeb((*Storage)[0].begin());
    auto newfin((*Storage)[0].end());
    if ((newdeb == ptr_deb) && (newfin == ptr_fin))
    {
        return;
    }
    for (Sphere* mysphere : list)
    {
        mysphere->setpointers();
    }
    ptr_deb = newdeb;
    ptr_fin = newfin;
}




/** Default constructor in local storage */
ListSphere::ListSphere():
    storage_list<9,Sphere>(),
    physicalmodel(nullptr),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (nullptr),
    lastSaved(0)
{}

ListSphere::ListSphere(PhysicalModel& _physicalmodel, const size_t _size) :
    storage_list<9,Sphere>(),
    physicalmodel(&_physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (new ThreadedIO(_physicalmodel,_size)),
    lastSaved(0)
{
    Init(_physicalmodel,_size);
}

/** Constructor with external storage */
ListSphere::ListSphere(ListSphere& parent,vector<size_t>& _index):
    storage_list<9,Sphere>(parent, _index),
    physicalmodel(parent.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (new ThreadedIO(*parent.physicalmodel,size())),
    lastSaved(0)
{
    setpointers();
}
/** Constructor with external storage */
ListSphere::ListSphere(ListSphere& parent,vector<size_t> _index):
    storage_list<9,Sphere>(parent, move(_index)),
    physicalmodel(parent.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (new ThreadedIO(*parent.physicalmodel,size())),
    lastSaved(0)
{
    setpointers();
}
/** Copy constructor */
ListSphere::ListSphere(const ListSphere& other):
    storage_list<9,Sphere>(other,*this),
    physicalmodel(other.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (new ThreadedIO(*physicalmodel,size())),
    lastSaved(other.lastSaved)
{
    for (Sphere* s: list)
    {
        s->physicalmodel = physicalmodel;
    }
    setpointers();
}

ListSphere::ListSphere(const ListSphere& other, ListSphere& _Storage):
    storage_list<9,Sphere>(other,*this, _Storage),
    physicalmodel(other.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (new ThreadedIO(*physicalmodel,size())),
    lastSaved(other.lastSaved)
{
    for (Sphere* s: list)
    {
        s->physicalmodel = physicalmodel;
    }
    setpointers();
}

/** Move constructor */
ListSphere::ListSphere (ListSphere&& other) noexcept: /* noexcept needed to enable optimizations in containers */
    storage_list<9,Sphere>(move(other)),
    physicalmodel(other.physicalmodel),
    ptr_deb(nullptr),
    ptr_fin(nullptr),
    Writer (new ThreadedIO(*other.physicalmodel,size())),
    lastSaved(other.lastSaved)
{
    delete other.Writer;
    setpointers();
}

/** Destructor */
ListSphere::~ListSphere() noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    if(physicalmodel && physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
    }
    if(Writer)
    {
        delete Writer;
    }
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
    if(physicalmodel && physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
    }
    delete Writer;

    physicalmodel = other.physicalmodel;
    Writer=other.Writer;

    lastSaved=other.lastSaved;

    other.physicalmodel = nullptr;
    other.Writer = nullptr;

    storage_list<9,Sphere>::operator=(move(other));

    setpointers();
    return *this;
}

void ListSphere::print() const
{
    cout << "Printing list of "<< size() << " Sphere" << endl;
    if (external_storage)
    {
        cout << "  With external Storage" << endl;
    }
    else
    {
        cout << "  Without external Storage" << endl;
    }
    for (const Sphere* s : list)
    {
        s->print();
    }
}


 __attribute__((pure)) bool operator==(const ListSphere& A, const ListSphere& B)
{
    if (A.physicalmodel != B.physicalmodel)
    {
        return false;
    }
    if (A.ptr_deb != B.ptr_deb)
    {
        return false;
    }
    if (A.ptr_fin != B.ptr_fin)
    {
        return false;
    }
    if (A.list != B.list)
    {
        return false;
    }
    if (A.Storage != B.Storage)
    {
        return false;
    }
    if (A.external_storage != B.external_storage)
    {
       return false;
    }
    return true;
}

 __attribute__((pure)) bool operator!=(const ListSphere& A, const ListSphere& B)
{
    return !(A==B);
}
}  // namespace MCAC

