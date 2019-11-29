
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
#include <utility>

using namespace std;

namespace MCAC{

/* Getters */

 __attribute__((pure)) double Sphere::Volume() const noexcept
{
    return *volume;
}

 __attribute__((pure)) double Sphere::Surface() const noexcept
{
    return *surface;
}

 __attribute__((pure)) double Sphere::Radius() const noexcept
{
    return *r;
}

 __attribute__((pure)) array<double, 3> Sphere::GetPosition() const noexcept
{
    array<double, 3> mypos{{*x,*y,*z}};
    return mypos;
}

/* Setters */

void Sphere::SetLabel(const long value) noexcept
{
    AggLabel = value;
}

void Sphere::DecreaseLabel() noexcept
{
    AggLabel--;
}

void Sphere::Translate(const double transx,const double transy,const double transz) noexcept
{
    SetPosition(*x + transx, *y + transy, *z + transz);
}

void Sphere::Translate(const array<double, 3> trans) noexcept
{
    SetPosition(*x + trans[0], *y + trans[1], *z + trans[2]);
}

void Sphere::RelativeTranslate(const double transx,const double transy,const double transz) noexcept
{
    *rx += transx;
    *ry += transy;
    *rz += transz;
}

/* Alias for different type of arguments*/

void Sphere::SetPosition(const array<double, 3> newposition) noexcept
{
    SetPosition(newposition[0],newposition[1],newposition[2]);
}

/* Basic initializer */

void Sphere::InitVal()
{
    InitVal(0, 0, 0, 0.);
}

void Sphere::InitVal(const double newx,const double newy,const double newz,const double newr)
{
    setpointers();
    SetPosition(newx, newy, newz);
    *r = newr;
    *rx=0.;
    *ry=0.;
    *rz=0.;
    UpdateVolAndSurf();
}

void Sphere::InitVal(const array<double, 3> newposition,const double newr)
{
    InitVal(newposition[0],newposition[1],newposition[2],newr);
}

/* format */

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



/* Storage specific */

void Sphere::setpointers()
{
    x = &(*Storage)[0][indexInStorage];
    y = &(*Storage)[1][indexInStorage];
    z = &(*Storage)[2][indexInStorage];
    r = &(*Storage)[3][indexInStorage];
    volume = &(*Storage)[4][indexInStorage];
    surface = &(*Storage)[5][indexInStorage];
    rx = &(*Storage)[6][indexInStorage];
    ry = &(*Storage)[7][indexInStorage];
    rz = &(*Storage)[8][indexInStorage];
}




/** Default constructor in local storage */
Sphere::Sphere():
    storage_elem<9,ListSphere>(),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(nullptr),
    AggLabel(-1)
{
    InitVal();
}

Sphere::Sphere(PhysicalModel& _physicalmodel):
    storage_elem<9,ListSphere>(),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(&_physicalmodel),
    AggLabel(-1)
{
    InitVal();
}


/** Constructor in local storage with initialization */
Sphere::Sphere(PhysicalModel& _physicalmodel,
    const double newx,
    const double newy,
    const double newz,
    const double newr) : Sphere(_physicalmodel)
{
    InitVal(newx, newy, newz, newr);
}
Sphere::Sphere(PhysicalModel& _physicalmodel,
    const array<double, 3> newposition,
    const double newr) :
    Sphere(_physicalmodel,newposition[0],newposition[1],newposition[2],newr){}

/** Constructor with external storage */
Sphere::Sphere(ListSphere& aggregat,const size_t id):
    storage_elem<9,ListSphere>(aggregat,id),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(aggregat.physicalmodel),
    AggLabel(-1)
{
    InitVal();
    external_storage->setpointers();
}

/** Copy constructor */
Sphere::Sphere(const Sphere& other) :
    storage_elem<9,ListSphere>(other),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(other.physicalmodel),
    AggLabel(other.AggLabel)
{
    setpointers();
}
Sphere::Sphere(const Sphere& other, ListSphere& aggregat, size_t id):
    storage_elem<9,ListSphere>(other,*this, aggregat),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
    volume(nullptr),
    surface(nullptr),
    physicalmodel(other.physicalmodel),
    AggLabel(long(id))
{
    setpointers();
}

/** Move constructor */
Sphere::Sphere (Sphere&& other) noexcept : /* noexcept needed to enable optimizations in containers */
    storage_elem<9,ListSphere>(std::move(other)),
    x(nullptr),
    y(nullptr),
    z(nullptr),
    r(nullptr),
    rx(nullptr),
    ry(nullptr),
    rz(nullptr),
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
    other.rx=nullptr;
    other.ry=nullptr;
    other.rz=nullptr;
    other.volume=nullptr;
    other.surface=nullptr;
    other.AggLabel=-1;
}

/** Destructor */
Sphere::~Sphere() noexcept /* explicitly specified destructors should be annotated noexcept as best-practice */
{
    if(physicalmodel && physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
    }
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
    if(physicalmodel && physicalmodel->toBeDestroyed)
    {
        delete physicalmodel;
    }


    physicalmodel = other.physicalmodel;
    AggLabel = other.AggLabel;
    other.AggLabel=-1;
    other.physicalmodel=nullptr;

    storage_elem<9,ListSphere>::operator=(move(other));

    setpointers();

    return *this;
}

void Sphere::print() const
{
    cout << "Printing Sphere " << (indexInStorage) << endl;
    if (!external_storage)
    {
        cout << "  With external Storage" << endl;
    }
    else
    {
        cout << "  Without external Storage" << endl;
    }
    if (AggLabel<0)
    {
        cout << "  This is a virtual Sphere" << endl;
    }
    else
    {
        cout << "  This Sphere is own by the aggregate " << AggLabel << endl;
    }
    cout << "    Position : " << *x << " "<< *y << " "<< *z << endl;
    cout << "    Radius   : " << *r << endl;
    cout << "    Volume   : " << *volume << endl;
    cout << "    Surface  : " << *surface << endl;
    /*
    double* rx;
    double* ry;
    double* rz;
    */


}

}  // namespace MCAC

