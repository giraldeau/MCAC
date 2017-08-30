#ifndef VERLET_H
#define VERLET_H

#include <list>
#include <array>
#include <vector>

using namespace std;

class Verlet
{
private:
    list<int>**** verletlist;
    int GridDiv;
    double L;

public:
    void Remove(const int id,const array<int, 4> Index);
    void Add(const int id,const array<int, 4> Index);

    list<int>* GetCell(const int i,const int j,const int k)const;
    void Init(const int GridDiv, const double L);
    void destroy(void);
    vector<int> GetSearchSpace(array<double, 4> sourceposition , const double mindist, const array<double, 4> Vector);

public:
    // Default constructor
    Verlet(void);

    // Copy constructor
    Verlet(const Verlet&);

    // Move constructor
    Verlet (Verlet&&) noexcept; // noexcept needed to enable optimizations in containers

    // Destructor
    ~Verlet(void) noexcept; // explicitly specified destructors should be annotated noexcept as best-practice

    // Copy assignment operator
    Verlet& operator= (const Verlet& other);

    // Move assignment operator
    Verlet& operator= (Verlet&& other) noexcept;
};

#endif // VERLET_H
