#ifndef VERLET_H
#define VERLET_H

#include <list>
#include <array>
#include <vector>

using namespace std;

class Verlet : public vector< vector< vector< list< int > > > >
{
private:
    int GridDiv{0};
    double L{0};

public:
    void Remove(const int id,const array<int, 4> Index);
    void Add(const int id,const array<int, 4> Index);

    void Init(const int GridDiv, const double L);
    vector<int> GetSearchSpace(const array<double, 4> sourceposition , const double mindist, const array<double, 4> Vector) const;
};

#endif // VERLET_H
