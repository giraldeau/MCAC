#ifndef VERLET_H
#define VERLET_H

#include <list>
#include <array>
#include <vector>

namespace DLCA{

class Verlet : public std::vector< std::vector< std::vector< std::list< int > > > >
{
private:
    int GridDiv{0};
    double L{0};

public:
    void Remove(const int id,const std::array<int, 3> Index);
    void Add(const int id,const std::array<int, 3> Index);
    void print(const std::array<int, 3> Index) const;
    void search(const int id) const;

    void Init(const int GridDiv, const double L);
    std::vector<int> GetSearchSpace(const std::array<double, 3> sourceposition , const double mindist, const std::array<double, 3> Vector) const;
};
}
#endif // VERLET_H
