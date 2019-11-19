#ifndef VERLET_H
#define VERLET_H 1

#include <array>
#include <list>
#include <vector>

namespace MCAC{

class Verlet
{
public:
    std::vector< std::vector< std::vector< std::list< size_t > > > > Cell;
private:
    size_t GridDiv{0};
    double L{0};

public:
    void Remove(size_t id, std::array<size_t, 3> Index);
    void Add(size_t id, std::array<size_t, 3> Index);
    void print(std::array<size_t, 3> Index) const;
    void search(size_t id) const;

    void Init(size_t _GridDiv, double _L);
    std::vector<size_t> GetSearchSpace(std::array<double, 3> sourceposition , double width, std::array<double, 3> Vector) const;
    std::vector<size_t> GetSearchSpace(std::array<double, 3> sourceposition , double width) const;

    /* Default constructor */
    Verlet();
};
}  // namespace MCAC

#endif // VERLET_H
