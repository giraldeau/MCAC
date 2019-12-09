#ifndef INCLUDE_VERLET_VERLET_HPP_
#define INCLUDE_VERLET_VERLET_HPP_ 1
#include <array>
#include <list>
#include <vector>


namespace MCAC {
class Verlet {
public:
    std::vector<std::vector<std::vector<std::list<size_t> > > > grid;
private:
    size_t n_div{0};
    double width{0};
public:
    void remove(size_t id, std::array<size_t, 3> index);
    void add(size_t id, std::array<size_t, 3> index);
    std::vector<size_t> get_search_space(std::array<double, 3> sourceposition,
                                         double distance,
                                         std::array<double, 3> direction) const;
    std::vector<size_t> get_search_space(std::array<double, 3> sourceposition, double distance) const;
    /* Default constructor */
    Verlet(size_t new_n_div, double new_width);
}; // class verlet

}  // namespace MCAC

#endif //INCLUDE_VERLET_VERLET_HPP_
