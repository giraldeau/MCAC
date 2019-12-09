#include "verlet/verlet.hpp"
#include "physical_model/physical_model.hpp"
#include "tools/tools.hpp"
#include <cmath>


using namespace std;
namespace MCAC {
void Verlet::remove(size_t id, array<size_t, 3> index) {
    grid[index[0]][index[1]][index[2]].remove(id);
}
void Verlet::add(size_t id, array<size_t, 3> index) {
    grid[index[0]][index[1]][index[2]].push_front(id);
}
/* Default constructor */
Verlet::Verlet(size_t new_n_div, double new_width) {
    grid.clear();
    n_div = new_n_div;
    width = new_width;
    grid.resize(n_div);
    for (vector<vector<list<size_t> > > &x_plane : grid) {
        x_plane.resize(n_div);
        for (vector<list<size_t> > &y_line : x_plane) {
            y_line.resize(n_div);
        }
    }
}
vector<size_t> Verlet::get_search_space(array<double, 3> sourceposition, double distance) const {
    return get_search_space(sourceposition, distance, {0,0,0});
}
vector<size_t> Verlet::get_search_space(array<double, 3> sourceposition,
                                        double distance,
                                        array<double, 3> direction) const {
    double xp{sourceposition[0] + distance + MAX(direction[0], 0)};
    double xm{sourceposition[0] - distance + MIN(direction[0], 0)};
    double yp{sourceposition[1] + distance + MAX(direction[1], 0)};
    double ym{sourceposition[1] - distance + MIN(direction[1], 0)};
    double zp{sourceposition[2] + distance + MAX(direction[2], 0)};
    double zm{sourceposition[2] - distance + MIN(direction[2], 0)};
    auto bornei1{static_cast<int>(floor(xm * static_cast<double>(n_div) / width))};
    auto bornei2{static_cast<int>(floor(xp * static_cast<double>(n_div) / width) + 1)};
    auto bornej1{static_cast<int>(floor(ym * static_cast<double>(n_div) / width))};
    auto bornej2{static_cast<int>(floor(yp * static_cast<double>(n_div) / width) + 1)};
    auto bornek1{static_cast<int>(floor(zm * static_cast<double>(n_div) / width))};
    auto bornek2{static_cast<int>(floor(zp * static_cast<double>(n_div) / width) + 1)};
    if (bornei2 - bornei1 >= static_cast<int>(n_div)) {
        bornei1 = 0;
        bornei2 = static_cast<int>(n_div) - 1;
    }
    if (bornej2 - bornej1 >= static_cast<int>(n_div)) {
        bornej1 = 0;
        bornej2 = static_cast<int>(n_div) - 1;
    }
    if (bornek2 - bornek1 >= static_cast<int>(n_div)) {
        bornek1 = 0;
        bornek2 = static_cast<int>(n_div) - 1;
    }
    list<size_t> tmp_search_space;

    // ///////
    for (int i = bornei1; i <= bornei2; i++) {
        for (int j = bornej1; j <= bornej2; j++) {
            for (int k = bornek1; k <= bornek2; k++) {
                // periodic
                auto ii = static_cast<size_t>(periodic_position(i, static_cast<int>(n_div)));
                auto jj = static_cast<size_t>(periodic_position(j, static_cast<int>(n_div)));
                auto kk = static_cast<size_t>(periodic_position(k, static_cast<int>(n_div)));
                tmp_search_space.insert(tmp_search_space.end(),
                                        grid[ii][jj][kk].begin(),
                                        grid[ii][jj][kk].end());
            }
        }
    }
    vector<size_t> search_space{make_move_iterator(tmp_search_space.begin()),
                                make_move_iterator(tmp_search_space.end())};
    return search_space;
}
}  // namespace MCAC
