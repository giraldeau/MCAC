/*
 * MCAC
 * Copyright (C) 2020 CORIA
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "physical_model/physical_model.hpp"
#include "tools/tools.hpp"
#include "verlet/verlet.hpp"
#include <cmath>


namespace mcac {
void Verlet::remove(size_t id, std::array<size_t, 3> index) {
    auto cell = &grid[index[0]][index[1]][index[2]];
    auto it = cell->find(id);
    if (it != cell->end()) {
        cell->erase(it);
    }
}
void Verlet::add(size_t id, std::array<size_t, 3> index) {
    auto cell = &grid[index[0]][index[1]][index[2]];
    cell->insert(id);
}
/* Default constructor */
Verlet::Verlet(size_t new_n_div, double new_width) :
    grid(new_n_div) {
    n_div = new_n_div;
    width = new_width;
    for (auto &x_plane : grid) {
        x_plane.resize(n_div);
        for (auto &y_line : x_plane) {
            y_line.resize(n_div);
        }
    }
}
std::vector<size_t> Verlet::get_search_space(std::array<double, 3> sourceposition, double distance) const {
    return get_search_space(sourceposition, distance, {0, 0, 0});
}
std::vector<size_t> Verlet::get_search_space(std::array<double, 3> sourceposition,
                                             double distance,
                                             std::array<double, 3> direction) const {
    double xp{sourceposition[0] + distance + MAX(direction[0], 0)};
    double xm{sourceposition[0] - distance + MIN(direction[0], 0)};
    double yp{sourceposition[1] + distance + MAX(direction[1], 0)};
    double ym{sourceposition[1] - distance + MIN(direction[1], 0)};
    double zp{sourceposition[2] + distance + MAX(direction[2], 0)};
    double zm{sourceposition[2] - distance + MIN(direction[2], 0)};
    auto bornei_1{static_cast<int>(floor(static_cast<double>(n_div) * xm / width))};
    auto bornei_2{static_cast<int>(floor(static_cast<double>(n_div) * xp / width) + 1)};
    auto bornej_1{static_cast<int>(floor(static_cast<double>(n_div) * ym / width))};
    auto bornej_2{static_cast<int>(floor(static_cast<double>(n_div) * yp / width) + 1)};
    auto bornek_1{static_cast<int>(floor(static_cast<double>(n_div) * zm / width))};
    auto bornek_2{static_cast<int>(floor(static_cast<double>(n_div) * zp / width) + 1)};
    if (bornei_2 - bornei_1 >= static_cast<int>(n_div)) {
        bornei_1 = 0;
        bornei_2 = static_cast<int>(n_div) - 1;
    }
    if (bornej_2 - bornej_1 >= static_cast<int>(n_div)) {
        bornej_1 = 0;
        bornej_2 = static_cast<int>(n_div) - 1;
    }
    if (bornek_2 - bornek_1 >= static_cast<int>(n_div)) {
        bornek_1 = 0;
        bornek_2 = static_cast<int>(n_div) - 1;
    }
    std::list<size_t> tmp_search_space;

    // ///////
    for (int i = bornei_1; i <= bornei_2; i++) {
        for (int j = bornej_1; j <= bornej_2; j++) {
            for (int k = bornek_1; k <= bornek_2; k++) {
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
    std::vector<size_t> search_space{make_move_iterator(tmp_search_space.begin()),
                                     make_move_iterator(tmp_search_space.end())};
    return search_space;
}
}  // namespace mcac
