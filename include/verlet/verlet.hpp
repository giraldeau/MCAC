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
#ifndef INCLUDE_VERLET_VERLET_HPP
#define INCLUDE_VERLET_VERLET_HPP 1
#include <array>
#include <list>
#include <vector>


namespace mcac {
class Verlet {
public:
    std::vector<std::vector<std::vector<std::list<size_t> > > > grid;
private:
    size_t n_div{0};
    double width{0};
public:
    void remove(size_t id, std::array<size_t, 3> index);
    void add(size_t id, std::array<size_t, 3> index);
    std::vector<size_t> get_neighborhood(const std::array<double, 3> &sourceposition,
                                         const std::array<double, 3> &direction,
                                         const double distance) const;
    std::vector<size_t> get_neighborhood(const std::array<double, 3> &sourceposition,
                                         const double distance) const;
    /* Default constructor */
    Verlet(size_t new_n_div, double new_width);
}; // class verlet

}  // namespace mcac

#endif //INCLUDE_VERLET_VERLET_HPP
