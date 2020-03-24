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
#ifndef INCLUDE_SPHERES_SPHERE_LIST_HPP
#define INCLUDE_SPHERES_SPHERE_LIST_HPP 1
#include "constants.hpp"
#include "io/threaded_io.hpp"
#include "list_storage/list_storage.hpp"
#include "physical_model/physical_model.hpp"
#include <gsl/gsl>

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
namespace mcac {
class Sphere;

class SphereList : public ListStorage<SpheresFields::SPHERE_NFIELDS, Sphere> {
    friend class Sphere;

    /* Generic */

private:
    std::vector<double>::iterator ptr_deb;
    std::vector<double>::iterator ptr_fin;
    gsl::owner<ThreadedIO *> writer;
    size_t last_saved;
public:
    const PhysicalModel *physicalmodel;
    void init(const PhysicalModel &physical_model, size_t size);
    void decrease_label() noexcept;
    void croissance_surface(double dt) noexcept;
    void print() const;
    void save() {
        save(false);
    }
    void save(bool finish);
    auto get_data() const;
    std::vector<double> format_position() const;
    std::vector<double> format_get_radius() const;
    std::vector<long> format_agg_label() const;

    /* Storage specific */
private:
    void setpointers();
public:
    /** Default constructor in local storage */
    SphereList() noexcept;
    SphereList(const PhysicalModel &physical_model, size_t size) noexcept;
    /** Constructor with external storage */
    SphereList(SphereList *parent, const std::vector<size_t> &index) noexcept;
    /** Copy constructor */
    explicit SphereList(const SphereList &other) noexcept = delete;
    SphereList(const SphereList &other, SphereList *sphere_list) noexcept; // TODO delete
    /** Move constructor */
    SphereList(SphereList &&) noexcept; // TODO delete
    /** Destructor */
    ~SphereList() noexcept;
    /** Copy assignment operator */
    SphereList &operator=(const SphereList &other) noexcept = delete;
    /** Move assignment operator */
    SphereList &operator=(SphereList &&other) noexcept;
};
}  // namespace mcac

#endif //INCLUDE_SPHERES_SPHERE_LIST_HPP
