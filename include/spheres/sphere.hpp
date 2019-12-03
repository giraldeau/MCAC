#ifndef INCLUDE_SPHERES_SPHERE_HPP_
#define INCLUDE_SPHERES_SPHERE_HPP_
#include "physical_model.hpp"
#include "storage.hpp"
#include "cst.hpp"
#include <array>
#include <vector>
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
namespace MCAC {
class Aggregate;

class SphereList;

class ListAggregat;

class Sphere : public storage_elem<SpheresFields::SPHERE_NFIELDS, SphereList> {
    friend class SphereList; // TODO remove (pointers)
    friend class ListAggregat; // TODO remove (pointers)

    /* Generic */

private:
    double *x;
    double *y;
    double *z;
    double *r;
    double *rx;
    double *ry;
    double *rz;
    double *volume;
    double *surface;
    void update_vol_and_surf() noexcept;
public:
    long agg_label;
    const PhysicalModel *physicalmodel;
    /* getters */
    double get_volume() const noexcept;
    double get_surface() const noexcept;
    double get_radius() const noexcept;
    std::array<double, 3> get_position() const noexcept;
    std::array<double, 3> get_relative_position() const noexcept;
    /* modifiers */
    void set_label(long value) noexcept;
    void decrease_label() noexcept;
    void set_position(std::array<double, 3> newposition) noexcept;
    void translate(const std::array<double, 3>& trans) noexcept;

    void relative_translate(std::array<double, 3> trans) noexcept;
    void init_val() noexcept;
    void init_val(std::array<double, 3> newposition, double newr) noexcept;
    void croissance_surface(double dt) noexcept;
    /* other */
    void print() const noexcept;

    /* Storage specific */
private:
    void setpointers();
public:
    /** Default constructor in local storage */
    Sphere() noexcept;
    explicit Sphere(const PhysicalModel &) noexcept;
    explicit Sphere(const Aggregate &) noexcept;
    /** Constructor in local storage with initialization */
    Sphere(const PhysicalModel &, const std::array<double, 3> &newposition, double newr) noexcept;
    /** Constructor with external storage */
    Sphere(SphereList &aggregat, size_t id) noexcept;
    /** Destructor */
    ~Sphere() noexcept;
    /** Copy constructor */
    Sphere(const Sphere &) noexcept = delete;
    Sphere(const Sphere &, SphereList &aggregat, size_t id) noexcept;
    /** Move constructor */
    Sphere(Sphere &&) noexcept = delete;
    /** Copy assignment operator */
    Sphere &operator=(const Sphere &other) noexcept = delete;
    /** Move assignment operator */
    Sphere &operator=(Sphere &&other) noexcept = delete;
};
}  // namespace MCAC

#endif //INCLUDE_SPHERES_SPHERE_HPP_
