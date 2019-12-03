#ifndef INCLUDE_SPHERES_SPHERE_LIST_HPP_
#define INCLUDE_SPHERES_SPHERE_LIST_HPP_ 1
#include "io/threaded_io.hpp"
#include "physical_model.hpp"
#include "storagelist.hpp"
#include "cst.hpp"

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
class Sphere;

class SphereList : public storage_list<SpheresFields::SPHERE_NFIELDS, Sphere> {
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
    };
    void save(bool finish);
    auto get_data() const;
    std::vector<double> Format_Position() const;
    std::vector<double> Format_get_radius() const;
    std::vector<long> Format_agg_label() const;

    /* Storage specific */
private:
    void setpointers();
public:
    /** Default constructor in local storage */
    SphereList() noexcept;
    SphereList(const PhysicalModel &physical_model, size_t size) noexcept;
    /** Constructor with external storage */
    SphereList(SphereList &parent, const std::vector<size_t> &index) noexcept;
    /** Copy constructor */
    explicit SphereList(const SphereList &other) noexcept = delete;
    SphereList(const SphereList &other, SphereList &storage) noexcept; // TODO delete
    /** Move constructor */
    explicit SphereList(SphereList &&) noexcept; // TODO delete
    /** Destructor */
    ~SphereList() noexcept;
    /** Copy assignment operator */
    SphereList &operator=(const SphereList &other) noexcept = delete;
    /** Move assignment operator */
    SphereList &operator=(SphereList &&other) noexcept;
};
}  // namespace MCAC

#endif //INCLUDE_SPHERES_SPHERE_LIST_HPP_
