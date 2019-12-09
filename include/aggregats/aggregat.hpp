#ifndef INCLUDE_AGGREGATS_AGGREGAT_HPP_
#define INCLUDE_AGGREGATS_AGGREGAT_HPP_ 1
#include "spheres/sphere_list.hpp"
#include "elem_storage/elem_storage.hpp"
#include "constants.hpp"
#include <list>


namespace MCAC {
class AggregatList;

class Verlet;

class Aggregate :
    public ElemStorage<AggregatesFields::AGGREGAT_NFIELDS, AggregatList> {
    friend class AggregatList;

private:
    double *rg;                     // Gyration Radius
    double *f_agg;                  // Friction coeff
    double *lpm;                    // Apparent Mean Free Path
    double *time_step;              // Time to move along lpm
    double *rmax;                   // Radius of the sphere containing Agg
    double *agregat_volume;         // Etimation of the aggregate's volume
    double *agregat_surface;        // Estimation of the aggregate's surface
    double *x, *y, *z;              // position of the gravity center
    double *rx, *ry, *rz;           // position of the gravity center
    double *time;                   // Proper time of the aggregate
    double *dp;
    double *dg_over_dp;
    size_t n_spheres;               // Number of spheres
    size_t label;                   // Uniq label of the aggregat

    std::vector<std::list<std::pair<size_t, double> > > distances;
    std::vector<double> distances_center;
    std::vector<double> volumes;
    std::vector<double> surfaces;
    /*
    bool padding2,padding3,padding4;
    int padding1;
    */
    void update_verlet_index() noexcept;
    void update_distances() noexcept;
    double internal_sphere_distance(size_t i, size_t j) const noexcept;
    double sphere_distance_center(size_t i) const noexcept;
public:
    const PhysicalModel *physicalmodel;
    Verlet *verlet;
    std::array<size_t, 3> index_verlet;
    SphereList myspheres;
    /* getters */
    double get_rg() const noexcept;
    double get_f_agg() const noexcept;
    double get_lpm() const noexcept;
    double get_time_step() const noexcept;
    double get_rmax() const noexcept;
    double get_agregat_volume() const noexcept;
    double get_agregat_surface() const noexcept;
    std::array<double, 3> get_position() const noexcept;
    std::array<double, 3> get_relative_position() const noexcept;
    std::array<size_t, 3> get_verlet_index() const noexcept;
    double get_time() const noexcept;
    size_t size() const noexcept;
    size_t get_label() const noexcept;
    /* modifiers */
    void decrease_label() noexcept;
    void set_verlet(Verlet &) noexcept;
    void unset_verlet() noexcept;
    void set_time(double newtime) noexcept;
    void time_forward(double deltatemps) noexcept;
    void set_position(const std::array<double, 3> &position) noexcept;
    void translate(std::array<double, 3> vector) noexcept;
    //    void init();
    void init(const PhysicalModel &,
              Verlet &,
              const std::array<double, 3> &position,
              size_t new_label,
              SphereList &,
              double sphere_diameter) noexcept;
    void update() noexcept;
    void compute_volume() noexcept;
    void compute_mass_center() noexcept;
    void compute_max_radius() noexcept;
    void compute_giration_radius() noexcept;
    /* other */
    void merge(Aggregate &) noexcept;
    void print() const noexcept;
    /* Storage specific */
private:
    void setpointers() noexcept;
public:
    /** Default constructor in local storage */
    Aggregate() noexcept;
    //    explicit Aggregate(PhysicalModel &);
    /** Constructor with external storage */
    Aggregate(AggregatList &, size_t label) noexcept;
    /** Destructor */
    ~Aggregate() noexcept;
    /** Copy constructor */
    Aggregate(const Aggregate &, AggregatList &) noexcept;
    Aggregate(const Aggregate &) noexcept = delete;
    /** Move constructor */
    Aggregate(Aggregate &&) noexcept = delete;
    /** Copy assignment operator */
    Aggregate &operator=(const Aggregate &other) noexcept = delete;
    /** Move assignment operator */
    Aggregate &operator=(Aggregate &&other) noexcept = delete;
};
}  // namespace MCAC


#endif //INCLUDE_AGGREGATS_AGGREGAT_HPP_
