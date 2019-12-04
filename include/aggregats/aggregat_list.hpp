#ifndef AGGREGATLIST_H
#define AGGREGATLIST_H 1
#include "spheres/sphere_list.hpp"
#include "aggregat.hpp"
#include "physical_model/physical_model.hpp"
#include "list_storage/list_storage.hpp"
#include "verlet/verlet.hpp"


namespace MCAC {
class AggregatList : public ListStorage<AggregatesFields::AGGREGAT_NFIELDS, Aggregate> {
    friend class Aggregate;

private:
    PhysicalModel *physicalmodel;
    double maxradius;
    double avg_npp;
    double max_time_step;
    std::vector<size_t> index_sorted_time_steps;
    std::vector<double> cumulative_time_steps;
    std::vector<double>::iterator ptr_deb;
    std::vector<double>::iterator ptr_fin;
    gsl::owner<ThreadedIO *> writer;
    size_t last_saved;
public:
    SphereList spheres;
    Verlet verlet;
    /* getters */
    double get_avg_npp() const;
    double get_max_time_step() const;
    double get_time_step(double max) const;
    size_t pick_random() const;
    size_t pick_last() const;
    /* modifiers */
    void refresh();
    void sort_time_steps(double factor);
    void duplication();
    size_t merge(size_t first, size_t second);
    /* other */
    std::tuple<bool, double, double, double> get_instantaneous_fractal_law() const;
    std::pair<size_t, double> distance_to_next_contact(size_t source,
                                                       std::array<double, 3> direction) const;
    std::vector<size_t> get_search_space(size_t source, std::array<double, 3> direction) const;
    std::vector<std::pair<size_t, double> > sort_search_space(size_t moving_aggregate,
                                                              std::array<double, 3> direction,
                                                              const std::vector<size_t> &search_space) const;
    bool test_free_space(std::array<double, 3> pos, double diameter) const;
    /* I/O */
    void save() {
        save(false);
    };
    void save(bool);
    auto get_data() const;
    std::vector<double> Format_Position() const;
    std::vector<double> Format_rg() const;
    std::vector<int> Format_n_spheres() const;
    std::vector<double> Format_f_agg() const;
    std::vector<double> Format_lpm() const;
    std::vector<double> Format_time_step() const;
    std::vector<double> Format_rmax() const;
    std::vector<double> Format_agregat_volume() const;
    std::vector<double> Format_agregat_surface() const;
    std::vector<int> Format_label() const;

    /* Storage specific */
private:
    void setpointers();
public:
    /** Default constructor in local storage */
    explicit AggregatList(PhysicalModel *) noexcept;
    /** Destructor */
    ~AggregatList() noexcept;
    /** Copy constructor */
    AggregatList(const AggregatList &other) noexcept = delete;
    /** Move constructor */
    AggregatList(AggregatList &&) noexcept = delete;
    /** Copy assignment operator */
    AggregatList &operator=(const AggregatList &other) noexcept = delete;
    /** Move assignment operator */
    AggregatList &operator=(AggregatList &&other) noexcept = delete;
};
}// namespace MCAC


#endif // AGGREGATLIST_H
