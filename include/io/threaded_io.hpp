#ifndef INCLUDE_IO_THREADED_IO_HPP_
#define INCLUDE_IO_THREADED_IO_HPP_ 1
#include "physical_model.hpp"


#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"
#include <array>
#include <string>
#include <thread>
#include <gsl/gsl>
#include <experimental/filesystem>


namespace MCAC {
class ThreadedIO {
public:
    static void wait();
    void create_file();
    void write(const std::experimental::filesystem::path &prefix,
               const shared_ptr<XdmfUnstructuredGrid> &data,
               bool all);
private:
    const PhysicalModel *physicalmodel;
    // Number of time steps per file
    const size_t _n_time_per_file;
    // estimate data size
    const size_t _n;
    static gsl::owner<std::thread *> writer;
    static ThreadedIO *writer_owner;
    bool current_thread;
    std::array<Writer_status, 2> status;
    std::array<shared_ptr<XdmfDomain>, 2> xmf_file;
    std::array<shared_ptr<XdmfGridCollection>, 2> time_collection;
    // Total number of time steps written (or ready to be saved)
    size_t step;
    // next file number
    int num_file;

    //    int padding;

public:
    ThreadedIO(const PhysicalModel &physicalmodel, size_t size);
    ~ThreadedIO() noexcept;
    /** Copy constructor */
    ThreadedIO(ThreadedIO const &) = delete;
    /** Copy assignment operator */
    void operator=(ThreadedIO const &) = delete;
    /** Move constructor */
    ThreadedIO(ThreadedIO &&) = delete;
    /** Move assignment operator */
    ThreadedIO &operator=(ThreadedIO &&other) = delete;
private:
    friend void write_task(const std::string &filename, shared_ptr<XdmfDomain> *data);
};
}  // namespace MCAC

#else //WITHOUT_HDF5

namespace MCAC{

class ThreadedIO{
public:
    ThreadedIO(PhysicalModel& _physicalmodel, size_t size)  {}
};
}  // namespace MCAC

#endif
#endif //INCLUDE_IO_THREADED_IO_HPP_

