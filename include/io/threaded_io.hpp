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
#ifndef INCLUDE_IO_THREADED_IO_HPP
#define INCLUDE_IO_THREADED_IO_HPP 1
#include "physical_model/physical_model.hpp"


#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"
#include <array>
#include <experimental/filesystem>
#include <gsl/gsl>
#include <string>
#include <thread>


namespace fs = std::experimental::filesystem;
namespace mcac {
class ThreadedIO {
public:
    const fs::path prefix;
    static void wait();
    void create_file();
    void write(const shared_ptr<XdmfUnstructuredGrid> &data);
private:
    const PhysicalModel *physicalmodel;
    // Number of time steps per file
    const size_t _n_time_per_file;
    // estimate data size
    const size_t _n;
    static std::unique_ptr<std::thread> writer;
    static ThreadedIO *writer_owner;
    bool current_thread;
    std::array<WriterStatus, 2> status;
    std::array<shared_ptr<XdmfDomain>, 2> xmf_file;
    std::array<shared_ptr<XdmfGridCollection>, 2> time_collection;
    // Total number of time steps written (or ready to be saved)
    size_t step;
    // next file number
    int num_file;

    //    int padding;

public:
    ThreadedIO(const std::experimental::filesystem::path &prefix,
               const PhysicalModel &new_physicalmodel, size_t size) noexcept;
    ~ThreadedIO() noexcept;
    /** Copy constructor */
    explicit ThreadedIO(ThreadedIO const &) = delete;
    /** Copy assignment operator */
    void operator=(ThreadedIO const &) = delete;
    /** Move constructor */
    explicit ThreadedIO(ThreadedIO &&) = delete;
    /** Move assignment operator */
    ThreadedIO &operator=(ThreadedIO &&other) = delete;
private:
    void _write();
    friend void write_task(const std::string &filename, const shared_ptr<XdmfDomain> *data);
};
}  // namespace mcac

#else //WITHOUT_HDF5

namespace mcac{

class ThreadedIO{
public:
    ThreadedIO(const std::experimental::filesystem::path &prefix,
               const PhysicalModel& _physicalmodel, size_t size)  {}
};
}  // namespace mcac

#endif
#endif //INCLUDE_IO_THREADED_IO_HPP
