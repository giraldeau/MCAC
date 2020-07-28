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
#ifdef WITH_HDF5
#include "constants.hpp"
#include "io/threaded_io.hpp"
#include <gsl/gsl>
#include <iostream>


namespace mcac {
std::unique_ptr<std::thread> ThreadedIO::writer{};
ThreadedIO *ThreadedIO::writer_owner = nullptr;
ThreadedIO::ThreadedIO(const std::experimental::filesystem::path &new_prefix,
                       const PhysicalModel &new_physicalmodel, size_t size) noexcept:
    prefix(new_prefix),
    physicalmodel(&new_physicalmodel),
    _n_time_per_file(new_physicalmodel.n_time_per_file),
    _n(size),
    current_thread(false),
    status({WriterStatus::IDLE, WriterStatus::IDLE}),
    xmf_file(),
    time_collection(),
    step(0),
    num_file(0) {
}
ThreadedIO::~ThreadedIO() noexcept {
    if (gsl::at(status, static_cast<long>(current_thread)) == WriterStatus::APPENDING) {
        _write();
    }
    if (ThreadedIO::writer_owner == this) {
        wait();
    }
}
}  // namespace mcac

#endif
