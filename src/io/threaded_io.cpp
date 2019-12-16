#ifdef WITH_HDF5
#include "constants.hpp"
#include "io/threaded_io.hpp"
#include <gsl/gsl>
#include <iostream>


namespace mcac {
gsl::owner<std::thread *> ThreadedIO::writer = nullptr;
ThreadedIO *ThreadedIO::writer_owner = nullptr;
ThreadedIO::ThreadedIO(const PhysicalModel &new_physicalmodel, size_t size) noexcept :
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
    if (gsl::at(status, static_cast<long>(current_thread)) == WriterStatus::APPENDING
        || gsl::at(status, static_cast<long>(!current_thread)) == WriterStatus::APPENDING) {
        std::cout << "And what should happen to the data ?" << std::endl;
        exit(ErrorCodes::IO_ERROR);
    }
    if (ThreadedIO::writer_owner == this) {
        wait();
    }
}
}  // namespace mcac

#endif
