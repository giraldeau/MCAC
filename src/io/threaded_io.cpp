#ifdef WITH_HDF5
#include "io/threaded_io.hpp"
#include "cst.hpp"
#include <gsl/gsl>
#include <iomanip>
#include <iostream>


namespace MCAC {
gsl::owner<std::thread *> ThreadedIO::writer = nullptr;
ThreadedIO *ThreadedIO::writer_owner = nullptr;
ThreadedIO::ThreadedIO(const PhysicalModel &physical_model, size_t size) noexcept :
    physicalmodel(&physical_model),
    _n_time_per_file(physical_model.DeltaSauve),
    _n(size),
    current_thread(false),
    status({WriterStatus::IDLE, WriterStatus::IDLE}),
    xmf_file(),
    time_collection(),
    step(0),
    num_file(0) {
}
ThreadedIO::~ThreadedIO() noexcept {
    if (gsl::at(status, static_cast<size_t>(current_thread)) == WriterStatus::APPENDING
        || gsl::at(status, static_cast<size_t>(!current_thread)) == WriterStatus::APPENDING) {
        std::cout << "And what should happen to the data ?" << std::endl;
        exit(5);
    }
    if (ThreadedIO::writer_owner == this) {
        wait();
    }
}
}  // namespace MCAC

#endif
