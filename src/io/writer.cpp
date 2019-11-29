#ifdef WITH_HDF5
#include "io/writer.hpp"
#include <io/threaded_io.hpp>
#include "io/format.hpp"
#include <io/xmf_includes.hpp>
#include <iomanip>
#include <iostream>
#include <gsl/gsl>
#include "XdmfGridCollection.hpp"


namespace fs = std::experimental::filesystem;
namespace MCAC {
// Actual write functions
shared_ptr<XdmfTopology> the_topology() {
    shared_ptr<XdmfTopology> particules = XdmfTopology::New();
    particules->setType(XdmfTopologyType::Polyvertex());
    return particules;
}
shared_ptr<XdmfGeometry> the_positions(const std::vector<double> &formated_positions) {
    const auto N = unsigned(int(formated_positions.size()));
    shared_ptr<XdmfGeometry> positions = XdmfGeometry::New();
    positions->setType(XdmfGeometryType::XYZ());
    positions->insert(0, formated_positions.data(), N, 1, 1);
    return positions;
}
template<class T>
shared_ptr<XdmfAttribute> scalar(const std::string &name, const std::vector<T> &formated_field) {
    const auto N = unsigned(int(formated_field.size()));
    shared_ptr<XdmfAttribute> xdmf_field = XdmfAttribute::New();
    xdmf_field->setName(name);
    xdmf_field->setType(XdmfAttributeType::Scalar());
    xdmf_field->setCenter(XdmfAttributeCenter::Node());
    xdmf_field->insert(0, formated_field.data(), N, 1, 1);
    return xdmf_field;
}
template<class T>
shared_ptr<XdmfAttribute> attribute(const std::string &name, const T &value) {
    shared_ptr<XdmfAttribute> xdmf_field = XdmfAttribute::New();
    xdmf_field->setName(name);
    xdmf_field->setType(XdmfAttributeType::NoAttributeType());
    xdmf_field->insert(0, value);
    return xdmf_field;
}
template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                          const std::vector<double> &formated_field);
template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                          const std::vector<int> &formated_field);
template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                          const std::vector<long> &formated_field);
template shared_ptr<XdmfAttribute> attribute(const std::string &name,
                                             const double &value);
template shared_ptr<XdmfAttribute> attribute(const std::string &name,
                                             const int &value);
void ThreadedIO::wait() {
    if (static_cast<bool>(writer)) {
        writer->join();
        delete writer;
        writer_owner = nullptr;
    }
}
void ThreadedIO::create_file() {
    if (gsl::at(status, static_cast<size_t>(current_thread)) == WriterStatus::WRITING) {
        if (writer_owner == this) {
            wait();
        }
        gsl::at(status, static_cast<size_t>(current_thread)) = WriterStatus::IDLE;
    }
    if (gsl::at(status, static_cast<size_t>(current_thread)) == WriterStatus::APPENDING) {
        std::cout << "And what should happen to the old file ?" << std::endl;
        exit(5);
    }
    // Prepare Xmf file
    gsl::at(xmf_file, static_cast<size_t>(current_thread)) = XdmfDomain::New();
    shared_ptr<XdmfInformation> xmf_info = XdmfInformation::New("Copyright", "Produced by MCAC");
    gsl::at(xmf_file, static_cast<size_t>(current_thread))->insert(xmf_info);

    // Save simulation properties in every xmf file
    gsl::at(xmf_file, static_cast<size_t>(current_thread))->insert(physicalmodel->xmf_write());
    gsl::at(time_collection, static_cast<size_t>(current_thread)) = XdmfGridCollection::New();
    gsl::at(time_collection, static_cast<size_t>(current_thread))->setType(XdmfGridCollectionType::Temporal());
    gsl::at(status, static_cast<size_t>(current_thread)) = WriterStatus::APPENDING;
}
void ThreadedIO::write(const fs::path &prefix, const shared_ptr<XdmfUnstructuredGrid> &data, bool all) {
    if (step % _n_time_per_file == 0) {
        create_file();
    }
    if (gsl::at(status, static_cast<size_t>(current_thread)) != WriterStatus::APPENDING) {
        std::cout << "File is not ready..." << std::endl;
        exit(5);
    }
    gsl::at(time_collection, static_cast<size_t>(current_thread))->insert(data);
    step++;
    if (step % _n_time_per_file == 0 || all) {
        gsl::at(xmf_file, static_cast<size_t>(current_thread))
            ->insert(gsl::at(time_collection, static_cast<size_t>(current_thread)));
        std::string file_name = std::string(prefix) + filename(num_file, _n);
        wait();

        //std::cout << "Writing " << file_name << std::endl;
        if (true) {
            writer =
                new std::thread(write_task, file_name, &gsl::at(xmf_file, static_cast<size_t>(current_thread)));
            writer_owner = this;
            gsl::at(status, static_cast<size_t>(current_thread)) = WriterStatus::WRITING;
        } else { // Sequential
            write_task(file_name, &gsl::at(xmf_file, static_cast<size_t>(current_thread)));
            gsl::at(status, static_cast<size_t>(current_thread)) = WriterStatus::IDLE;
        }
        gsl::at(status, static_cast<size_t>(current_thread)) = WriterStatus::IDLE;
        current_thread = !current_thread;
        num_file++;
        if (all) {
            if (gsl::at(status, static_cast<size_t>(!current_thread)) == WriterStatus::WRITING) {
                if (ThreadedIO::writer_owner == this) {
                    wait();
                }
                gsl::at(status, static_cast<size_t>(!current_thread)) = WriterStatus::IDLE;
            }
            step = 0;
        }
    }
}
void write_task(const std::string &filename, const shared_ptr<XdmfDomain> *data) {
    shared_ptr<XdmfHDF5Writer> hdf5_file = XdmfHDF5Writer::New(filename + ".h5", false);
    shared_ptr<XdmfWriter> xmf_file = XdmfWriter::New(filename + ".xmf", hdf5_file);
    hdf5_file->setUseDeflate(false); // do not use compression (too slow)
    hdf5_file->setDeflateFactor(0);  // 0 to 6, 6 being the most compressed

    xmf_file->setLightDataLimit(0); // everything go to the hdf5

    // Write data
    (*data)->accept(xmf_file);
}
}  // namespace MCAC

#endif
