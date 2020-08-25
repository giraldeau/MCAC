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
#include "io/format.hpp"
#include "io/threaded_io.hpp"
#include "io/writer.hpp"
#include "io/xmf_includes.hpp"
#include "exceptions.hpp"
#include <XdmfGridCollection.hpp>
#include <gsl/gsl>
#include <iostream>


namespace fs = std::experimental::filesystem;
namespace mcac {
// Actual write functions
shared_ptr<XdmfTopology> the_topology() {
    shared_ptr<XdmfTopology> particules = XdmfTopology::New();
    particules->setType(XdmfTopologyType::Polyvertex());
    return particules;
}
shared_ptr<XdmfGeometry> the_positions(const std::vector<double> &formated_positions) {
    const auto _n = unsigned(int(formated_positions.size()));
    shared_ptr<XdmfGeometry> positions = XdmfGeometry::New();
    positions->setType(XdmfGeometryType::XYZ());
    positions->insert(0, formated_positions.data(), _n, 1, 1);
    return positions;
}
template<class T>
shared_ptr<XdmfAttribute> scalar(const std::string &name, const std::vector<T> &formated_field) {
    const auto _n = unsigned(int(formated_field.size()));
    shared_ptr<XdmfAttribute> xdmf_field = XdmfAttribute::New();
    xdmf_field->setName(name);
    xdmf_field->setType(XdmfAttributeType::Scalar());
    xdmf_field->setCenter(XdmfAttributeCenter::Node());
    xdmf_field->insert(0, formated_field.data(), _n, 1, 1);
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
void ThreadedIO::create_file() {
    if (gsl::at(status, static_cast<long>(current_thread)) == WriterStatus::WRITING) {
        if (writer_owner == this) {
            wait();
        }
        gsl::at(status, static_cast<long>(current_thread)) = WriterStatus::IDLE;
    }
    if (gsl::at(status, static_cast<long>(current_thread)) == WriterStatus::APPENDING) {
        throw IOError("File is already created");
    }
    // Prepare Xmf file
    gsl::at(xmf_file, static_cast<long>(current_thread)) = XdmfDomain::New();
    shared_ptr<XdmfInformation> xmf_info = XdmfInformation::New("Copyright", "Produced by MCAC");
    gsl::at(xmf_file, static_cast<long>(current_thread))->insert(xmf_info);

    // Save simulation properties in every xmf file
    gsl::at(xmf_file, static_cast<long>(current_thread))->insert(physicalmodel->xmf_write());
    gsl::at(time_collection, static_cast<long>(current_thread)) = XdmfGridCollection::New();
    gsl::at(time_collection, static_cast<long>(current_thread))->setType(XdmfGridCollectionType::Temporal());
    gsl::at(status, static_cast<long>(current_thread)) = WriterStatus::APPENDING;
}
void ThreadedIO::write(const shared_ptr<XdmfUnstructuredGrid> &data) {
    if (step % _n_time_per_file == 0) {
        create_file();
    }
    if (gsl::at(status, static_cast<long>(current_thread)) != WriterStatus::APPENDING) {
        throw IOError("File is not ready...");
    }
    gsl::at(time_collection, static_cast<long>(current_thread))->insert(data);
    step++;
    if (step % _n_time_per_file == 0) {
        _write();
        current_thread = !current_thread;
        num_file++;
    }
}
void ThreadedIO::_write() {
    gsl::at(xmf_file, static_cast<long>(current_thread))
        ->insert(gsl::at(time_collection, static_cast<long>(current_thread)));
    std::string file_name = std::string(prefix) + filename(num_file, _n);

    //std::cout << "Writing " << file_name << std::endl;
    if (true) {
        wait();
        writer = std::make_unique<std::thread>(write_task,
                                               file_name, &gsl::at(xmf_file, static_cast<long>(current_thread)));
        writer_owner = this;
        gsl::at(status, static_cast<long>(current_thread)) = WriterStatus::WRITING;
    } else { // Sequential
        write_task(file_name, &gsl::at(xmf_file, static_cast<long>(current_thread)));
        gsl::at(status, static_cast<long>(current_thread)) = WriterStatus::IDLE;
    }
}
void write_task(const std::string &filename, const shared_ptr<XdmfDomain> *data) {
    shared_ptr<XdmfHDF5Writer> hdf_5_file = XdmfHDF5Writer::New(filename + ".h5", false);
    shared_ptr<XdmfWriter> xmf_file = XdmfWriter::New(filename + ".xmf", hdf_5_file);
    hdf_5_file->setUseDeflate(false); // do not use compression (too slow)
    hdf_5_file->setDeflateFactor(0);  // 0 to 6, 6 being the most compressed

    xmf_file->setLightDataLimit(0); // everything go to the hdf5

    // Write data
    (*data)->accept(xmf_file);
}
void ThreadedIO::wait() {
    if (static_cast<bool>(writer)) {
        writer->join();
        writer.reset();
    }
    writer_owner = nullptr;
}
}  // namespace mcac

#endif
