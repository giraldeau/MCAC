#ifndef INCLUDE_IO_WRITER_HPP_
#define INCLUDE_IO_WRITER_HPP_ 1
#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"


namespace MCAC {
void write_task(const std::string &filename, const shared_ptr<XdmfDomain> *data);
shared_ptr<XdmfTopology> the_topology();
shared_ptr<XdmfGeometry> the_positions(const std::vector<double> &formated_positions);
template<class T>
shared_ptr<XdmfAttribute> scalar(const std::string &name, const std::vector<T> &formated_field);
template<class T>
shared_ptr<XdmfAttribute> attribute(const std::string &name, const T &value);
}  // namespace MCAC

#endif // WITH_HDF5
#endif //INCLUDE_IO_WRITER_HPP_

