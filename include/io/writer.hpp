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
#ifndef INCLUDE_IO_WRITER_HPP
#define INCLUDE_IO_WRITER_HPP 1
#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"


namespace mcac {
void write_task(const std::string &filename, const shared_ptr<XdmfDomain> *data);
shared_ptr<XdmfTopology> the_topology();
shared_ptr<XdmfGeometry> the_positions(const std::vector<double> &formated_positions);
template<class T>
shared_ptr<XdmfAttribute> scalar(const std::string &name, const std::vector<T> &formated_field);
template<class T>
shared_ptr<XdmfAttribute> attribute(const std::string &name, const T &value);
}  // namespace mcac

#endif // WITH_HDF5
#endif //INCLUDE_IO_WRITER_HPP

