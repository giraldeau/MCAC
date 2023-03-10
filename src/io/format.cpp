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
#include "io/xmf_includes.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>


namespace mcac {
// Usefull tool
template<class T>
std::string to_string(const T &value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}
// Format for Physical Model
template<class T>
shared_ptr<XdmfInformation> xmf_format(const std::string &name, const T &number) {
    shared_ptr<XdmfInformation> info = XdmfInformation::New(name, to_string(number));
    return info;
}
shared_ptr<XdmfInformation> xmf_format(const std::string &name, const bool &active) {
    shared_ptr<XdmfInformation> info;
    if (active) {
        info = XdmfInformation::New(name, "True");
    } else {
        info = XdmfInformation::New(name, "False");
    }
    return info;
}
template shared_ptr<XdmfInformation> xmf_format(const std::string &name,
                                                const double &number);
template shared_ptr<XdmfInformation> xmf_format(const std::string &name,
                                                const int &number);
template shared_ptr<XdmfInformation> xmf_format(const std::string &name,
                                                const size_t &number);
shared_ptr<XdmfTime> format_time(const double &value) {
    shared_ptr<XdmfTime> thetime = XdmfTime::New(value);
    return thetime;
}
// Format filename
std::string filename(int step, size_t n) {
    int width = int(std::ceil(std::log10(static_cast<float>(n)))) + 4;
    std::ostringstream filename_stream;
    filename_stream << "_" << std::setfill('0') << std::setw(width) << step;
    return filename_stream.str();
}
}  // namespace mcac

#endif
