#ifdef WITH_HDF5
#include "io/xmf_includes.hpp"
#include "io/format.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>


namespace MCAC {
// Usefull tool
std::string to_string(const double &value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}
// Format for Physical Model
shared_ptr<XdmfInformation> xmf_format_double(const std::string &name, const double &number) {
    shared_ptr<XdmfInformation> info = XdmfInformation::New(name, to_string(number));
    return info;
}
shared_ptr<XdmfInformation> xmf_format_integer(const std::string &name, const long &number) {
    shared_ptr<XdmfInformation> info = XdmfInformation::New(name, std::to_string(number));
    return info;
}
shared_ptr<XdmfInformation> xmf_format_bool(const std::string &name, const bool &active) {
    shared_ptr<XdmfInformation> info;
    if (active) {
        info = XdmfInformation::New(name, "True");
    } else {
        info = XdmfInformation::New(name, "False");
    }
    return info;
}
shared_ptr<XdmfTime> format_time(const double &value) {
    shared_ptr<XdmfTime> thetime = XdmfTime::New(value);
    return thetime;
}
// Format filename
std::string filename(int step, size_t N) {
    int witdh = int(ceil(log10(double(N)))) + 4;
    std::ostringstream filename_stream;
    filename_stream << "_" << std::setfill('0') << std::setw(witdh) << step;
    return filename_stream.str();
}
}  // namespace MCAC

#endif
