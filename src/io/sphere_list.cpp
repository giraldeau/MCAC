#include "spheres/Spherelist.hpp"


#ifdef WITH_HDF5
#include "io/format.hpp"
#include "io/writer.hpp"


namespace MCAC {
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<double> &formated_field);
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<long> &formated_field);
extern template shared_ptr<XdmfAttribute> attribute(const std::string &name,
                                                    const double &value);
auto ListSphere::get_data() const {
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> spheres_data = XdmfUnstructuredGrid::New();
    spheres_data->setName("Spheres");
    spheres_data->setTopology(the_topology());

    // Set time
    spheres_data->setTime(format_time(physicalmodel->Time));
    spheres_data->insert(attribute("Time", physicalmodel->Time));
    spheres_data->insert(attribute("BoxSize", physicalmodel->L));

    // Set Positions
    spheres_data->setGeometry(the_positions(Format_Position()));

    // Set Radius
    spheres_data->insert(scalar("Radius", Format_get_radius()));
    spheres_data->insert(scalar("Label", Format_agg_label()));
    return spheres_data;
}
void ListSphere::save(const bool _finish) {
    //get everything
    auto data = get_data();
    Writer->write(physicalmodel->CheminSauve / "Spheres", data, _finish);
}
DEF_FORMATER_POSITION(ListSphere);
DEF_FORMATER(ListSphere, agg_label, long);
DEF_FORMATER_FUNC(ListSphere, get_radius, double);
}  // namespace MCAC


#else

namespace MCAC{

void ListSphere::save(const bool finish){}

}  // namespace MCAC


#endif
