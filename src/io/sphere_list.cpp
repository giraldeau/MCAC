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
#include "spheres/sphere.hpp"
#include "spheres/sphere_list.hpp"


#ifdef WITH_HDF5
#include "io/format.hpp"
#include "io/writer.hpp"


namespace mcac {
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<double> &formated_field);
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<long> &formated_field);
extern template shared_ptr<XdmfAttribute> attribute(const std::string &name,
                                                    const double &value);
auto SphereList::get_data() const {
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> spheres_data = XdmfUnstructuredGrid::New();
    spheres_data->setName("Spheres");
    spheres_data->setTopology(the_topology());

    // Set time
    spheres_data->setTime(format_time(physicalmodel->time));
    spheres_data->insert(attribute("Time", physicalmodel->time));
    spheres_data->insert(attribute("BoxSize", physicalmodel->box_lenght));

    // Set Positions
    spheres_data->setGeometry(the_positions(format_position()));

    // Set Radius
    spheres_data->insert(scalar("Radius", format_get_radius()));
    spheres_data->insert(scalar("Label", format_agg_label()));
    return spheres_data;
}
void SphereList::save(const bool _finish) {
    //get everything
    auto data = get_data();
    writer->write(physicalmodel->output_dir / "Spheres", data, _finish);
}
DEF_FORMATER_POSITION(SphereList)
DEF_FORMATER(SphereList, agg_label, long)
DEF_FORMATER_FUNC(SphereList, get_radius, double)
}  // namespace mcac


#else

namespace MCAC{

void SphereList::save(const bool finish){}

}  // namespace mcac


#endif
