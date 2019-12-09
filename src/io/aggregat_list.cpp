#include "aggregats/aggregat_list.hpp"
#include "io/format.hpp"
#include "io/writer.hpp"
#include "io/xmf_includes.hpp"


namespace mcac {
#ifdef WITH_HDF5
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<double> &formated_field);
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<int> &formated_field);
extern template shared_ptr<XdmfAttribute> attribute(const std::string &name,
                                                    const double &value);
auto AggregatList::get_data() const {
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> aggregats_data = XdmfUnstructuredGrid::New();
    aggregats_data->setName("Aggregats");
    aggregats_data->setTopology(the_topology());

    // Set time
    shared_ptr<XdmfTime> time = XdmfTime::New(physicalmodel->time);
    aggregats_data->setTime(time);
    aggregats_data->insert(attribute("Time", physicalmodel->time));
    aggregats_data->insert(attribute("BoxSize", physicalmodel->box_lenght));

    // Set Positions
    aggregats_data->setGeometry(the_positions(format_position()));

    // Set Radius
    aggregats_data->insert(scalar("Rg", format_rg()));
    aggregats_data->insert(scalar("Np", format_n_spheres()));
    aggregats_data->insert(scalar("f_agg", format_f_agg()));
    aggregats_data->insert(scalar("lpm", format_lpm()));
    aggregats_data->insert(scalar("Deltat", format_time_step()));
    aggregats_data->insert(scalar("Rmax", format_rmax()));
    aggregats_data->insert(scalar("compute_volume", format_agregat_volume()));
    aggregats_data->insert(scalar("Surface", format_agregat_surface()));
    aggregats_data->insert(scalar("Label", format_label()));
    return aggregats_data;
}
void AggregatList::save(const bool _finish) {
    auto data = get_data();
    writer->write(physicalmodel->output_dir / "Aggregats", data, _finish);
}
DEF_FORMATER_POSITION(AggregatList)
DEF_FORMATER_PTR(AggregatList, rg, double)
DEF_FORMATER_PTR(AggregatList, f_agg, double)
DEF_FORMATER_PTR(AggregatList, lpm, double)
DEF_FORMATER_PTR(AggregatList, time_step, double)
DEF_FORMATER_PTR(AggregatList, rmax, double)
DEF_FORMATER_PTR(AggregatList, agregat_volume, double)
DEF_FORMATER_PTR(AggregatList, agregat_surface, double)
DEF_FORMATER(AggregatList, n_spheres, int)
DEF_FORMATER(AggregatList, label, int)
#else


void AggregatList::save(const bool finish){}

#endif
}  // namespace mcac

