#include "aggregats/aggregat_list.hpp"
#include "io/xmf_includes.hpp"
#include "io/format.hpp"
#include "io/writer.hpp"


namespace MCAC {
#ifdef WITH_HDF5
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<double> &formated_field);
extern template shared_ptr<XdmfAttribute> scalar(const std::string &name,
                                                 const std::vector<int> &formated_field);
extern template shared_ptr<XdmfAttribute> attribute(const std::string &name,
                                                    const double &value);
auto ListAggregat::get_data() const {
    // Set geometry
    shared_ptr<XdmfUnstructuredGrid> aggregats_data = XdmfUnstructuredGrid::New();
    aggregats_data->setName("Aggregats");
    aggregats_data->setTopology(the_topology());

    // Set time
    shared_ptr<XdmfTime> time = XdmfTime::New(physicalmodel->Time);
    aggregats_data->setTime(time);
    aggregats_data->insert(attribute("Time", physicalmodel->Time));
    aggregats_data->insert(attribute("BoxSize", physicalmodel->L));

    // Set Positions
    aggregats_data->setGeometry(the_positions(Format_Position()));

    // Set Radius
    aggregats_data->insert(scalar("Rg", Format_rg()));
    aggregats_data->insert(scalar("Np", Format_n_spheres()));
    aggregats_data->insert(scalar("f_agg", Format_f_agg()));
    aggregats_data->insert(scalar("lpm", Format_lpm()));
    aggregats_data->insert(scalar("Deltat", Format_time_step()));
    aggregats_data->insert(scalar("Rmax", Format_rmax()));
    aggregats_data->insert(scalar("compute_volume", Format_agregat_volume()));
    aggregats_data->insert(scalar("Surface", Format_agregat_surface()));
    aggregats_data->insert(scalar("Label", Format_label()));
    return aggregats_data;
}
void ListAggregat::save(const bool _finish) {
    auto data = get_data();
    writer->write(physicalmodel->CheminSauve / "Aggregats", data, _finish);
}
DEF_FORMATER_POSITION(ListAggregat);
DEF_FORMATER_PTR(ListAggregat, rg, double);
DEF_FORMATER_PTR(ListAggregat, f_agg, double);
DEF_FORMATER_PTR(ListAggregat, lpm, double);
DEF_FORMATER_PTR(ListAggregat, time_step, double);
DEF_FORMATER_PTR(ListAggregat, rmax, double);
DEF_FORMATER_PTR(ListAggregat, agregat_volume, double);
DEF_FORMATER_PTR(ListAggregat, agregat_surface, double);
DEF_FORMATER(ListAggregat, n_spheres, int);
DEF_FORMATER(ListAggregat, label, int);
#else


void ListAggregat::save(const bool finish){}

#endif
}  // namespace MCAC

