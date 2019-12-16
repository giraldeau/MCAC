#include "physical_model/physical_model.hpp"


#ifdef WITH_HDF5
#include "io/format.hpp"


namespace mcac {
extern template shared_ptr<XdmfInformation> xmf_format(const std::string &name,
                                                       const double &number);
extern template shared_ptr<XdmfInformation> xmf_format(const std::string &name,
                                                       const int &number);
shared_ptr<XdmfInformation> PhysicalModel::xmf_write() const {
    shared_ptr<XdmfInformation> info = XdmfInformation::New("Physics", "Physical properties of the simulation");
    info->insert(xmf_format("Asurfgrowth", a_surfgrowth));
    info->insert(xmf_format("dfe", fractal_dimension));
    info->insert(xmf_format("kfe", fractal_prefactor));
    info->insert(xmf_format("xsurfgrowth", x_surfgrowth));
    info->insert(xmf_format("coeffB", coeff_b));
    info->insert(xmf_format("lambda", gaz_mean_free_path));
    info->insert(xmf_format("rpeqmass", mean_massic_radius));
    info->insert(xmf_format("gamma_", friction_exponnant));
    info->insert(xmf_format("P [Pa]", pressure));
    info->insert(xmf_format("T [K]", temperature));
    info->insert(xmf_format("Mu", viscosity));
    info->insert(xmf_format("Rho [kg/m3]", density));
    info->insert(xmf_format("Dpm [nm]", mean_diameter));
    info->insert(xmf_format("sigmaDpm [nm]", dispersion_diameter));
    info->insert(xmf_format("FV [ppt]", volume_fraction));
    info->insert(xmf_format("L", box_lenght));
    info->insert(xmf_format("N []", static_cast<int>(n_monomeres)));
    return info;
}
} // namespace mcac

#else

auto PhysicalModel::xmfWrite() const
{
    return false;
}

#endif
