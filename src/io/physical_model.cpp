#include "physical_model/physical_model.hpp"


#ifdef WITH_HDF5
#include "io/format.hpp"


namespace MCAC {
shared_ptr<XdmfInformation> PhysicalModel::xmf_write() const {
    shared_ptr<XdmfInformation> info = XdmfInformation::New("Physics", "Physical properties of the simulation");
    info->insert(xmf_format_double("Asurfgrowth", a_surfgrowth));
    info->insert(xmf_format_double("dfe", dfe));
    info->insert(xmf_format_double("kfe", kfe));
    info->insert(xmf_format_double("xsurfgrowth", x_surfgrowth));
    info->insert(xmf_format_double("coeffB", coeff_b));
    info->insert(xmf_format_double("lambda", gaz_mean_free_path));
    info->insert(xmf_format_double("rpeqmass", mean_massic_radius));
    info->insert(xmf_format_double("gamma_", friction_exponnant));
    info->insert(xmf_format_double("P [Pa]", pressure));
    info->insert(xmf_format_double("T [K]", temperature));
    info->insert(xmf_format_double("Mu", viscosity));
    info->insert(xmf_format_double("Rho [kg/m3]", density));
    info->insert(xmf_format_double("Dpm [nm]", mean_diameter));
    info->insert(xmf_format_double("sigmaDpm [nm]", dispersion_diameter));
    info->insert(xmf_format_double("FV [ppt]", volume_fraction));
    info->insert(xmf_format_double("L", box_lenght));
    info->insert(xmf_format_integer("N []", n_monomeres));
    return info;
}
} // namespace MCAC

#else

auto PhysicalModel::xmfWrite() const
{
    return false;
}

#endif
