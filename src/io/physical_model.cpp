#include "physical_model.hpp"


#ifdef WITH_HDF5
#include "io/format.hpp"


namespace MCAC {
shared_ptr<XdmfInformation> PhysicalModel::xmf_write() const {
    shared_ptr<XdmfInformation> info = XdmfInformation::New("Physics", "Physical properties of the simulation");
    info->insert(xmf_format_double("Asurfgrowth", Asurfgrowth));
    info->insert(xmf_format_double("dfe", dfe));
    info->insert(xmf_format_double("kfe", kfe));
    info->insert(xmf_format_double("xsurfgrowth", xsurfgrowth));
    info->insert(xmf_format_double("coeffB", coeffB));
    info->insert(xmf_format_double("lambda", lambda));
    info->insert(xmf_format_double("Dpeqmass", Dpeqmass));
    info->insert(xmf_format_double("rpeqmass", rpeqmass));
    info->insert(xmf_format_double("gamma_", gamma_));
    info->insert(xmf_format_double("P [Pa]", P));
    info->insert(xmf_format_double("T [K]", T));
    info->insert(xmf_format_double("Mu", Mu));
    info->insert(xmf_format_double("K", K));
    info->insert(xmf_format_double("Rho [kg/m3]", Rho));
    info->insert(xmf_format_double("Dpm [nm]", Dpm));
    info->insert(xmf_format_double("sigmaDpm [nm]", sigmaDpm));
    info->insert(xmf_format_double("X []", X));
    info->insert(xmf_format_double("FV [ppt]", FV));
    info->insert(xmf_format_double("L", L));
    info->insert(xmf_format_integer("N []", N));
    info->insert(xmf_format_bool("ActiveModulephysique", ActiveModulephysique));
    return info;
}
} // namespace MCAC

#else

auto PhysicalModel::xmfWrite() const
{
    return false;
}

#endif
