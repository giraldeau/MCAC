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
    info->insert(xmf_format("flux_surfgrowth", flux_surfgrowth));
    info->insert(xmf_format("u_sg", u_sg));
    info->insert(xmf_format("dfe", fractal_dimension));
    info->insert(xmf_format("kfe", fractal_prefactor));
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

namespace mcac {
XMF_OUTPUT PhysicalModel::xmf_write() const
{
    return ;
}
} // namespace mcac

#endif
