/**********************************************************************\
 @file cluster_properties.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015 brg

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

\**********************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/log/trivial.hpp>

#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/default_values.hpp"

#include "IceBRG_main/math/random/random_functions.hpp"
#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_main/units/units.hpp"

#include "IceBRG_physics/cluster_visibility.hpp"
#include "IceBRG_physics/mass_function.hpp"

namespace SHE_SIM {

using namespace IceBRG;

flt_t generate_cluster_mass( flt_t const & cluster_redshift, gen_t & rng )
{
	flt_t l10_min_mass = std::log10(min_cluster_mass(cluster_redshift)/(unitconv::Msuntokg*kg));

	auto l10_mass_function_at_z = [&] (flt_t const & l10_m)
	{
		return value_of(log10_mass_function(l10_m,cluster_redshift));
	};

	flt_t l10_mass = rand_from_pdf(l10_mass_function_at_z,40,l10_min_mass,mass_func_l10_max,rng);

	return std::pow(10.,l10_mass)*unitconv::kgtoMsun;
}

flt_t get_cluster_richness( flt_t const & cluster_mass, flt_t const & cluster_redshift )
{
	flt_t richness = cluster_richness(cluster_mass*unitconv::Msuntokg*kg,
			cluster_redshift);

	return richness;
}

} // namespace SHE_SIM
