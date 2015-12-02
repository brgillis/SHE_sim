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
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM {

flt_t generate_cluster_mass( flt_t const & cluster_redshift, gen_t & rng )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_cluster_mass used.";

	return log10Gaus_rand(dv::cluster_mass_l10_mean, dv::cluster_mass_l10_stddev, rng);
}

flt_t generate_cluster_redshift( gen_t & rng )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_cluster_redshift used.";

	return drand(dv::cluster_redshift_min, dv::cluster_redshift_max, rng);
}

flt_t get_cluster_richness( flt_t const & cluster_mass, flt_t const & cluster_redshift )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_cluster_richness used.";

	return dv::cluster_richness;
}

} // namespace SHE_SIM
