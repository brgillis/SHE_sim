/**********************************************************************\
 @file galaxy_observed_properties.cpp
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
#include <stdexcept>

#include <SHE_SIM_gal_params/dependency_functions/misc_dependencies.hpp>
#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/default_values.hpp"

#include "IceBRG_main/math/random/random_functions.hpp"

#include "IceBRG_physics/luminosity.hpp"

namespace SHE_SIM {

flt_t get_apparent_mag_vis( flt_t const & absolute_mag_vis, flt_t const & redshift )
{
	return IceBRG::get_app_mag_from_abs_mag(absolute_mag_vis,redshift);
}

flt_t generate_apparent_size_bulge( flt_t const & apparent_mag_vis, gen_t & rng )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_apparent_size_bulge used.";

	return dv::apparent_size_bulge;
}

flt_t generate_apparent_size_disk( flt_t const & apparent_mag_vis, gen_t & rng )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_apparent_size_disk used.";

	return dv::apparent_size_disk;
}

flt_t generate_shear_angle( flt_t const & xp, flt_t const & yp, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_shear_angle used.";

	return IceBRG::drand( dv::shear_angle_min, dv::shear_angle_max, rng );
}

flt_t generate_shear_magnitude( flt_t const & xp, flt_t const & yp, flt_t const & redshift, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_shear_magnitude used.";

	return IceBRG::contracted_Rayleigh_rand( dv::shear_magnitude_sigma,
			dv::shear_magnitude_max, dv::shear_magnitude_p , rng );
}


} // namespace SHE_SIM
