/**********************************************************************\
 @file galaxy_intrinsic_properties.cpp
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

#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/default_values.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM {

flt_t generate_bulge_fraction( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, flt_t const & morphology, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_bulge_fraction used.";

	return trunc_log10Gaus_rand( dv::bulge_fraction_l10_mean, dv::bulge_fraction_l10_stddev,
			dv::bulge_fraction_l10_min, dv::bulge_fraction_l10_max, rng );
}

flt_t generate_bulge_fraction( flt_t const & apparent_mag_vis, flt_t const & morphology, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_bulge_fraction (alt) used.";

	return trunc_log10Gaus_rand( dv::bulge_fraction_l10_mean, dv::bulge_fraction_l10_stddev,
			dv::bulge_fraction_l10_min, dv::bulge_fraction_l10_max, rng );
}

flt_t generate_morphology( flt_t const & galaxy_type, flt_t const & redshift, flt_t const & stellar_mass, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_morphology used.";

	return drand( dv::morphology_min, dv::morphology_max, rng );
}

flt_t generate_morphology( flt_t const & apparent_mag_vis, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_morphology (alt) used.";

	return drand( dv::morphology_min, dv::morphology_max, rng );
}

flt_t generate_physical_size_bulge( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_physical_size_bulge used.";

	return dv::physical_size_bulge;
}

flt_t generate_physical_size_disk( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_physical_size_disk used.";

	return dv::physical_size_disk;
}

flt_t generate_rotation( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & morphology, flt_t const & stellar_mass, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_rotation used.";

	return drand( dv::rotation_min, dv::rotation_max, rng );
}

flt_t generate_tilt( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & morphology, flt_t const & stellar_mass, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_tilt used.";

	return drand( dv::tilt_min, dv::tilt_max, rng );
}

flt_t generate_stellar_mass( flt_t const & galaxy_type, flt_t const & redshift, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_stellar_mass used.";

	return dv::stellar_mass;
}

} // namespace SHE_SIM


