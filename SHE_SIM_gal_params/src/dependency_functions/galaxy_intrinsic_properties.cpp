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
#include <cassert>
#include <stdexcept>

#include "IceBRG_main/math/random/random_functions.hpp"
#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_physics/abundance_matching.hpp"
#include "IceBRG_physics/luminosity.hpp"

#include "SHE_SIM_gal_params/dependency_functions/galaxy_type.hpp"
#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/default_values.hpp"

namespace SHE_SIM {

using namespace IceBRG;

flt_t generate_rotation( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & sersic_index, flt_t const & stellar_mass, gen_t & rng  )
{
	return IceBRG::drand( dv::rotation_min, dv::rotation_max, rng );
}

flt_t generate_tilt( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & sersic_index, flt_t const & stellar_mass, gen_t & rng  )
{
	return IceBRG::drand( dv::tilt_min, dv::tilt_max, rng );
}

flt_t get_central_abs_mag_vis( flt_t const & cluster_mass, flt_t const & redshift )
{
	flt_t central_abs_mag_B = get_abs_mag_B_from_mass(
			units_cast<mass_type>(cluster_mass*unitconv::Msuntokg*kg), redshift);
	flt_t central_abs_mag_vis = estimate_abs_mag_vis_from_abs_mag_g(central_abs_mag_B);

	return central_abs_mag_vis;
}

constexpr flt_t abs_mag_vis_scatter = 0.5; // FIXME Get realistic value for this

flt_t generate_central_abs_mag_vis( flt_t const & redshift, flt_t const & cluster_mass, gen_t & rng  )
{
	flt_t central_abs_mag_vis = get_central_abs_mag_vis( cluster_mass, redshift );

	flt_t max_visible_abs_mag_i = get_abs_mag_from_app_mag( faint_app_mag_i_max, redshift);
	flt_t max_visible_abs_mag_vis = estimate_abs_mag_vis_from_abs_mag_i(max_visible_abs_mag_i);

	if( central_abs_mag_vis > max_visible_abs_mag_vis ) central_abs_mag_vis = max_visible_abs_mag_vis;

	return central_abs_mag_vis;
}

flt_t generate_field_abs_mag_vis( flt_t const & redshift, gen_t & rng  )
{
	flt_t max_visible_abs_mag_i = get_abs_mag_from_app_mag( faint_app_mag_i_max, redshift);
	flt_t max_visible_abs_mag_g = estimate_abs_mag_g_from_abs_mag_i(max_visible_abs_mag_i);

	flt_t abs_mag_g = rand_from_pdf(&differential_luminosity_function,40,lum_func_min_abs_mag_B,
			max_visible_abs_mag_g,rng);

	flt_t abs_mag_vis = estimate_abs_mag_vis_from_abs_mag_g(abs_mag_g);

	return abs_mag_vis;
}

flt_t generate_satellite_abs_mag_vis( flt_t const & redshift, flt_t const & cluster_mass, gen_t & rng  )
{
	flt_t central_abs_mag_vis = get_central_abs_mag_vis( redshift, cluster_mass );

	flt_t satellite_abs_mag_vis;

	int i=0;
	do
	{
		satellite_abs_mag_vis = generate_field_abs_mag_vis(redshift,rng);
		++i;
	}
	while((satellite_abs_mag_vis<central_abs_mag_vis) and (i<100));

	if(i<100)
		return satellite_abs_mag_vis;
	else
		return central_abs_mag_vis/2.; // Fallback, to make sure central is brightest
}

flt_t generate_abs_mag_vis(flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & cluster_mass, gen_t & rng)
{
	flt_t abs_mag_vis = 0.;

	if(is_field_galaxy(galaxy_type))
	{
		abs_mag_vis = generate_field_abs_mag_vis(redshift,rng);
	}
	else if(is_satellite_galaxy(galaxy_type))
	{
		abs_mag_vis = generate_satellite_abs_mag_vis(redshift,cluster_mass,rng);
	}
	else if(is_central_galaxy(galaxy_type))
	{
		abs_mag_vis = generate_central_abs_mag_vis(redshift,cluster_mass,rng);
	}
	else
	{
		assert(false); // to catch errors if something changes with galaxy type
	}

	return abs_mag_vis;
}

flt_t get_stellar_mass( flt_t const & abs_mag_vis )
{
	flt_t stellar_mass_msun = estimate_stellar_mass_from_abs_mag_i(abs_mag_vis) /
			(unitconv::Msuntokg*kg);

	return stellar_mass_msun;
}

} // namespace SHE_SIM


