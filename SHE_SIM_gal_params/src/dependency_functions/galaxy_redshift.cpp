/**********************************************************************\
 @file galaxy_redshift.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2016 brg

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

#include <cmath>

#include <boost/log/trivial.hpp>

#include "SHE_SIM_gal_params/dependency_functions/galaxy_redshift.hpp"
#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "IceBRG_main/math/misc_math.hpp"
#include "IceBRG_main/math/random/random_functions.hpp"
#include "IceBRG_main/units/units.hpp"
#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_physics/astro.h"
#include "IceBRG_physics/astro_caches.h"

namespace SHE_SIM {

flt_t get_cluster_pz( flt_t const & z )
{
	IceBRG::square_angle_type sq_arcmin = 1.*IceBRG::square(IceBRG::unitconv::amintorad*IceBRG::rad);
	IceBRG::custom_unit_type<0,0,0,-2,0> cluster_angular_density = IceBRG::cluster_angular_density_at_z(z);
	return cluster_angular_density * sq_arcmin;
}

flt_t get_total_pz( flt_t const & z )
{
	IceBRG::square_angle_type sq_arcmin = 1.*IceBRG::square(IceBRG::unitconv::amintorad*IceBRG::rad);
	IceBRG::custom_unit_type<0,0,0,-2,0> total_angular_density = IceBRG::galaxy_angular_density_at_z(z);
	return total_angular_density * sq_arcmin;
}

flt_t get_field_pz( flt_t const & z, flt_t const & total_enhancement,
		flt_t const & cluster_enhancement )
{
	flt_t p = total_enhancement*get_total_pz(z) - cluster_enhancement*get_cluster_pz(z);
	if(p < 0) return 0;
	return p;
}

flt_t generate_cluster_z( flt_t const & z_min, flt_t const & z_max, gen_t & rng )
{
	auto cluster_pz = [] (flt_t const & z) {return get_cluster_pz(z);};
	return IceBRG::rand_from_pdf(cluster_pz,40,z_min,z_max,rng);
}

flt_t generate_field_z( flt_t const & total_enhancement,
		flt_t const & cluster_enhancement,
		flt_t const & z_min, flt_t const & z_max, gen_t & rng )
{
	auto get_field_pz_for_zm = [&] (flt_t const & z)
		{
			return get_field_pz(z, total_enhancement, cluster_enhancement);
		};

	// Use a fallback here in case the cluster distribution is bad
	try
	{
		return IceBRG::rand_from_pdf(get_field_pz_for_zm,40,z_min,z_max,rng);
	}
	catch(std::runtime_error & e)
	{
		BOOST_LOG_TRIVIAL(warning) << "Cluster N(z) is uniformly greater than total N(z).";
		return IceBRG::rand_from_pdf(&get_total_pz,40,z_min,z_max,rng);
	}
}

flt_t get_cluster_enhancement(flt_t const & cluster_density, flt_t const & z_min, flt_t const & z_max)
{
	IceBRG::square_angle_type sq_arcmin = 1.*IceBRG::square(IceBRG::unitconv::amintorad*IceBRG::rad);
	flt_t num_clusters_per_sq_arcmin = IceBRG::visible_clusters(sq_arcmin,z_min,z_max);

	return cluster_density/num_clusters_per_sq_arcmin;
}

flt_t get_total_enhancement(flt_t const & total_density, flt_t const & z_min, flt_t const & z_max)
{
	IceBRG::square_angle_type sq_arcmin = 1.*IceBRG::square(IceBRG::unitconv::amintorad*IceBRG::rad);
	flt_t num_galaxies_per_sq_arcmin = IceBRG::visible_galaxies(sq_arcmin,z_min,z_max);

	return total_density/num_galaxies_per_sq_arcmin;
}

flt_t get_ex_num_cluster_galaxies(flt_t const & num_clusters, flt_t const & z_min, flt_t const & z_max)
{
	flt_t mean_richness = IceBRG::mean_cluster_richness(z_min,z_max);

	return num_clusters*mean_richness;
}

} // namespace SHE_SIM
