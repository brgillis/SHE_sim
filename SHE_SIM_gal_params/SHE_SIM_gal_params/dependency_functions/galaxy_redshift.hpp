/**********************************************************************\
 @file galaxy_redshift.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_GALAXY_REDSHIFT_HPP_
#define SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_GALAXY_REDSHIFT_HPP_

#include <cmath>

#include <boost/log/trivial.hpp>

#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "IceBRG_main/math/random/random_functions.hpp"

namespace SHE_SIM {

inline flt_t get_cluster_pz( flt_t const & z, flt_t const & z_m )
{
	return square(z)*exp(-pow(1.412*z/z_m,1.5));
}

inline flt_t get_total_pz( flt_t const & z, flt_t const & z_m )
{
	return square(z)*exp(-pow(1.412*z/z_m,1.5));
}

inline flt_t get_field_pz( flt_t const & z, flt_t const & total_scale, flt_t const & total_z_m,
		flt_t const & cluster_scale, flt_t const & cluster_z_m )
{
	flt_t p = total_scale*get_total_pz(z,total_z_m) - cluster_scale*get_cluster_pz(z,cluster_z_m);
	if(p < 0) return 0;
	return p;
}

inline flt_t generate_cluster_z( flt_t const & z_m, flt_t const & z_min, flt_t const & z_max, gen_t & rng )
{
	auto get_cluster_pz_for_zm = [&z_m] (flt_t const & z)
		{
			return get_total_pz(z,z_m);
		};
	return IceBRG::rand_from_pdf(get_cluster_pz_for_zm,40,z_min,z_max,rng);
}

inline flt_t generate_field_z( flt_t const & total_scale, flt_t const & total_z_m,
		flt_t const & cluster_scale, flt_t const & cluster_z_m,
		flt_t const & z_min, flt_t const & z_max, gen_t & rng )
{
	auto get_field_pz_for_zm = [&] (flt_t const & z)
		{
			return get_field_pz(z, total_scale, total_z_m, cluster_scale, cluster_z_m);
		};

	// Use a fallback here in case the cluster distribution is bad
	try
	{
		return IceBRG::rand_from_pdf(get_field_pz_for_zm,40,z_min,z_max,rng);
	}
	catch(std::runtime_error & e)
	{
		BOOST_LOG_TRIVIAL(warning) << "Cluster N(z) is uniformly greater than total N(z).";
		auto get_total_pz_for_zm = [&total_z_m] (flt_t const & z)
			{
				return get_total_pz(z,total_z_m);
			};
		return IceBRG::rand_from_pdf(get_total_pz_for_zm,40,z_min,z_max,rng);
	}
}

} // namespace SHE_SIM



#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_GALAXY_REDSHIFT_HPP_
