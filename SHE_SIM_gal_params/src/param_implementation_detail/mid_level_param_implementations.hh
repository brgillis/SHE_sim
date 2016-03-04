/**********************************************************************\
 @file mid_level_param_implementations.hh
 ------------------

 Implementations for Cluster and Field-level parameters, separated out
 for readability.

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

#ifndef SRC_PARAM_IMPLEMENTATION_DETAIL_MID_LEVEL_PARAM_IMPLEMENTATIONS_HH_
#define SRC_PARAM_IMPLEMENTATION_DETAIL_MID_LEVEL_PARAM_IMPLEMENTATIONS_HH_

// Cluster level

IMPLEMENT_PARAM(cluster_mass, dv::cluster_level, Calculated
	,
		_cached_value = generate_cluster_mass(REQUEST(cluster_redshift), get_rng());
	,
		_cached_value = generate_cluster_mass(REQUEST(cluster_redshift), get_rng());
	);
IMPLEMENT_PARAM(cluster_redshift, dv::cluster_level, IndClusterRedshift(
		dv::cluster_redshift_enhancement,
		dv::cluster_redshift_min,
		dv::cluster_redshift_max)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(cluster_xp, dv::cluster_level, IndUniform(dv::cluster_xp_min,
														dv::cluster_xp_max)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(cluster_yp, dv::cluster_level, IndUniform(dv::cluster_yp_min,
														dv::cluster_yp_max)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);

IMPLEMENT_PARAM(cluster_num_satellites, dv::cluster_level, Calculated
	,
		_cached_value = generate_count( get_cluster_richness(REQUEST(cluster_mass),
				REQUEST(cluster_redshift)) - 1, get_rng());
	,
		_cached_value = generate_count( get_cluster_richness(REQUEST(cluster_mass),
				REQUEST(cluster_redshift)) - 1, get_rng());
	);

// Field level

IMPLEMENT_PARAM(num_field_galaxies, dv::field_level, Calculated
	,
		const IndClusterRedshift * p_redshift_pp = dynamic_cast<const IndClusterRedshift *>(
				_request_param(cluster_redshift_name)->get_p_params());
		flt_t z_min;
		flt_t z_max;
		if(p_redshift_pp==nullptr)
		{
			z_min = dv::cluster_redshift_min;
			z_max = dv::cluster_redshift_max;
		}
		else
		{
			z_min = p_redshift_pp->get_z_min();
			z_max = p_redshift_pp->get_z_max();
		}
		 _cached_value = generate_count( (REQUEST(image_area) * REQUEST(galaxy_density)) -
				 get_ex_num_cluster_galaxies(REQUEST(image_area) * REQUEST(cluster_density),
						 z_min, z_max), get_rng());
	,
		const IndClusterRedshift * p_redshift_pp = dynamic_cast<const IndClusterRedshift *>(
				_request_param(cluster_redshift_name)->get_p_params());
		flt_t z_min;
		flt_t z_max;
		if(p_redshift_pp==nullptr)
		{
			z_min = dv::cluster_redshift_min;
			z_max = dv::cluster_redshift_max;
		}
		else
		{
			z_min = p_redshift_pp->get_z_min();
			z_max = p_redshift_pp->get_z_max();
		}
		 _cached_value = generate_count( (REQUEST(image_area) * REQUEST(galaxy_density)) -
				 get_ex_num_cluster_galaxies(REQUEST(image_area) * REQUEST(cluster_density),
						 z_min, z_max), get_rng());
	);



#endif // SRC_PARAM_IMPLEMENTATION_DETAIL_MID_LEVEL_PARAM_IMPLEMENTATIONS_HH_
