/**********************************************************************\
 @file dependent_params.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <SHE_SIM_gal_params/param_params/DepFieldRedshift.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/dependency_functions/galaxy_type.hpp"
#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "IceBRG_main/math/random/random_functions.hpp"

namespace SHE_SIM
{

// Define a macro for each param

#define DEPENDENT_PARAM( param_name, dependent_generation ) \
class param_name##_obj : public ParamGenerator \
{ \
private: \
\
	virtual void _generate() override \
	{ \
		if(_p_params->get_mode()==ParamParam::DEPENDENT) \
		{ \
			dependent_generation; \
		} \
		else if(_p_params->get_mode()==ParamParam::INDEPENDENT) \
		{ \
			_cached_value = _p_params->get_independently(_rng); \
		} \
		else \
		{ \
			throw bad_mode_error(_p_params->get_mode_name()); \
		} \
	} \
\
public: \
	param_name##_obj( owner_t & owner) \
	: ParamGenerator(owner) \
	{ \
		/* See if we can get generation level and params from the parent */ \
		auto p_parent_version = _p_parent_version(); \
		if(p_parent_version) \
		{ \
			_p_generation_level = p_parent_version->get_p_generation_level(); \
			_p_params = p_parent_version->get_p_params(); \
		} \
		else \
		{ \
			_p_params = default_param_params_map.at(name()).get(); \
			_p_generation_level = default_generation_levels_map.at(name()).get(); \
		} \
	} \
\
	virtual ~param_name##_obj() \
	{ \
	} \
\
	virtual name_t name() const override \
	{ \
		return param_name##_name; \
	} \
\
	virtual ParamGenerator * clone() const override \
	{ \
		return new param_name##_obj(*this); \
	} \
};

// Define a macro to request a parameter
#define REQUEST(param) _request_param_value(param##_name)

// Define each param

// Survey level

// Image level

DEPENDENT_PARAM(background_noise,
		_cached_value = get_background_noise(REQUEST(subtracted_background), REQUEST(unsubtracted_background),
				REQUEST(read_noise), REQUEST(gain), REQUEST(pixel_scale) ));

DEPENDENT_PARAM(image_area,
		 _cached_value = get_image_area(REQUEST(image_size_xp),REQUEST(image_size_yp),
				 REQUEST(pixel_scale)));

DEPENDENT_PARAM(mag_i_zp,  _cached_value = get_zp(REQUEST(mag_i_inst_zp),
		REQUEST(exp_time)));

DEPENDENT_PARAM(mag_vis_zp,  _cached_value = get_zp(REQUEST(mag_vis_inst_zp),
		REQUEST(exp_time)));

DEPENDENT_PARAM(num_background_galaxies,
		 _cached_value = generate_count( REQUEST(image_area) *
				 REQUEST(background_galaxy_density) , _rng) );

DEPENDENT_PARAM(num_clusters,
		 _cached_value = generate_count( REQUEST(image_area) *
				 REQUEST(cluster_density) , _rng) );

DEPENDENT_PARAM(num_field_galaxies,
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
						 z_min, z_max), _rng) );

DEPENDENT_PARAM(num_stars,
		 _cached_value = generate_count( REQUEST(image_area) *
				 REQUEST(star_density) , _rng) );

// Cluster level

DEPENDENT_PARAM(cluster_mass,
		_cached_value = generate_cluster_mass(REQUEST(cluster_redshift), _rng));

DEPENDENT_PARAM(cluster_num_satellites,
		_cached_value = generate_count( get_cluster_richness(REQUEST(cluster_mass), REQUEST(cluster_redshift)) - 1, _rng) );

// Galaxy level

DEPENDENT_PARAM(apparent_mag_vis,
		_cached_value = get_apparent_mag_vis(REQUEST(stellar_mass), REQUEST(redshift)));

DEPENDENT_PARAM(physical_size_bulge,
		_cached_value = generate_physical_size_bulge(REQUEST(galaxy_type), REQUEST(redshift),
				REQUEST(stellar_mass), _rng));

DEPENDENT_PARAM(physical_size_disk,
		_cached_value = generate_physical_size_disk(REQUEST(galaxy_type), REQUEST(redshift),
				REQUEST(stellar_mass), _rng));

DEPENDENT_PARAM(redshift,
		if(is_field_galaxy(REQUEST(galaxy_type)))
		{
			const DepFieldRedshift * p_redshift_pp = dynamic_cast<const DepFieldRedshift *>(_p_params);
			if(p_redshift_pp==nullptr)
			{
				_cached_value = _p_params->get_independently(_rng);
			}
			else
			{
				_cached_value = p_redshift_pp->get_dependently(
						_request_param(cluster_redshift_name)->get_p_params(),
						REQUEST(cluster_density), REQUEST(galaxy_density),
						_rng);
			}
		}
		else
		{
			_cached_value = REQUEST(cluster_redshift);
		})

DEPENDENT_PARAM(rotation,
		if(is_satellite_galaxy(REQUEST(galaxy_type)))
			_cached_value = generate_rotation( REQUEST(xp), REQUEST(yp), REQUEST(cluster_xp), REQUEST(cluster_yp),
		         REQUEST(morphology), REQUEST(stellar_mass), _rng  );
		else
			_cached_value = _p_params->get_independently(_rng);)

DEPENDENT_PARAM(rp,
		_cached_value = generate_rp(REQUEST(galaxy_type), REQUEST(cluster_mass), REQUEST(cluster_redshift),
				REQUEST(pixel_scale), _rng));

DEPENDENT_PARAM(shear_angle,
		_cached_value = generate_shear_angle(REQUEST(xp), REQUEST(yp), _rng));

DEPENDENT_PARAM(shear_magnitude,
		_cached_value = generate_shear_magnitude(REQUEST(xp), REQUEST(yp), REQUEST(redshift), _rng));

DEPENDENT_PARAM(stellar_mass,
		_cached_value = generate_stellar_mass(REQUEST(galaxy_type), REQUEST(redshift), _rng));

DEPENDENT_PARAM(tilt,
		if(is_satellite_galaxy(REQUEST(galaxy_type)))
			_cached_value = generate_tilt( REQUEST(xp), REQUEST(yp), REQUEST(cluster_xp), REQUEST(cluster_yp),
		         REQUEST(morphology), REQUEST(stellar_mass), _rng  );
		else
			_cached_value = _p_params->get_independently(_rng);)

DEPENDENT_PARAM(xp,
		if(is_field_galaxy(REQUEST(galaxy_type)))
			_cached_value = _p_params->get_independently(_rng);
		if(is_central_galaxy(REQUEST(galaxy_type)))
			_cached_value = REQUEST(cluster_xp);
		else
			_cached_value = generate_xp(REQUEST(rp), REQUEST(theta_sat),
					REQUEST(cluster_xp),  _rng));

DEPENDENT_PARAM(yp,
		if(is_field_galaxy(REQUEST(galaxy_type)))
			_cached_value = _p_params->get_independently(_rng);
		if(is_central_galaxy(REQUEST(galaxy_type)))
			_cached_value = REQUEST(cluster_yp);
		else
			_cached_value = generate_yp(REQUEST(rp), REQUEST(theta_sat),
				REQUEST(cluster_yp), _rng));

// Undef the macro
#undef DEPENDENT_PARAM
#undef REQUEST

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_PARAMS_HPP_
