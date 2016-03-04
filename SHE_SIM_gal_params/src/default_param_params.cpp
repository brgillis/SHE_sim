/**********************************************************************\
 @file random_functions.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015  Bryan R. Gillis

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

#include <utility>

#include "SHE_SIM_gal_params/default_values.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/param_names.hpp"

#include <SHE_SIM_gal_params/params/alt_dependent_params.hpp>
#include <SHE_SIM_gal_params/params/dependent_params.hpp>
#include <SHE_SIM_gal_params/params/independent_params.hpp>

// Include all needed param param headers here
#include "SHE_SIM_gal_params/param_params/Calculated.hpp"
#include <SHE_SIM_gal_params/param_params/DepFieldRedshift.hpp>
#include <SHE_SIM_gal_params/param_params/IndClusterRedshift.hpp>
#include <SHE_SIM_gal_params/param_params/IndContRayleigh.hpp>
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/param_params/IndLogNormalMean.hpp"
#include "SHE_SIM_gal_params/param_params/IndTruncLogNormalMean.hpp"
#include "SHE_SIM_gal_params/param_params/IndUniform.hpp"

#include <SHE_SIM_gal_params/default_param_params.hpp>

namespace SHE_SIM {

// Implement default maps

params_t default_params_map;
param_params_t default_param_params_map;
generation_level_map_t default_generation_levels_map;

// Functions to help insert objects into maps

template<typename T_in, typename... Args>
void insert_default_param_param(const name_t & param_name,
		level_t const & gen_level,
		Args... args)
{
	typename param_params_t::mapped_type new_ptr(new T_in(args...));

	default_param_params_map.insert(std::make_pair(std::move(param_name),std::move(new_ptr)));
	default_generation_levels_map.insert(std::make_pair(param_name,
			level_ptr_t(new level_t(gen_level))));
}

// Macros to simplify adding param params to the default maps, using the attach by initialization
// idiom.
#define INSERT_PARAM(param,level,param_params) \
 \
const name_t param##_name = #param; \
 \
struct param##_initializer \
{ \
	param##_initializer() \
	{ \
		default_generation_levels_map.insert(std::make_pair(param##_name, \
				level_ptr_t(new level_t(level)))); \
	     \
		typename param_params_t::mapped_type new_ptr(new param_params); \
		default_param_params_map.insert(std::make_pair(param##_name, \
				std::move(new_ptr))); \
		 \
		default_params_map.insert(std::make_pair(param##_name, \
			param_ptr_t(new param##_obj(nullptr)))); \
		 \
	} \
}; \
\
param##_initializer param##_initializer_instance;

// Survey level

INSERT_PARAM(num_images, dv::survey_level, IndFixed(dv::num_images));
INSERT_PARAM(pixel_scale, dv::survey_level, IndFixed(dv::pixel_scale));

// Image level

INSERT_PARAM(exp_time, dv::image_level, IndFixed(dv::exp_time));
INSERT_PARAM(cluster_density, dv::image_level, IndFixed(dv::cluster_density));
INSERT_PARAM(galaxy_density, dv::image_level, IndFixed(dv::galaxy_density));
INSERT_PARAM(image_area, dv::image_level, Calculated);
INSERT_PARAM(image_size_xp, dv::image_level, IndFixed(dv::image_size_xp));
INSERT_PARAM(image_size_yp, dv::image_level, IndFixed(dv::image_size_xp));
INSERT_PARAM(num_clusters, dv::image_level, Calculated);
INSERT_PARAM(num_fields, dv::image_level, IndFixed(dv::num_fields));
INSERT_PARAM(subtracted_background, dv::image_level,
		IndLogNormalMean(dv::subtracted_background_l10_mean,dv::subtracted_background_l10_stddev));
INSERT_PARAM(unsubtracted_background, dv::image_level, IndFixed(dv::unsubtracted_background));

// Cluster level

INSERT_PARAM(cluster_mass, dv::cluster_level, Calculated);
INSERT_PARAM(cluster_redshift, dv::cluster_level, IndClusterRedshift(
		dv::cluster_redshift_enhancement,
		dv::cluster_redshift_min,
		dv::cluster_redshift_max));
INSERT_PARAM(cluster_xp, dv::cluster_level, IndUniform(dv::cluster_xp_min,
														dv::cluster_xp_max));
INSERT_PARAM(cluster_yp, dv::cluster_level, IndUniform(dv::cluster_yp_min,
														dv::cluster_yp_max));

INSERT_PARAM(cluster_num_satellites, dv::cluster_level, Calculated);

// Field level

INSERT_PARAM(num_field_galaxies, dv::field_level, Calculated);

// Galaxy level

INSERT_PARAM(absolute_mag_vis, dv::galaxy_level, Calculated);
INSERT_PARAM(apparent_mag_vis, dv::galaxy_level, Calculated);
INSERT_PARAM(apparent_size_bulge, dv::galaxy_level, Calculated);
INSERT_PARAM(apparent_size_disk, dv::galaxy_level, Calculated);
INSERT_PARAM(bulge_class, dv::galaxy_level, Calculated);
INSERT_PARAM(bulge_fraction, dv::galaxy_level, Calculated);
INSERT_PARAM(galaxy_type, dv::galaxy_level, IndFixed(dv::galaxy_type));
INSERT_PARAM(physical_size_bulge, dv::galaxy_level, Calculated);
INSERT_PARAM(physical_size_disk, dv::galaxy_level, Calculated);
INSERT_PARAM(redshift, dv::galaxy_level, DepFieldRedshift(dv::galaxy_redshift_enhancement,
		dv::galaxy_redshift_min, dv::galaxy_redshift_max));
INSERT_PARAM(rotation, dv::galaxy_level, IndUniform(dv::rotation_min, dv::rotation_max));
INSERT_PARAM(rp, dv::galaxy_level, Calculated);
INSERT_PARAM(sersic_index, dv::galaxy_level, Calculated);
INSERT_PARAM(shear_angle, dv::galaxy_level, IndUniform(dv::shear_angle_min,
															dv::shear_angle_max));
INSERT_PARAM(shear_magnitude, dv::galaxy_level, IndContRayleigh(dv::shear_magnitude_sigma,
		dv::shear_magnitude_max, \
		dv::shear_magnitude_p));
INSERT_PARAM(stellar_mass, dv::galaxy_level, Calculated);
INSERT_PARAM(theta_sat, dv::galaxy_level, IndUniform(dv::theta_sat_min, dv::theta_sat_max));
INSERT_PARAM(tilt, dv::galaxy_level, IndUniform(dv::tilt_min, dv::tilt_max));
INSERT_PARAM(xp, dv::galaxy_level, IndUniform(dv::xp_min, dv::yp_max));
INSERT_PARAM(yp, dv::galaxy_level, IndUniform(dv::yp_min, dv::yp_max));

} // namespace SHE_SIM
