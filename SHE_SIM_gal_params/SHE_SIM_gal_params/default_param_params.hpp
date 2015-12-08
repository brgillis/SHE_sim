/**********************************************************************\
 @file default_param_params.h
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

#ifndef SHE_SIM_GAL_PARAMS_DEFAULT_PARAM_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_DEFAULT_PARAM_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_values.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <memory>
#include <unordered_map>
#include <utility>

#include "SHE_SIM_gal_params/ParamParam.hpp"

// Include all needed param param headers here
#include "SHE_SIM_gal_params/param_params/Calculated.hpp"
#include <SHE_SIM_gal_params/param_params/IndContRayleigh.hpp>
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/param_params/IndLogNormalMean.hpp"
#include "SHE_SIM_gal_params/param_params/IndUniform.hpp"

namespace SHE_SIM {

extern const param_params_t default_param_params_map;
extern const generation_level_map_t default_generation_levels_map;

template<typename T_in, typename... Args>
void insert_default_param_param(param_params_t & res, const name_t & param_name, Args... args)
{
	typename param_params_t::mapped_type new_ptr(new T_in(args...));

	res.insert(std::make_pair(std::move(param_name),std::move(new_ptr)));
}

// Function to get a list of all params
inline param_params_t make_default_param_params_map()
{
	param_params_t res;

#define INSERT_FIXED_PARAM(param) insert_default_param_param<IndFixed>(res, param##_name, dv::param);
#define INSERT_LOGNORMAL_PARAM(param) insert_default_param_param<IndLogNormalMean>(res, param##_name, dv::param##_l10_mean, dv::param##_l10_stddev);
#define INSERT_UNIFORM_PARAM(param) insert_default_param_param<IndUniform>(res, param##_name, dv::param##_min, dv::param##_max);
#define INSERT_CONTRAYLEIGH_PARAM(param) insert_default_param_param<IndContRayleigh>(res, param##_name, dv::param##_sigma, dv::param##_max, dv::param##_p);
#define INSERT_CALCULATED_PARAM(param) insert_default_param_param<Calculated>(res, param##_name);

	// Insert all defaults  here

	// Survey level

	INSERT_FIXED_PARAM(gain);
	INSERT_FIXED_PARAM(num_images);
	INSERT_FIXED_PARAM(mag_i_inst_zp);
	INSERT_CALCULATED_PARAM(mag_i_zp);
	INSERT_FIXED_PARAM(mag_vis_inst_zp);
	INSERT_CALCULATED_PARAM(mag_vis_zp);
	INSERT_FIXED_PARAM(pixel_scale);
	INSERT_FIXED_PARAM(read_noise);


	// Image level

	INSERT_CALCULATED_PARAM(background_noise);
	INSERT_FIXED_PARAM(exp_time);
	INSERT_FIXED_PARAM(background_galaxy_density);
	INSERT_FIXED_PARAM(cluster_density);
	INSERT_FIXED_PARAM(field_galaxy_density);
	INSERT_CALCULATED_PARAM(image_area);
	INSERT_FIXED_PARAM(image_size_xp);
	INSERT_FIXED_PARAM(image_size_yp);
	INSERT_CALCULATED_PARAM(num_background_galaxies);
	INSERT_CALCULATED_PARAM(num_clusters);
	INSERT_FIXED_PARAM(num_fields);
	INSERT_CALCULATED_PARAM(num_stars);
	INSERT_FIXED_PARAM(star_density);
	INSERT_FIXED_PARAM(psf_params);
	INSERT_LOGNORMAL_PARAM(subtracted_background);
	INSERT_FIXED_PARAM(unsubtracted_background);


	// Cluster level

	INSERT_LOGNORMAL_PARAM(cluster_mass);
	INSERT_UNIFORM_PARAM(cluster_redshift);
	INSERT_UNIFORM_PARAM(cluster_xp);
	INSERT_UNIFORM_PARAM(cluster_yp);

	INSERT_CALCULATED_PARAM(cluster_num_satellites);

	// Field level

	INSERT_CALCULATED_PARAM(num_field_galaxies);

	// Galaxy level

	INSERT_CALCULATED_PARAM(apparent_mag_vis);
	INSERT_LOGNORMAL_PARAM(apparent_size);
	INSERT_CALCULATED_PARAM(galaxy_type);
	INSERT_UNIFORM_PARAM(morphology);
	INSERT_CALCULATED_PARAM(physical_size);
	INSERT_CALCULATED_PARAM(psf_model);
	INSERT_UNIFORM_PARAM(redshift);
	INSERT_UNIFORM_PARAM(rotation);
	INSERT_CALCULATED_PARAM(rp);
	INSERT_UNIFORM_PARAM(shear_angle);
	INSERT_CONTRAYLEIGH_PARAM(shear_magnitude);
	INSERT_CALCULATED_PARAM(stellar_mass);
	INSERT_UNIFORM_PARAM(theta_sat);
	INSERT_UNIFORM_PARAM(tilt);
	INSERT_UNIFORM_PARAM(xp);
	INSERT_UNIFORM_PARAM(yp);


	// GalaxyDither level

	INSERT_FIXED_PARAM(dither_xp_shift);
	INSERT_FIXED_PARAM(dither_yp_shift);

#undef INSERT_CALCULATED_PARAM
#undef INSERT_CONTRAYLEIGH_PARAM
#undef INSERT_FIXED_PARAM
#undef INSERT_LOGNORMAL_PARAM
#undef INSERT_UNIFORM_PARAM

	return res;

} // param_params_t get_full_param_params_map()

// Function to get a list of all params
inline generation_level_map_t make_default_generation_levels_map()
{
	generation_level_map_t res;

#define INSERT_LEVEL(param,level) res.insert(std::make_pair(param##_name,level_ptr_t(new level_t(level))));

	// Insert all defaults  here

	// Survey level

	INSERT_LEVEL(gain, dv::survey_level);
	INSERT_LEVEL(mag_i_inst_zp, dv::survey_level);
	INSERT_LEVEL(mag_vis_inst_zp, dv::survey_level);
	INSERT_LEVEL(num_images, dv::survey_level);
	INSERT_LEVEL(pixel_scale, dv::survey_level);
	INSERT_LEVEL(read_noise, dv::survey_level);

	// Image level

	INSERT_LEVEL(exp_time, dv::image_level);
	INSERT_LEVEL(background_galaxy_density, dv::image_level);
	INSERT_LEVEL(cluster_density, dv::image_level);
	INSERT_LEVEL(field_galaxy_density, dv::image_level);
	INSERT_LEVEL(image_size_xp, dv::image_level);
	INSERT_LEVEL(image_size_yp, dv::image_level);
	INSERT_LEVEL(num_fields, dv::image_level);
	INSERT_LEVEL(star_density, dv::image_level);
	INSERT_LEVEL(psf_params, dv::image_level);
	INSERT_LEVEL(subtracted_background, dv::image_level);
	INSERT_LEVEL(unsubtracted_background, dv::image_level);

	INSERT_LEVEL(background_noise, dv::image_level);
	INSERT_LEVEL(image_area, dv::image_level);
	INSERT_LEVEL(mag_i_zp, dv::image_level);
	INSERT_LEVEL(mag_vis_zp, dv::image_level);
	INSERT_LEVEL(num_background_galaxies, dv::image_level);
	INSERT_LEVEL(num_clusters, dv::image_level);
	INSERT_LEVEL(num_field_galaxies, dv::image_level);
	INSERT_LEVEL(num_stars, dv::image_level);

	// Cluster level

	INSERT_LEVEL(cluster_mass, dv::cluster_level);
	INSERT_LEVEL(cluster_redshift, dv::cluster_level);
	INSERT_LEVEL(cluster_xp, dv::cluster_level);
	INSERT_LEVEL(cluster_yp, dv::cluster_level);

	INSERT_LEVEL(cluster_num_satellites, dv::cluster_level);

	// Field level

	INSERT_LEVEL(num_field_galaxies, dv::field_level);

	// Galaxy level

	INSERT_LEVEL(apparent_mag_vis, dv::galaxy_level);
	INSERT_LEVEL(apparent_size, dv::galaxy_level);
	INSERT_LEVEL(galaxy_type, dv::galaxy_level);
	INSERT_LEVEL(morphology, dv::galaxy_level);
	INSERT_LEVEL(physical_size, dv::galaxy_level);
	INSERT_LEVEL(redshift, dv::galaxy_level);
	INSERT_LEVEL(rotation, dv::galaxy_level);
	INSERT_LEVEL(rp, dv::galaxy_level);
	INSERT_LEVEL(shear_angle, dv::galaxy_level);
	INSERT_LEVEL(shear_magnitude, dv::galaxy_level);
	INSERT_LEVEL(stellar_mass, dv::galaxy_level);
	INSERT_LEVEL(theta_sat, dv::galaxy_level);
	INSERT_LEVEL(tilt, dv::galaxy_level);
	INSERT_LEVEL(xp, dv::galaxy_level);
	INSERT_LEVEL(yp, dv::galaxy_level);


	// GalaxyDither level

	INSERT_LEVEL(dither_xp_shift, dv::galaxy_dither_level);
	INSERT_LEVEL(dither_yp_shift, dv::galaxy_dither_level);

#undef INSERT_LEVEL

	return res;

} // param_params_t get_full_param_params_map()

} // namespace SHE_SIM



#endif // SHE_SIM_GAL_PARAMS_DEFAULT_PARAM_PARAMS_HPP_

