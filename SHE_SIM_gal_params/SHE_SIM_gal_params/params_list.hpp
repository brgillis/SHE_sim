/**********************************************************************\
 @file params_list.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_LIST_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_LIST_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/params/alt_dependent_params.hpp>
#include <SHE_SIM_gal_params/params/dependent_params.hpp>
#include <SHE_SIM_gal_params/params/independent_params.hpp>
#include <memory>
#include <unordered_map>
#include <utility>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"


namespace SHE_SIM {

class ParamGenerator;
class ParamHierarchyLevel;

// Function to get a list of all params

template<typename T_in, typename T_map>
void insert_param(T_map & res, ParamHierarchyLevel & owner)
{
	typename T_map::mapped_type new_ptr(new T_in(owner));
	auto name(new_ptr->name());

	res.insert(std::make_pair(std::move(name),std::move(new_ptr)));
}

#define INSERT_PARAM(param) insert_param<param##_obj>(res,owner);

inline params_t get_full_params_map(ParamHierarchyLevel & owner)
{
	params_t res;

	// Survey level
	INSERT_PARAM(exp_time);
	INSERT_PARAM(gain);
	INSERT_PARAM(num_images);
	INSERT_PARAM(mag_i_inst_zp);
	INSERT_PARAM(mag_i_zp);
	INSERT_PARAM(mag_vis_inst_zp);
	INSERT_PARAM(mag_vis_zp);
	INSERT_PARAM(pixel_scale);
	INSERT_PARAM(read_noise);


	// Image level
	INSERT_PARAM(background_galaxy_density);
	INSERT_PARAM(background_noise);
	INSERT_PARAM(cluster_density);
	INSERT_PARAM(galaxy_density);
	INSERT_PARAM(image_area);
	INSERT_PARAM(image_size_xp);
	INSERT_PARAM(image_size_yp);
	INSERT_PARAM(num_background_galaxies);
	INSERT_PARAM(num_clusters);
	INSERT_PARAM(num_fields);
	INSERT_PARAM(num_stars);
	INSERT_PARAM(psf_params);
	INSERT_PARAM(star_density);
	INSERT_PARAM(subtracted_background);
	INSERT_PARAM(unsubtracted_background);


	// Cluster level
	INSERT_PARAM(cluster_mass);
	INSERT_PARAM(cluster_redshift);
	INSERT_PARAM(cluster_num_satellites);
	INSERT_PARAM(cluster_xp);
	INSERT_PARAM(cluster_yp);

	// Field level
	INSERT_PARAM(num_field_galaxies);

	// Galaxy-level params here

	INSERT_PARAM(apparent_mag_vis);
	INSERT_PARAM(apparent_size_bulge);
	INSERT_PARAM(apparent_size_disk);
	INSERT_PARAM(bulge_fraction);
	INSERT_PARAM(galaxy_type);
	INSERT_PARAM(morphology);
	INSERT_PARAM(physical_size_bulge);
	INSERT_PARAM(physical_size_disk);
	INSERT_PARAM(redshift);
	INSERT_PARAM(rotation);
	INSERT_PARAM(rp);
	INSERT_PARAM(shear_angle);
	INSERT_PARAM(shear_magnitude);
	INSERT_PARAM(stellar_mass);
	INSERT_PARAM(theta_sat);
	INSERT_PARAM(tilt);
	INSERT_PARAM(xp);
	INSERT_PARAM(yp);

	return res;

} // params_t get_full_params_map()

#undef INSERT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_PARAMS_LIST_HPP_
