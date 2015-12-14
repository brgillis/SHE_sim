/**********************************************************************\
 @file independent_params.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"

namespace SHE_SIM
{

// Define a macro for each param

#define INDEPENDENT_PARAM(param_name) \
class param_name##_obj : public ParamGenerator \
{ \
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

// Define each param

// Survey level
INDEPENDENT_PARAM(exp_time);
INDEPENDENT_PARAM(gain);
INDEPENDENT_PARAM(num_images);
INDEPENDENT_PARAM(mag_i_inst_zp);
INDEPENDENT_PARAM(mag_vis_inst_zp);
INDEPENDENT_PARAM(pixel_scale);
INDEPENDENT_PARAM(read_noise);

// Image level
INDEPENDENT_PARAM(background_galaxy_density);
INDEPENDENT_PARAM(cluster_density);
INDEPENDENT_PARAM(field_galaxy_density);
INDEPENDENT_PARAM(num_fields);
INDEPENDENT_PARAM(psf_params);
INDEPENDENT_PARAM(star_density);
INDEPENDENT_PARAM(subtracted_background);
INDEPENDENT_PARAM(unsubtracted_background);
INDEPENDENT_PARAM(image_size_xp);
INDEPENDENT_PARAM(image_size_yp);

// Cluster level
INDEPENDENT_PARAM(cluster_redshift);
INDEPENDENT_PARAM(cluster_xp);
INDEPENDENT_PARAM(cluster_yp);

// Galaxy level

INDEPENDENT_PARAM(galaxy_type);
INDEPENDENT_PARAM(theta_sat);

// Undef the macro
#undef INDEPENDENT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_PARAMS_HPP_
