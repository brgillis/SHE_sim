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
		_params = default_param_params_map.at(name()).get(); \
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
INDEPENDENT_PARAM(mag_i_inst_zp);
INDEPENDENT_PARAM(mag_vis_inst_zp);
INDEPENDENT_PARAM(pixel_scale);
INDEPENDENT_PARAM(read_noise);

// Image level
INDEPENDENT_PARAM(background_galaxy_density);
INDEPENDENT_PARAM(cluster_density);
INDEPENDENT_PARAM(field_galaxy_density);
INDEPENDENT_PARAM(galaxy_type);
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

INDEPENDENT_PARAM(theta_sat);

// GalaxyDither level
INDEPENDENT_PARAM(dither_xp_shift);
INDEPENDENT_PARAM(dither_yp_shift);

// Undef the macro
#undef INDEPENDENT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_PARAMS_HPP_
