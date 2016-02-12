/**********************************************************************\
 @file alt_dependent_params.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_ALT_DEPENDENT_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_ALT_DEPENDENT_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/dependency_functions/cosmology.hpp"
#include "SHE_SIM_gal_params/dependency_functions/galaxy_type.hpp"
#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "IceBRG_main/math/random/random_functions.hpp"

namespace SHE_SIM
{

// Define a macro for setting up params

#define ALT_DEPENDENT_PARAM( param_name, dependent_generation, alt_dependent_generation ) \
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
		else if(_p_params->get_mode()==ParamParam::ALT_DEPENDENT) \
		{ \
			alt_dependent_generation; \
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

// Cluster level

// Galaxy level

ALT_DEPENDENT_PARAM(apparent_size_bulge,
		_cached_value = get_angle_from_distance(REQUEST(physical_size_bulge), REQUEST(redshift)),
		_cached_value = generate_apparent_size_bulge(REQUEST(apparent_mag_vis), _rng));

ALT_DEPENDENT_PARAM(apparent_size_disk,
		_cached_value = get_angle_from_distance(REQUEST(physical_size_bulge), REQUEST(redshift)),
		_cached_value = generate_apparent_size_disk(REQUEST(apparent_mag_vis), _rng));

ALT_DEPENDENT_PARAM(bulge_fraction,
		_cached_value = generate_bulge_fraction(REQUEST(galaxy_type), REQUEST(redshift),
				REQUEST(stellar_mass), REQUEST(morphology), _rng),
		_cached_value = generate_bulge_fraction(REQUEST(apparent_mag_vis), REQUEST(morphology), _rng));

ALT_DEPENDENT_PARAM(morphology,
		_cached_value = generate_morphology(REQUEST(galaxy_type), REQUEST(redshift), REQUEST(stellar_mass), _rng),
		_cached_value = generate_morphology(REQUEST(apparent_mag_vis), _rng));

// Undef the macro
#undef ALT_DEPENDENT_PARAM
#undef REQUEST

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_ALT_DEPENDENT_PARAMS_HPP_
