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

// Include all needed param param headers here
#include "SHE_SIM_gal_params/param_params/Calculated.hpp"
#include <SHE_SIM_gal_params/param_params/DepFieldRedshift.hpp>
#include <SHE_SIM_gal_params/param_params/IndClusterRedshift.hpp>
#include <SHE_SIM_gal_params/param_params/IndContRayleigh.hpp>
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/param_params/IndLogNormalMean.hpp"
#include "SHE_SIM_gal_params/param_params/IndTruncLogNormalMean.hpp"
#include "SHE_SIM_gal_params/param_params/IndUniform.hpp"

#include "SHE_SIM_gal_params/dependency_functions/cosmology.hpp"
#include "SHE_SIM_gal_params/dependency_functions/galaxy_redshift.hpp"
#include "SHE_SIM_gal_params/dependency_functions/galaxy_type.hpp"
#include "SHE_SIM_gal_params/dependency_functions/halos.hpp"
#include "SHE_SIM_gal_params/dependency_functions/morphology.hpp"
#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"

#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_declarations.hpp>

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
#define IMPLEMENT_PARAM(param, \
			            level, \
                        param_params,\
                        dependent_generation, \
						alt_dependent_generation) \
 \
const name_t param##_name = #param; \
 \
void param##_obj::_generate() \
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
		_cached_value = _p_params->get_independently(get_rng()); \
	} \
	else \
	{ \
		throw bad_mode_error(_p_params->get_mode_name()); \
	} \
} \
param##_obj::param##_obj( owner_t * const & p_owner) \
: ParamGenerator(p_owner) \
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

#define REQUEST(param) _request_param_value(param##_name)

#include "src/param_implementation_detail/high_level_param_implementations.hh"

#include "src/param_implementation_detail/mid_level_param_implementations.hh"

#include "src/param_implementation_detail/galaxy_level_param_implementations.hh"

} // namespace SHE_SIM
