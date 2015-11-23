/**********************************************************************\
 @file param_params_list.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_LIST_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_LIST_HPP_

#include <memory>
#include <unordered_map>
#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

// Include all param param headers here
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/param_params/IndContRayleigh.hpp"
#include "SHE_SIM_gal_params/param_params/IndGaussian.hpp"
#include <SHE_SIM_gal_params/param_params/IndLogNormalPeak.hpp>
#include <SHE_SIM_gal_params/param_params/IndLogNormalMean.hpp>
#include "SHE_SIM_gal_params/param_params/IndPoisson.hpp"
#include "SHE_SIM_gal_params/param_params/IndRayleigh.hpp"
#include "SHE_SIM_gal_params/param_params/IndTruncGaussian.hpp"
#include "SHE_SIM_gal_params/param_params/IndTruncRayleigh.hpp"
#include "SHE_SIM_gal_params/param_params/IndUniform.hpp"

namespace SHE_SIM {

extern const ParamHierarchyLevel::param_params_t param_params_map;

template<typename T_in, typename T_map>
void insert_param_param(T_map & res)
{
	typename T_map::mapped_type new_ptr(new T_in);
	auto name(new_ptr->name());

	res.insert(std::make_pair(std::move(name),std::move(new_ptr)));
}

// Function to get a list of all params
inline ParamHierarchyLevel::param_params_t get_full_param_params_map()
{
	ParamHierarchyLevel::param_params_t res;

	// Insert all param_params here
	insert_param_param<IndFixed>(res);
	insert_param_param<IndContRayleigh>(res);
	insert_param_param<IndGaussian>(res);
	insert_param_param<IndLogNormalPeak>(res);
	insert_param_param<IndLogNormalMean>(res);
	insert_param_param<IndPoisson>(res);
	insert_param_param<IndRayleigh>(res);
	insert_param_param<IndTruncGaussian>(res);
	insert_param_param<IndTruncRayleigh>(res);
	insert_param_param<IndUniform>(res);

	return res;

} // param_params_t get_full_param_params_map()

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_LIST_HPP_
