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

#include <memory>
#include <unordered_map>
#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/default_values.h"
#include "SHE_SIM_gal_params/param_names.h"
#include "SHE_SIM_gal_params/ParamParam.hpp"

// Include all needed param param headers here
#include "SHE_SIM_gal_params/param_params/Calculated.hpp"
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"

namespace SHE_SIM {

extern const param_params_t default_param_params_map;

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

	// Insert all defaults  here
	insert_default_param_param<IndFixed>(res, exp_time_name, dv::exp_time);
	insert_default_param_param<IndFixed>(res, mag_vis_inst_zp_name, dv::mag_vis_inst_zp);
	insert_default_param_param<IndFixed>(res, mag_i_inst_zp_name, dv::mag_i_inst_zp);
	insert_default_param_param<Calculated>(res, mag_vis_zp_name);
	insert_default_param_param<Calculated>(res, mag_i_zp_name);

	return res;

} // param_params_t get_full_param_params_map()

} // namespace SHE_SIM



#endif // SHE_SIM_GAL_PARAMS_DEFAULT_PARAM_PARAMS_HPP_

