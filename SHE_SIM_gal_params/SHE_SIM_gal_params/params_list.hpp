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

#include <memory>
#include <unordered_map>
#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"

#include "SHE_SIM_gal_params/params/ExposureTime.hpp"
#include "SHE_SIM_gal_params/params/MagVisInstZp.hpp"
#include "SHE_SIM_gal_params/params/MagVisZp.hpp"

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

inline ParamHierarchyLevel::params_t get_full_params_map(ParamHierarchyLevel & owner)
{
	ParamHierarchyLevel::params_t res;

	// Insert all params here
	insert_param<ExposureTime>(res, owner);
	insert_param<MagVisInstZp>(res, owner);
	insert_param<MagVisZp>(res, owner);

	return std::move(res);

} // params_t get_full_params_map()

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_PARAMS_LIST_HPP_
