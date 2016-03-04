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
#include <utility>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"


namespace SHE_SIM {

class ParamGenerator;
class ParamHierarchyLevel;

extern params_t default_params_map;

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

	for( auto const & name_and_ptr : default_params_map )
	{
		name_t name = name_and_ptr.first;
		param_ptr_t new_ptr(name_and_ptr.second->clone());
		new_ptr->set_owner(owner);

		res.insert(std::make_pair(std::move(name),std::move(new_ptr)));
	}

	return res;

} // params_t get_full_params_map()

#undef INSERT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_PARAMS_LIST_HPP_
