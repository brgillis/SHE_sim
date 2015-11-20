/**********************************************************************\
 @file Survey.cpp
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <utility>

#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/params_list.hpp"

namespace SHE_SIM
{

Survey::Survey()
: ParamHierarchyLevel(nullptr,
		&_survey_generation_level_map,
		get_full_params_map(*this))
{
}

Survey::~Survey()
{
}

const generation_level_map_t & Survey::get_survey_generation_level_map() const noexcept
{
	return _survey_generation_level_map;
}

void Survey::set_survey_generation_level_map(
		const generation_level_map_t & survey_generation_level_map)
{
	_survey_generation_level_map = survey_generation_level_map;
}

void Survey::set_survey_generation_level_map(
		generation_level_map_t && survey_generation_level_map)
{
	_survey_generation_level_map = std::move(survey_generation_level_map);
}

void Survey::set_generation_level( const Survey::param_name_t & name, const int_t & generation_level )
{
	_survey_generation_level_map[name] = generation_level;
}

ParamHierarchyLevel * Survey::clone() const
{
	return new Survey(*this);
}

} // namespace SHE_SIM
