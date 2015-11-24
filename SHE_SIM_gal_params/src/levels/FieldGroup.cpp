/**********************************************************************\
 @file FieldGroup.cpp
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

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/params_list.hpp"
#include "SHE_SIM_gal_params/levels/Field.hpp"
#include "SHE_SIM_gal_params/levels/FieldGroup.hpp"

namespace SHE_SIM
{

FieldGroup::FieldGroup(ParamHierarchyLevel * const & p_parent)
: ParamHierarchyLevel(p_parent,
		get_full_params_map(*this))
{
}

FieldGroup::~FieldGroup()
{
}

// Methods to add children
#if(1)

Field * FieldGroup::add_field()
{
	return static_cast<Field *>(ParamHierarchyLevel::spawn_child<Field>());
}

void FieldGroup::add_fields(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<Field>(N);
}

#endif

ParamHierarchyLevel * FieldGroup::clone() const
{
	return new FieldGroup(*this);
}

} // namespace SHE_SIM
