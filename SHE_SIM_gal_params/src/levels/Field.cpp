/**********************************************************************\
 @file Field.cpp
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

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/params_list.hpp"
#include "SHE_SIM_gal_params/levels/Field.hpp"
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"
#include "SHE_SIM_gal_params/levels/GalaxyGroup.hpp"

namespace SHE_SIM
{

Field::Field(ParamHierarchyLevel * const & p_parent)
: ParamHierarchyLevel(p_parent,
		get_full_params_map(*this))
{
}

Field::~Field()
{
}

// Methods to add children
#if(1)

GalaxyGroup * Field::add_galaxy_group()
{
	return static_cast<GalaxyGroup *>(ParamHierarchyLevel::spawn_child<GalaxyGroup>());
}

void Field::add_galaxy_groups(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<GalaxyGroup>(N);
}

Galaxy * Field::add_galaxy()
{
	return static_cast<Galaxy *>(ParamHierarchyLevel::spawn_child<Galaxy>());
}

void Field::add_galaxies(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<Galaxy>(N);
}

#endif

ParamHierarchyLevel * Field::clone() const
{
	return new Field(*this);
}

} // namespace SHE_SIM
