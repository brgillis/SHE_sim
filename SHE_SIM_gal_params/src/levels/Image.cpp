/**********************************************************************\
 @file Image.cpp
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
#include "SHE_SIM_gal_params/levels/Cluster.hpp"
#include "SHE_SIM_gal_params/levels/ClusterGroup.hpp"
#include "SHE_SIM_gal_params/levels/Field.hpp"
#include "SHE_SIM_gal_params/levels/FieldGroup.hpp"
#include "SHE_SIM_gal_params/levels/Image.hpp"

namespace SHE_SIM
{

Image::Image(ParamHierarchyLevel * const & p_parent)
: ParamHierarchyLevel(p_parent,
		get_full_params_map(*this))
{
}

Image::~Image()
{
}

// Methods to add children
#if(1)

ClusterGroup * Image::add_cluster_group()
{
	return static_cast<ClusterGroup *>(ParamHierarchyLevel::spawn_child<ClusterGroup>());
}

void Image::add_cluster_groups(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<ClusterGroup>(N);
}

Cluster * Image::add_cluster()
{
	return static_cast<Cluster *>(ParamHierarchyLevel::spawn_child<Cluster>());
}

void Image::add_clusters(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<Cluster>(N);
}

FieldGroup * Image::add_field_group()
{
	return static_cast<FieldGroup *>(ParamHierarchyLevel::spawn_child<FieldGroup>());
}

void Image::add_field_groups(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<FieldGroup>(N);
}

Field * Image::add_field()
{
	return static_cast<Field *>(ParamHierarchyLevel::spawn_child<Field>());
}

void Image::add_fields(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<Field>(N);
}

#endif

ParamHierarchyLevel * Image::clone() const
{
	return new Image(*this);
}

} // namespace SHE_SIM
