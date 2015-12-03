/**********************************************************************\
 @file Cluster.cpp
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
#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/dependency_functions/galaxy_type.hpp"
#include "SHE_SIM_gal_params/levels/Cluster.hpp"
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"
#include "SHE_SIM_gal_params/levels/GalaxyGroup.hpp"

namespace SHE_SIM
{

Cluster::Cluster(ParamHierarchyLevel * const & p_parent)
: ParamHierarchyLevel(p_parent,
		get_full_params_map(*this))
{
}

Cluster::~Cluster()
{
}

// Methods to add children
#if(1)

GalaxyGroup * Cluster::add_galaxy_group()
{
	return static_cast<GalaxyGroup *>(ParamHierarchyLevel::spawn_child<GalaxyGroup>());
}

void Cluster::add_galaxy_groups(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<GalaxyGroup>(N);
}

Galaxy * Cluster::add_galaxy()
{
	return static_cast<Galaxy *>(ParamHierarchyLevel::spawn_child<Galaxy>());
}

void Cluster::add_galaxies(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<Galaxy>(N);
}

Galaxy * Cluster::add_central_galaxy()
{
	Galaxy * gal = add_galaxy();
	gal->set_param_params(galaxy_type_name,"fixed",central_galaxy_type);

	return gal;
}

Galaxy * Cluster::add_satellite_galaxy()
{
	Galaxy * gal = add_galaxy();
	gal->set_param_params(galaxy_type_name,"fixed",satellite_galaxy_type);

	return gal;
}

void Cluster::add_satellite_galaxies(int_t const & N)
{
	for(int i=0; i<N; ++i) add_satellite_galaxy();
}

#endif

// Methods to automatically add children
#if(1)

void Cluster::fill_children()
{
	fill_galaxies();
}

void Cluster::fill_galaxies()
{
	add_central_galaxy();
	add_satellite_galaxies( round_int( get_param_value( cluster_num_satellites_name ) ) );
}

#endif

ParamHierarchyLevel * Cluster::clone() const
{
	return new Cluster(*this);
}

} // namespace SHE_SIM
