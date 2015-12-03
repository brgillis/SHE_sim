/**********************************************************************\
 @file Galaxy.cpp
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
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"
#include "SHE_SIM_gal_params/levels/GalaxyDither.hpp"

namespace SHE_SIM
{

Galaxy::Galaxy(ParamHierarchyLevel * const & p_parent)
: ParamHierarchyLevel(p_parent,
		get_full_params_map(*this))
{
}

Galaxy::~Galaxy()
{
}

// Methods to add children
#if(1)

GalaxyDither * Galaxy::add_galaxy_dither()
{
	return static_cast<GalaxyDither *>(ParamHierarchyLevel::spawn_child<GalaxyDither>());
}

void Galaxy::add_galaxy_dithers(int_t const & N)
{
	return ParamHierarchyLevel::spawn_children<GalaxyDither>(N);
}

#endif

void Galaxy::set_as_background_galaxy()
{
	set_param_params(galaxy_type_name,"fixed",field_galaxy_type);
	set_param_params(apparent_mag_vis_name,"alt_calculated");
	set_param_params(apparent_size_name,"alt_calculated");
}

ParamHierarchyLevel * Galaxy::clone() const
{
	return new Galaxy(*this);
}

} // namespace SHE_SIM
