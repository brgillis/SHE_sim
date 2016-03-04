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

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/params_list.hpp"
#include <SHE_SIM_gal_params/param_declarations.hpp>
#include "SHE_SIM_gal_params/dependency_functions/galaxy_type.hpp"
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"

namespace SHE_SIM
{

Galaxy::Galaxy(ParamHierarchyLevel * const & p_parent)
: ParamHierarchyLevel(p_parent),
  _is_background(false)
{
}

void Galaxy::set_as_background_galaxy()
{
	set_param_params(galaxy_type_name,"fixed",field_galaxy_type);
	set_param_params(apparent_mag_vis_name,"alt_calculated");
	set_param_params(apparent_size_bulge_name,"alt_calculated");
	set_param_params(apparent_size_disk_name,"alt_calculated");
	_is_background = true;
}

void Galaxy::set_as_foreground_galaxy()
{
	set_param_params(apparent_mag_vis_name,"calculated");
	set_param_params(apparent_size_bulge_name,"calculated");
	set_param_params(apparent_size_disk_name,"calculated");
}

bool Galaxy::is_central_galaxy()
{
	return ::SHE_SIM::is_central_galaxy(get_param_value(galaxy_type_name));
}
bool Galaxy::is_field_galaxy()
{
	return ::SHE_SIM::is_field_galaxy(get_param_value(galaxy_type_name));
}
bool Galaxy::is_satellite_galaxy()
{
	return ::SHE_SIM::is_satellite_galaxy(get_param_value(galaxy_type_name));
}

} // namespace SHE_SIM
