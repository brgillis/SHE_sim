/**********************************************************************\
 @file GalaxyGroup.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_LEVELS_GALAXYGROUP_HPP_
#define SHE_SIM_GAL_PARAMS_LEVELS_GALAXYGROUP_HPP_

#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/default_values.h"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class GalaxyGroup: public ParamHierarchyLevel
{

public:
	GalaxyGroup(ParamHierarchyLevel * const & parent = nullptr,
			generation_level_map_t * const & p_generation_level_map = nullptr);
	virtual ~GalaxyGroup();

	/**
	 * Get the hierarchy level for this class.
	 * @return The hierachy level. 0 = highest, 1 = just below 0, etc.
	 */
	virtual int_t get_hierarchy_level() const override {return dv::galaxy_group_level;}

	// Methods to add children
#if(1)

	Image * add_galaxy();

	void add_galaxies(int_t const & N);

#endif

	virtual ParamHierarchyLevel * clone() const override;

};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_LEVELS_GALAXYGROUP_HPP_