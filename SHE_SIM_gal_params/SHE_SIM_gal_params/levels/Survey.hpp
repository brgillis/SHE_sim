/**********************************************************************\
 @file Survey.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_LEVELS_SURVEY_HPP_
#define SHE_SIM_GAL_PARAMS_LEVELS_SURVEY_HPP_

#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/default_values.h"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"

namespace SHE_SIM
{

// Forward-declare children
class ImageGroup;
class Image;

/**
 * TODO Auto-generated comment stub
 */
class Survey: public ParamHierarchyLevel
{
private:

	generation_level_map_t _survey_generation_level_map;

public:
	Survey();
	virtual ~Survey();

	/**
	 * Get the hierarchy level for this class.
	 * @return The hierachy level. 0 = highest, 1 = just below 0, etc.
	 */
	virtual int_t get_hierarchy_level() const override {return dv::survey_level;}

	// Methods for interacting with the generation level map
#if(1)

	const generation_level_map_t & get_survey_generation_level_map() const noexcept;

	void set_survey_generation_level_map(
			const generation_level_map_t & survey_generation_level_map);

	void set_survey_generation_level_map(
			generation_level_map_t && survey_generation_level_map);

	void set_generation_level( const name_t & name, const int_t & generation_level );

#endif

	// Methods to add children
#if(1)

	ImageGroup * add_image_group();

	void add_image_groups(int_t const & N);

	Image * add_image();

	void add_images(int_t const & N);

#endif

	virtual ParamHierarchyLevel * clone() const override;

};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_LEVELS_SURVEY_HPP_
