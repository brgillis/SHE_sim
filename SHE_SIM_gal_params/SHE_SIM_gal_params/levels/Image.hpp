/**********************************************************************\
 @file Image.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_LEVELS_IMAGE_HPP_
#define SHE_SIM_GAL_PARAMS_LEVELS_IMAGE_HPP_

#include <vector>

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_values.hpp>
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/level_names.hpp"

namespace SHE_SIM
{

// Forward-declare children
class ClusterGroup;
class Cluster;
class FieldGroup;
class Field;
class GalaxyGroup;
class Galaxy;

/**
 * TODO Auto-generated comment stub
 */
class Image: public ParamHierarchyLevel
{

public:
	Image(ParamHierarchyLevel * const & parent = nullptr);
	virtual ~Image() {}

	/**
	 * Get the hierarchy level for this class.
	 * @return The hierachy level. 0 = highest, 1 = just below 0, etc.
	 */
	virtual int_t get_hierarchy_level() const override {return dv::image_level;}

	virtual name_t get_name() const override {return image_name;}

	// Methods to add children
#if(1)

	ClusterGroup * add_cluster_group();

	void add_cluster_groups(int_t const & N);

	Cluster * add_cluster();

	void add_clusters(int_t const & N);

	FieldGroup * add_field_group();

	void add_field_groups(int_t const & N);

	Field * add_field();

	void add_fields(int_t const & N);

	GalaxyGroup * add_galaxy_group();

	void add_galaxy_groups(int_t const & N);

	Galaxy * add_galaxy();

	void add_galaxies(int_t const & N);

	Galaxy * add_background_galaxy();

	void add_background_galaxies(int_t const & N);

	Galaxy * add_foreground_galaxy();

	void add_foreground_galaxies(int_t const & N);

#endif

	// Methods to automatically add children
#if(1)

	virtual void fill_children() override;

	void fill_clusters();
	void autofill_clusters();

	void fill_field();
	void autofill_field();

	void fill_background_galaxies();
	void autofill_background_galaxies();

#endif

	// Methods to get children of specific types
#if(1)

	std::vector<ClusterGroup *> get_cluster_groups();

	std::vector<Cluster *> get_clusters();

	std::vector<FieldGroup *> get_field_groups();

	std::vector<Field *> get_fields();

	std::vector<GalaxyGroup *> get_galaxy_groups();

	std::vector<Galaxy *> get_galaxies();

	std::vector<Galaxy *> get_background_galaxies();

	std::vector<Galaxy *> get_foreground_galaxies();

#endif

	virtual ParamHierarchyLevel * clone() const override { return new Image(*this); }

};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_LEVELS_IMAGE_HPP_
