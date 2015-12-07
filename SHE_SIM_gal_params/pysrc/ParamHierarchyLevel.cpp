/**********************************************************************\
 @file ParamHierarchyLevel.cpp
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

#include <boost/python.hpp>

#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/levels/ClusterGroup.hpp"
#include "SHE_SIM_gal_params/levels/Cluster.hpp"
#include "SHE_SIM_gal_params/levels/FieldGroup.hpp"
#include "SHE_SIM_gal_params/levels/Field.hpp"
#include "SHE_SIM_gal_params/levels/GalaxyGroup.hpp"
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"
#include "SHE_SIM_gal_params/levels/ImageGroup.hpp"
#include "SHE_SIM_gal_params/levels/Image.hpp"
#include "SHE_SIM_gal_params/levels/Survey.hpp"

using namespace boost::python;
using namespace SHE_SIM;

// Set up wrappers for functions we can't use directly for one reason or another
struct ParamHierarchyLevelWrap : ParamHierarchyLevel, wrapper<ParamHierarchyLevel>
{
    void fill_children()
    {
        if (override fill_children = this->get_override("fill_children"))
            fill_children();
        ParamHierarchyLevel::fill_children();
    }

    void default_fill_children() { return this->ParamHierarchyLevel::fill_children(); }

	std::vector<child_t *> get_children( name_t const & type_name = "" )
	{
		return get_children(type_name);
	}

	child_t * get_child(int const & i) { return ParamHierarchyLevel::get_child(i); }

	flt_t get_generation_level(name_t const & name) { return ParamHierarchyLevel::get_generation_level(name); }

	int_t get_hierarchy_level() const { return this->get_override("get_hierarchy_level")(); }

	name_t get_name() const { return this->get_override("fill_children")(); }

	int_t get_local_ID() const { return ParamHierarchyLevel::get_local_ID(); }

	flt_t get_param_value(name_t const & name) { return ParamHierarchyLevel::get_param_value(name); }

	void set_param_params(name_t const & name, str_t const & param_type )
	{ ParamHierarchyLevel::set_param_params(name, param_type); }
	void set_param_params(name_t const & name, str_t const & param_type,
			flt_t const & p1 )
	{ ParamHierarchyLevel::set_param_params(name, param_type, p1); }
	void set_param_params(name_t const & name, str_t const & param_type,
			flt_t const & p1, flt_t const & p2 )
	{ ParamHierarchyLevel::set_param_params(name, param_type, p1, p2); }
	void set_param_params(name_t const & name, str_t const & param_type,
			flt_t const & p1, flt_t const & p2, flt_t const & p3 )
	{ ParamHierarchyLevel::set_param_params(name, param_type, p1, p2, p3); }
	void set_param_params(name_t const & name, str_t const & param_type,
			flt_t const & p1, flt_t const & p2, flt_t const & p3, flt_t const & p4 )
	{ ParamHierarchyLevel::set_param_params(name, param_type, p1, p2, p3, p4); }

	int_t get_seed() const { return ParamHierarchyLevel::get_seed(); }
	void set_seed( int_t const & seed ) { return ParamHierarchyLevel::set_seed(seed); }
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PHL_get_children_overloads, ParamHierarchyLevelWrap::get_children, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PHL_set_param_param_overloads, ParamHierarchyLevelWrap::set_param_params, 2, 6)

BOOST_PYTHON_MODULE(SHE_SIM)
{
	ParamHierarchyLevel::children_t const & (ParamHierarchyLevel::*gc0)() = &ParamHierarchyLevel::get_children;

    class_<ParamHierarchyLevelWrap, boost::noncopyable>("ParamHierarchyLevel", no_init)

        .def("num_children", &ParamHierarchyLevel::num_children)

        .def("clear_children", &ParamHierarchyLevel::clear_children)
        .def("get_children", &ParamHierarchyLevelWrap::get_children, PHL_get_children_overloads( args("type_name") ))
        .def("children", gc0, return_value_policy<reference_existing_object>())
        .def("fill_children", &ParamHierarchyLevel::fill_children, &ParamHierarchyLevelWrap::default_fill_children)
        .def("autofill_children", &ParamHierarchyLevel::autofill_children)

		.def("get_child", &ParamHierarchyLevelWrap::get_child, return_value_policy<reference_existing_object>() )

		.def("get_generation_level", &ParamHierarchyLevelWrap::get_generation_level )
		.def("set_generation_level", &ParamHierarchyLevel::set_generation_level )

		.def("get_param_value", &ParamHierarchyLevelWrap::get_param_value )
		.def("set_param_params", (void(ParamHierarchyLevelWrap::*)(name_t const &, name_t const &,
				flt_t const &, flt_t const &, flt_t const &, flt_t const & ))0, PHL_set_param_param_overloads() )

		.def("clear", &ParamHierarchyLevel::clear)

        .add_property("hierarchy_level", &ParamHierarchyLevelWrap::get_hierarchy_level)
        .add_property("name", &ParamHierarchyLevelWrap::get_name)
	    .add_property("local_ID", &ParamHierarchyLevelWrap::get_local_ID)
	    .add_property("seed", &ParamHierarchyLevelWrap::get_seed, &ParamHierarchyLevelWrap::set_seed)

        .enable_pickling();

	class_<Survey, bases<ParamHierarchyLevel> >("Survey")
		.def("add_image_group", &Survey::add_image_group, return_value_policy<reference_existing_object>())
		.def("add_image_groups", &Survey::add_image_groups)
		.def("get_image_groups", &Survey::get_image_groups)

		.def("add_image", &Survey::add_image, return_value_policy<reference_existing_object>())
		.def("add_images", &Survey::add_images)
		.def("get_images", &Survey::get_images)
		.def("fill_images", &Survey::fill_images)
		.def("autofill_images", &Survey::autofill_images)

		.enable_pickling();

	class_<ImageGroup, bases<ParamHierarchyLevel> >("ImageGroup", no_init)
		.def("add_image", &ImageGroup::add_image, return_value_policy<reference_existing_object>())
		.def("add_images", &ImageGroup::add_images)
		.def("get_images", &ImageGroup::get_images)
		.enable_pickling();

	class_<Image, bases<ParamHierarchyLevel> >("Image", no_init)
		.def("add_cluster_group", &Image::add_cluster_group, return_value_policy<reference_existing_object>())
		.def("add_cluster_groups", &Image::add_cluster_groups)
		.def("get_cluster_groups", &Image::get_cluster_groups)

		.def("add_cluster", &Image::add_cluster, return_value_policy<reference_existing_object>())
		.def("add_clusters", &Image::add_clusters)
		.def("get_clusters", &Image::get_clusters)
		.def("fill_clusters", &Image::fill_clusters)
		.def("autofill_clusters", &Image::autofill_clusters)

		.def("add_field_group", &Image::add_field_group, return_value_policy<reference_existing_object>())
		.def("add_field_groups", &Image::add_field_groups)
		.def("get_field_groups", &Image::get_field_groups)

		.def("add_field", &Image::add_field, return_value_policy<reference_existing_object>())
		.def("add_fields", &Image::add_fields)
		.def("get_fields", &Image::get_fields)
		.def("fill_field", &Image::fill_field)
		.def("autofill_field", &Image::autofill_field)

		.def("add_galaxy_group", &Image::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Image::add_galaxy_groups)
		.def("get_galaxy_groups", &Image::get_galaxy_groups)

		.def("add_galaxy", &Image::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Image::add_galaxies)
		.def("get_galaxies", &Image::get_galaxies)

		.def("add_background_galaxy", &Image::add_background_galaxy, return_value_policy<reference_existing_object>())
		.def("add_background_galaxies", &Image::add_background_galaxies)
		.def("get_background_galaxies", &Image::get_background_galaxies)
		.def("fill_background_galaxies", &Image::fill_background_galaxies)
		.def("autofill_background_galaxies", &Image::autofill_background_galaxies)

		.def("add_foreground_galaxy", &Image::add_foreground_galaxy, return_value_policy<reference_existing_object>())
		.def("add_foreground_galaxies", &Image::add_foreground_galaxies)
		.def("get_foreground_galaxies", &Image::get_foreground_galaxies)

		.enable_pickling();

	class_<ClusterGroup, bases<ParamHierarchyLevel> >("ClusterGroup", no_init)
		.def("add_cluster", &ClusterGroup::add_cluster, return_value_policy<reference_existing_object>())
		.def("add_clusters", &ClusterGroup::add_clusters)
		.def("get_clusters", &ClusterGroup::get_clusters)

		.enable_pickling();

	class_<Cluster, bases<ParamHierarchyLevel> >("Cluster", no_init)

		.def("add_galaxy_group", &Cluster::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Cluster::add_galaxy_groups)
		.def("get_galaxy_groups", &Cluster::get_galaxy_groups)

		.def("add_galaxy", &Cluster::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Cluster::add_galaxies)
		.def("get_galaxies", &Cluster::get_galaxies)

		.def("add_central_galaxy", &Cluster::add_central_galaxy, return_value_policy<reference_existing_object>())
		.def("get_central_galaxy", &Cluster::get_central_galaxy, return_value_policy<reference_existing_object>())

		.def("add_satellite_galaxy", &Cluster::add_satellite_galaxy, return_value_policy<reference_existing_object>())
		.def("add_satellite_galaxies", &Cluster::add_satellite_galaxies)
		.def("get_satellite_galaxies", &Cluster::get_satellite_galaxies)

		.enable_pickling();

	class_<FieldGroup, bases<ParamHierarchyLevel> >("FieldGroup", no_init)
		.def("add_field", &FieldGroup::add_field, return_value_policy<reference_existing_object>())
		.def("add_fields", &FieldGroup::add_fields)
		.def("get_fields", &FieldGroup::get_fields)

		.enable_pickling();

	class_<Field, bases<ParamHierarchyLevel> >("Field", no_init)

		.def("add_galaxy_group", &Field::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Field::add_galaxy_groups)
		.def("get_galaxy_groups", &Field::get_galaxy_groups)

		.def("add_galaxy", &Field::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Field::add_galaxies)
		.def("get_galaxies", &Field::get_galaxies)

		.enable_pickling();

	class_<GalaxyGroup, bases<ParamHierarchyLevel> >("GalaxyGroup", no_init)
		.def("add_galaxy", &GalaxyGroup::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &GalaxyGroup::add_galaxies)
		.def("get_galaxies", &GalaxyGroup::get_galaxies)

		.enable_pickling();
}
