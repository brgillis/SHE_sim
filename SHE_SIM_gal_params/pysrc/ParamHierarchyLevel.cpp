/**********************************************************************\
 @file name.cpp
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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

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

#define PHL_WRAPPER(name) \
struct name##Wrap : name, wrapper<name> \
{ \
    void wrapped_fill_children() { this->fill_children(); } \
 \
	std::vector<child_t *> wrapped_get_children( str const & type_name = "" ) \
	{ \
		return name::get_children(extract<name_t>(type_name)); \
	} \
 \
	child_t * wrapped_get_child(int const & i) { return name::get_child(i); } \
 \
	flt_t wrapped_get_generation_level(str const & name) \
	{ return name::get_generation_level(extract<name_t>(name)); } \
 \
	int_t wrapped_get_hierarchy_level() const { return this->get_hierarchy_level(); } \
 \
	name_t wrapped_get_name() const { return this->get_name(); } \
 \
	int_t wrapped_get_local_ID() const { return name::get_local_ID(); } \
 \
	flt_t wrapped_get_param_value(str const & name) { return name::get_param_value(extract<name_t>(name)); } \
 \
	void wrapped_set_param_params(str const & name, str const & param_type ) \
	{ name::set_param_params(extract<name_t>(name), extract<name_t>(param_type)); } \
	void wrapped_set_param_params(str const & name, str const & param_type, \
			flt_t const & p1 ) \
	{ name::set_param_params(extract<name_t>(name), extract<name_t>(param_type), p1); } \
	void wrapped_set_param_params(str const & name, str const & param_type, \
			flt_t const & p1, flt_t const & p2 ) \
	{ name::set_param_params(extract<name_t>(name), extract<name_t>(param_type), p1, p2); } \
	void wrapped_set_param_params(str const & name, str const & param_type, \
			flt_t const & p1, flt_t const & p2, flt_t const & p3 ) \
	{ name::set_param_params(extract<name_t>(name), extract<name_t>(param_type), p1, p2, p3); } \
	void wrapped_set_param_params(str const & name, str const & param_type, \
			flt_t const & p1, flt_t const & p2, flt_t const & p3, flt_t const & p4 ) \
	{ name::set_param_params(extract<name_t>(name), extract<name_t>(param_type), p1, p2, p3, p4); } \
 \
	int_t wrapped_get_seed() const { return name::get_seed(); } \
	void wrapped_set_seed( int_t const & seed ) { return name::set_seed(seed); } \
};

PHL_WRAPPER(ParamHierarchyLevel)
PHL_WRAPPER(Survey)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ParamHierarchyLevel_get_children_overloads, ParamHierarchyLevelWrap::wrapped_get_children, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ParamHierarchyLevel_set_param_param_overloads, ParamHierarchyLevelWrap::wrapped_set_param_params, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Survey_get_children_overloads, SurveyWrap::wrapped_get_children, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Survey_set_param_param_overloads, SurveyWrap::wrapped_set_param_params, 2, 6)

BOOST_PYTHON_MODULE(SHE_SIM)
{
    class_<std::vector<ParamHierarchyLevel *>>("ChildList")
        .def(vector_indexing_suite<std::vector<ParamHierarchyLevel *>>() );
    class_<std::vector<ImageGroup *>>("ImageGroupList")
        .def(vector_indexing_suite<std::vector<ImageGroup *>>() );
    class_<std::vector<Image *>>("ImageList")
        .def(vector_indexing_suite<std::vector<Image *>>() );
    class_<std::vector<ClusterGroup *>>("ClusterGroupList")
        .def(vector_indexing_suite<std::vector<ClusterGroup *>>() );
    class_<std::vector<Cluster *>>("ClusterList")
        .def(vector_indexing_suite<std::vector<Cluster *>>() );
    class_<std::vector<FieldGroup *>>("FieldGroupList")
        .def(vector_indexing_suite<std::vector<FieldGroup *>>() );
    class_<std::vector<Field *>>("FieldList")
        .def(vector_indexing_suite<std::vector<Field *>>() );
    class_<std::vector<GalaxyGroup *>>("GalaxyGroupList")
        .def(vector_indexing_suite<std::vector<GalaxyGroup *>>() );
    class_<std::vector<Galaxy *>>("GalaxyList")
        .def(vector_indexing_suite<std::vector<Galaxy *>>() );

    ParamHierarchyLevel::children_t const & (ParamHierarchyLevel::*gc0)() = &ParamHierarchyLevel::get_children;

#define PHL_DEFS(name) \
\
    .add_property("num_children", &name::num_children) \
\
    .def("clear_children", &name::clear_children) \
    .def("get_children", &name##Wrap::wrapped_get_children, name##_get_children_overloads( args("type_name") )) \
    .def("children", gc0, return_value_policy<reference_existing_object>()) \
    .def("fill_children", &name::fill_children) \
    .def("autofill_children", &name::autofill_children) \
\
	.def("get_child", &name##Wrap::wrapped_get_child, return_value_policy<reference_existing_object>() ) \
\
	.def("get_generation_level", &name##Wrap::wrapped_get_generation_level ) \
	.def("set_generation_level", &name::set_generation_level ) \
\
	.def("get_param_value", &name##Wrap::wrapped_get_param_value ) \
	.def("set_param_params", (void(name##Wrap::*)(str const &, str const &, \
			flt_t const &, flt_t const &, flt_t const &, flt_t const & ))0, name##_set_param_param_overloads() ) \
\
	.def("clear", &name::clear) \
\
    .add_property("hierarchy_level", &name##Wrap::wrapped_get_hierarchy_level) \
    .add_property("name", &name##Wrap::wrapped_get_name) \
    .add_property("local_ID", &name##Wrap::wrapped_get_local_ID) \
    .add_property("seed", &name##Wrap::wrapped_get_seed, &name##Wrap::wrapped_set_seed)

    class_<ParamHierarchyLevelWrap, boost::noncopyable>("name", no_init)

		PHL_DEFS(ParamHierarchyLevel)

        .enable_pickling();

	class_<SurveyWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Survey")

		PHL_DEFS(Survey)

		.def("add_image_group", &Survey::add_image_group, return_value_policy<reference_existing_object>())
		.def("add_image_groups", &Survey::add_image_groups)
		.def("get_image_groups", &Survey::get_image_groups)

		.def("add_image", &Survey::add_image, return_value_policy<reference_existing_object>())
		.def("add_images", &Survey::add_images)
		.def("get_images", &Survey::get_images)
		.def("fill_images", &Survey::fill_images)
		.def("autofill_images", &Survey::autofill_images)

		.enable_pickling();

	class_<ImageGroup, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("ImageGroup", no_init)
		.def("add_image", &ImageGroup::add_image, return_value_policy<reference_existing_object>())
		.def("add_images", &ImageGroup::add_images)
		.def("get_images", &ImageGroup::get_images)
		.enable_pickling();

	class_<Image, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Image", no_init)
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

	class_<ClusterGroup, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("ClusterGroup", no_init)
		.def("add_cluster", &ClusterGroup::add_cluster, return_value_policy<reference_existing_object>())
		.def("add_clusters", &ClusterGroup::add_clusters)
		.def("get_clusters", &ClusterGroup::get_clusters)

		.enable_pickling();

	class_<Cluster, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Cluster", no_init)

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

	class_<FieldGroup, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("FieldGroup", no_init)
		.def("add_field", &FieldGroup::add_field, return_value_policy<reference_existing_object>())
		.def("add_fields", &FieldGroup::add_fields)
		.def("get_fields", &FieldGroup::get_fields)

		.enable_pickling();

	class_<Field, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Field", no_init)

		.def("add_galaxy_group", &Field::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Field::add_galaxy_groups)
		.def("get_galaxy_groups", &Field::get_galaxy_groups)

		.def("add_galaxy", &Field::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Field::add_galaxies)
		.def("get_galaxies", &Field::get_galaxies)

		.enable_pickling();

	class_<GalaxyGroup, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("GalaxyGroup", no_init)
		.def("add_galaxy", &GalaxyGroup::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &GalaxyGroup::add_galaxies)
		.def("get_galaxies", &GalaxyGroup::get_galaxies)

		.enable_pickling();
}
