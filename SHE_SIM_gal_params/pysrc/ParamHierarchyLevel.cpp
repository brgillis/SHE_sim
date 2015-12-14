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

// Define a coersion function
template <typename Tout, typename Tin>
Tout coerce( Tin const & vin )
{
	Tout res;
	res.reserve(vin.size());

	for( auto & v : vin )
	{
		res.push_back(static_cast<typename Tout::value_type>(v));
	}

	return res;
}

// Forward-declare wrappers

class ParamHierarchyLevelWrap;
class ImageGroupWrap;
class ImageWrap;
class ClusterGroupWrap;
class ClusterWrap;
class FieldGroupWrap;
class FieldWrap;
class GalaxyGroupWrap;
class GalaxyWrap;

// Set up wrappers for functions we can't use directly for one reason or another

#define PHL_WRAPPER(name) \
	struct name##Wrap : name, wrapper<name> \
	{ \
		void wrapped_fill_children() { this->fill_children(); } \
		 \
		std::vector<ParamHierarchyLevelWrap *> wrapped_get_children( ) \
		{ \
			auto children = name::get_children(""); \
			return coerce<std::vector<ParamHierarchyLevelWrap *>>(children); \
		} \
	 \
		std::vector<ParamHierarchyLevelWrap *> wrapped_get_children( str const & type_name ) \
		{ \
			auto children = name::get_children(extract<name_t>(type_name)); \
			return coerce<std::vector<ParamHierarchyLevelWrap *>>(children); \
		} \
	 \
		std::vector<ImageGroupWrap *> wrapped_get_image_groups( ) \
		{ \
			auto children = name::get_image_groups(); \
			return coerce<std::vector<ImageGroupWrap *>>(children); \
		} \
	 \
		std::vector<ImageWrap *> wrapped_get_images( ) \
		{ \
			auto children = name::get_images(); \
			return coerce<std::vector<ImageWrap *>>(children); \
		} \
	 \
		std::vector<ClusterGroupWrap *> wrapped_get_cluster_groups( ) \
		{ \
			auto children = name::get_cluster_groups(); \
			return coerce<std::vector<ClusterGroupWrap *>>(children); \
		} \
	 \
		std::vector<ClusterWrap *> wrapped_get_clusters( ) \
		{ \
			auto children = name::get_clusters(); \
			return coerce<std::vector<ClusterWrap *>>(children); \
		} \
	 \
		std::vector<FieldGroupWrap *> wrapped_get_field_groups( ) \
		{ \
			auto children = name::get_field_groups(); \
			return coerce<std::vector<FieldGroupWrap *>>(children); \
		} \
	 \
		std::vector<FieldWrap *> wrapped_get_fields( ) \
		{ \
			auto children = name::get_fields(); \
			return coerce<std::vector<FieldWrap *>>(children); \
		} \
	 \
		std::vector<GalaxyGroupWrap *> wrapped_get_galaxy_groups( ) \
		{ \
			auto children = name::get_galaxy_groups(); \
			return coerce<std::vector<GalaxyGroupWrap *>>(children); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_galaxies( ) \
		{ \
			auto children = name::get_galaxies(); \
			return coerce<std::vector<GalaxyWrap *>>(children); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_background_galaxies( ) \
		{ \
			auto children = name::get_background_galaxies(); \
			return coerce<std::vector<GalaxyWrap *>>(children); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_foreground_galaxies( ) \
		{ \
			auto children = name::get_foreground_galaxies(); \
			return coerce<std::vector<GalaxyWrap *>>(children); \
		} \
	 \
		GalaxyWrap * wrapped_get_central_galaxy( ) \
		{ \
			return static_cast<GalaxyWrap *>(name::get_central_galaxy()); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_satellite_galaxies( ) \
		{ \
			auto children = name::get_satellite_galaxies(); \
			return coerce<std::vector<GalaxyWrap *>>(children); \
		} \
	 \
		std::vector<ParamHierarchyLevelWrap *> wrapped_get_descendants( ) \
		{ \
			auto descendants = name::get_descendants(""); \
			return coerce<std::vector<ParamHierarchyLevelWrap *>>(descendants); \
		} \
	 \
		std::vector<ParamHierarchyLevelWrap *> wrapped_get_descendants( str const & type_name ) \
		{ \
			auto descendants = name::get_descendants(extract<name_t>(type_name)); \
			return coerce<std::vector<ParamHierarchyLevelWrap *>>(descendants); \
		} \
	 \
		std::vector<ImageGroupWrap *> wrapped_get_image_group_descendants( ) \
		{ \
			auto descendants = name::get_image_group_descendants(); \
			return coerce<std::vector<ImageGroupWrap *>>(descendants); \
		} \
	 \
		std::vector<ImageWrap *> wrapped_get_image_descendants( ) \
		{ \
			auto descendants = name::get_image_descendants(); \
			return coerce<std::vector<ImageWrap *>>(descendants); \
		} \
	 \
		std::vector<ClusterGroupWrap *> wrapped_get_cluster_group_descendants( ) \
		{ \
			auto descendants = name::get_cluster_group_descendants(); \
			return coerce<std::vector<ClusterGroupWrap *>>(descendants); \
		} \
	 \
		std::vector<ClusterWrap *> wrapped_get_cluster_descendants( ) \
		{ \
			auto descendants = name::get_cluster_descendants(); \
			return coerce<std::vector<ClusterWrap *>>(descendants); \
		} \
	 \
		std::vector<FieldGroupWrap *> wrapped_get_field_group_descendants( ) \
		{ \
			auto descendants = name::get_field_group_descendants(); \
			return coerce<std::vector<FieldGroupWrap *>>(descendants); \
		} \
	 \
		std::vector<FieldWrap *> wrapped_get_field_descendants( ) \
		{ \
			auto descendants = name::get_field_descendants(); \
			return coerce<std::vector<FieldWrap *>>(descendants); \
		} \
	 \
		std::vector<GalaxyGroupWrap *> wrapped_get_galaxy_group_descendants( ) \
		{ \
			auto descendants = name::get_galaxy_group_descendants(); \
			return coerce<std::vector<GalaxyGroupWrap *>>(descendants); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_galaxy_descendants( ) \
		{ \
			auto descendants = name::get_galaxy_descendants(); \
			return coerce<std::vector<GalaxyWrap *>>(descendants); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_background_galaxy_descendants( ) \
		{ \
			auto descendants = name::get_background_galaxy_descendants(); \
			return coerce<std::vector<GalaxyWrap *>>(descendants); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_foreground_galaxy_descendants( ) \
		{ \
			auto descendants = name::get_foreground_galaxy_descendants(); \
			return coerce<std::vector<GalaxyWrap *>>(descendants); \
		} \
	 \
	 	 std::vector<GalaxyWrap *> wrapped_get_central_galaxy_descendants( ) \
		{ \
			auto descendants = name::get_central_galaxy_descendants(); \
			return coerce<std::vector<GalaxyWrap *>>(descendants); \
		} \
	 \
		std::vector<GalaxyWrap *> wrapped_get_satellite_galaxy_descendants( ) \
		{ \
			auto descendants = name::get_satellite_galaxy_descendants(); \
			return coerce<std::vector<GalaxyWrap *>>(descendants); \
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
	}; \
	 \
	std::vector<ParamHierarchyLevelWrap *> (name##Wrap::*name##_gc0)(void) = &name##Wrap::wrapped_get_children; \
	std::vector<ParamHierarchyLevelWrap *> (name##Wrap::*name##_gc1)(str const &) = &name##Wrap::wrapped_get_children; \
	 \
	void (name##Wrap::*name##_spp0)(str const &, str const &) = &name##Wrap::wrapped_set_param_params; \
	void (name##Wrap::*name##_spp1)(str const &, str const &, \
			flt_t const & p1 ) = &name##Wrap::wrapped_set_param_params; \
	void (name##Wrap::*name##_spp2)(str const &, str const &, \
			flt_t const & p1, flt_t const & p2 ) = &name##Wrap::wrapped_set_param_params; \
	void (name##Wrap::*name##_spp3)(str const &, str const &, \
			flt_t const & p1, flt_t const & p2, flt_t const & p3 ) = &name##Wrap::wrapped_set_param_params; \
	void (name##Wrap::*name##_spp4)(str const &, str const &, \
			flt_t const & p1, flt_t const & p2, flt_t const & p3, flt_t const & p4 ) = &name##Wrap::wrapped_set_param_params;


PHL_WRAPPER(Galaxy) // Implement first due to the special get_central_galaxy() function

PHL_WRAPPER(ParamHierarchyLevel)
PHL_WRAPPER(Survey)
PHL_WRAPPER(ImageGroup)
PHL_WRAPPER(Image)
PHL_WRAPPER(ClusterGroup)
PHL_WRAPPER(Cluster)
PHL_WRAPPER(FieldGroup)
PHL_WRAPPER(Field)
PHL_WRAPPER(GalaxyGroup)

BOOST_PYTHON_MODULE(SHE_SIM)
{
    class_<std::vector<ParamHierarchyLevelWrap *>>("ChildList")
        .def(vector_indexing_suite<std::vector<ParamHierarchyLevelWrap *>>() );
    class_<std::vector<ImageGroupWrap *>>("ImageGroupList")
        .def(vector_indexing_suite<std::vector<ImageGroupWrap *>>() );
    class_<std::vector<ImageWrap *>>("ImageList")
        .def(vector_indexing_suite<std::vector<ImageWrap *>>() );
    class_<std::vector<ClusterGroupWrap *>>("ClusterGroupList")
        .def(vector_indexing_suite<std::vector<ClusterGroupWrap *>>() );
    class_<std::vector<ClusterWrap *>>("ClusterList")
        .def(vector_indexing_suite<std::vector<ClusterWrap *>>() );
    class_<std::vector<FieldGroupWrap *>>("FieldGroupList")
        .def(vector_indexing_suite<std::vector<FieldGroupWrap *>>() );
    class_<std::vector<FieldWrap *>>("FieldList")
        .def(vector_indexing_suite<std::vector<FieldWrap *>>() );
    class_<std::vector<GalaxyGroupWrap *>>("GalaxyGroupList")
        .def(vector_indexing_suite<std::vector<GalaxyGroupWrap *>>() );
    class_<std::vector<GalaxyWrap *>>("GalaxyList")
        .def(vector_indexing_suite<std::vector<GalaxyWrap *>>() );

#define PHL_DEFS(name) \
\
    .add_property("num_children", &name::num_children) \
\
    .def("clear_children", &name::clear_children) \
    .def("get_children", name##_gc0 ) \
    .def("get_children", name##_gc1 ) \
    .def("fill_children", &name::fill_children) \
    .def("autofill_children", &name::autofill_children) \
\
	.def("get_image_groups", &name##Wrap::wrapped_get_image_groups ) \
	.def("get_images", &name##Wrap::wrapped_get_images ) \
	.def("get_cluster_groups", &name##Wrap::wrapped_get_cluster_groups ) \
	.def("get_clusters", &name##Wrap::wrapped_get_clusters ) \
	.def("get_field_groups", &name##Wrap::wrapped_get_field_groups ) \
	.def("get_fields", &name##Wrap::wrapped_get_fields ) \
	.def("get_galaxy_groups", &name##Wrap::wrapped_get_galaxy_groups ) \
	.def("get_galaxies", &name##Wrap::wrapped_get_galaxies ) \
	.def("get_background_galaxies", &name##Wrap::wrapped_get_background_galaxies ) \
	.def("get_foreground_galaxies", &name##Wrap::wrapped_get_foreground_galaxies ) \
	.def("get_central_galaxy", &name##Wrap::wrapped_get_central_galaxy, return_value_policy<reference_existing_object>() ) \
	.def("get_satellite_galaxies", &name##Wrap::wrapped_get_satellite_galaxies ) \
\
	.def("get_image_group_descendants", &name##Wrap::wrapped_get_image_group_descendants ) \
	.def("get_image_descendants", &name##Wrap::wrapped_get_image_descendants ) \
	.def("get_cluster_group_descendants", &name##Wrap::wrapped_get_cluster_group_descendants ) \
	.def("get_cluster_descendants", &name##Wrap::wrapped_get_cluster_descendants ) \
	.def("get_field_group_descendants", &name##Wrap::wrapped_get_field_group_descendants ) \
	.def("get_field_descendants", &name##Wrap::wrapped_get_field_descendants ) \
	.def("get_galaxy_group_descendants", &name##Wrap::wrapped_get_galaxy_group_descendants ) \
	.def("get_galaxy_descendants", &name##Wrap::wrapped_get_galaxy_descendants ) \
	.def("get_background_galaxy_descendants", &name##Wrap::wrapped_get_background_galaxy_descendants ) \
	.def("get_foreground_galaxy_descendants", &name##Wrap::wrapped_get_foreground_galaxy_descendants ) \
	.def("get_central_galaxy_descendants", &name##Wrap::wrapped_get_foreground_galaxy_descendants ) \
	.def("get_satellite_galaxy_descendants", &name##Wrap::wrapped_get_satellite_galaxy_descendants ) \
\
	.def("get_generation_level", &name##Wrap::wrapped_get_generation_level ) \
	.def("set_generation_level", &name::set_generation_level ) \
\
	.def("get_param_value", &name##Wrap::wrapped_get_param_value ) \
	.def("set_param_params", name##_spp0 ) \
	.def("set_param_params", name##_spp1 ) \
	.def("set_param_params", name##_spp2 ) \
	.def("set_param_params", name##_spp3 ) \
	.def("set_param_params", name##_spp4 ) \
	.def("generate_parameters", &name::generate_parameters ) \
\
	.def("clear", &name::clear) \
\
    .add_property("hierarchy_level", &name##Wrap::wrapped_get_hierarchy_level) \
    .add_property("name", &name##Wrap::wrapped_get_name) \
    .add_property("local_ID", &name##Wrap::wrapped_get_local_ID) \
    .add_property("full_ID", &name::get_full_ID) \
    .add_property("seed", &name##Wrap::wrapped_get_seed, &name##Wrap::wrapped_set_seed) \
    .add_property("full_seed", &name::get_full_seed)

    class_<ParamHierarchyLevelWrap, boost::noncopyable>("name", no_init)

		PHL_DEFS(ParamHierarchyLevel)

        .enable_pickling();

	class_<SurveyWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Survey")

		PHL_DEFS(Survey)

		.def("add_image_group", &Survey::add_image_group, return_value_policy<reference_existing_object>())
		.def("add_image_groups", &Survey::add_image_groups)

		.def("add_image", &Survey::add_image, return_value_policy<reference_existing_object>())
		.def("add_images", &Survey::add_images)
		.def("fill_images", &Survey::fill_images)
		.def("autofill_images", &Survey::autofill_images)

		.enable_pickling();

	class_<Survey, bases<SurveyWrap>, boost::noncopyable >("SurveyBase").enable_pickling();

	class_<ImageGroupWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("ImageGroup", no_init)

		PHL_DEFS(ImageGroup)

		.def("add_image", &ImageGroup::add_image, return_value_policy<reference_existing_object>())
		.def("add_images", &ImageGroup::add_images)

		.enable_pickling();

	class_<ImageWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("ImageT", no_init)

		PHL_DEFS(Image)

		.def("add_cluster_group", &Image::add_cluster_group, return_value_policy<reference_existing_object>())
		.def("add_cluster_groups", &Image::add_cluster_groups)

		.def("add_cluster", &Image::add_cluster, return_value_policy<reference_existing_object>())
		.def("add_clusters", &Image::add_clusters)
		.def("fill_clusters", &Image::fill_clusters)
		.def("autofill_clusters", &Image::autofill_clusters)

		.def("add_field_group", &Image::add_field_group, return_value_policy<reference_existing_object>())
		.def("add_field_groups", &Image::add_field_groups)

		.def("add_field", &Image::add_field, return_value_policy<reference_existing_object>())
		.def("add_fields", &Image::add_fields)
		.def("fill_field", &Image::fill_field)
		.def("autofill_field", &Image::autofill_field)

		.def("add_galaxy_group", &Image::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Image::add_galaxy_groups)

		.def("add_galaxy", &Image::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Image::add_galaxies)

		.def("add_background_galaxy", &Image::add_background_galaxy, return_value_policy<reference_existing_object>())
		.def("add_background_galaxies", &Image::add_background_galaxies)
		.def("fill_background_galaxies", &Image::fill_background_galaxies)
		.def("autofill_background_galaxies", &Image::autofill_background_galaxies)

		.def("add_foreground_galaxy", &Image::add_foreground_galaxy, return_value_policy<reference_existing_object>())
		.def("add_foreground_galaxies", &Image::add_foreground_galaxies)

		.enable_pickling();

	class_<ClusterGroupWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("ClusterGroup", no_init)

		PHL_DEFS(ClusterGroup)

		.def("add_cluster", &ClusterGroup::add_cluster, return_value_policy<reference_existing_object>())
		.def("add_clusters", &ClusterGroup::add_clusters)

		.enable_pickling();

	class_<ClusterWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Cluster", no_init)

		PHL_DEFS(Cluster)

		.def("add_galaxy_group", &Cluster::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Cluster::add_galaxy_groups)

		.def("add_galaxy", &Cluster::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Cluster::add_galaxies)

		.def("add_central_galaxy", &Cluster::add_central_galaxy, return_value_policy<reference_existing_object>())

		.def("add_satellite_galaxy", &Cluster::add_satellite_galaxy, return_value_policy<reference_existing_object>())
		.def("add_satellite_galaxies", &Cluster::add_satellite_galaxies)

		.enable_pickling();

	class_<FieldGroupWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("FieldGroup", no_init)

		PHL_DEFS(FieldGroup)

		.def("add_field", &FieldGroup::add_field, return_value_policy<reference_existing_object>())
		.def("add_fields", &FieldGroup::add_fields)

		.enable_pickling();

	class_<FieldWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Field", no_init)

		PHL_DEFS(Field)

		.def("add_galaxy_group", &Field::add_galaxy_group, return_value_policy<reference_existing_object>())
		.def("add_galaxy_groups", &Field::add_galaxy_groups)

		.def("add_galaxy", &Field::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &Field::add_galaxies)

		.enable_pickling();

	class_<GalaxyGroupWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("GalaxyGroup", no_init)

		PHL_DEFS(GalaxyGroup)

		.def("add_galaxy", &GalaxyGroup::add_galaxy, return_value_policy<reference_existing_object>())
		.def("add_galaxies", &GalaxyGroup::add_galaxies)

		.enable_pickling();

	class_<GalaxyWrap, bases<ParamHierarchyLevelWrap>, boost::noncopyable >("Galaxy", no_init)

		PHL_DEFS(Galaxy)

		.enable_pickling();
}
