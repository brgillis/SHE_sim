/**********************************************************************\
 @file Hierarchy_test.cpp
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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/param_names.h"
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/levels/ImageGroup.hpp"
#include "SHE_SIM_gal_params/levels/Image.hpp"

namespace SHE_SIM
{

struct hierarchy_fixture {

	Survey survey;

	const IndFixed mag_vis_inst_zp1 = IndFixed(25.0);
	const IndFixed mag_vis_inst_zp2 = IndFixed(26.0);

	const IndFixed exp_time0 = IndFixed(1);
	const IndFixed exp_time1 = IndFixed(1234.5);
	const IndFixed exp_time2 = IndFixed(2468.0);
	const IndFixed exp_time3 = IndFixed(3451.1);

	const flt_t expected_mag_vis_zp10 = mag_vis_inst_zp1.get_independently() + 2.5 * std::log10(exp_time0.get_independently());
	const flt_t expected_mag_vis_zp11 = mag_vis_inst_zp1.get_independently() + 2.5 * std::log10(exp_time1.get_independently());
	const flt_t expected_mag_vis_zp12 = mag_vis_inst_zp1.get_independently() + 2.5 * std::log10(exp_time2.get_independently());
	const flt_t expected_mag_vis_zp13 = mag_vis_inst_zp1.get_independently() + 2.5 * std::log10(exp_time3.get_independently());
	const flt_t expected_mag_vis_zp20 = mag_vis_inst_zp2.get_independently() + 2.5 * std::log10(exp_time0.get_independently());
	const flt_t expected_mag_vis_zp21 = mag_vis_inst_zp2.get_independently() + 2.5 * std::log10(exp_time1.get_independently());
	const flt_t expected_mag_vis_zp22 = mag_vis_inst_zp2.get_independently() + 2.5 * std::log10(exp_time2.get_independently());
	const flt_t expected_mag_vis_zp23 = mag_vis_inst_zp2.get_independently() + 2.5 * std::log10(exp_time3.get_independently());

};


BOOST_AUTO_TEST_SUITE (Hierarchy_Test)

BOOST_FIXTURE_TEST_CASE(test_hierarchy, hierarchy_fixture) {

	// Setup
	survey.set_generation_level(exp_time_name, dv::survey_level);
	survey.set_generation_level(mag_vis_inst_zp_name, dv::survey_level);
	survey.set_generation_level(mag_vis_zp_name, dv::survey_level);

	survey.set_p_param_params(exp_time_name,&exp_time0);
	survey.set_p_param_params(mag_vis_inst_zp_name,&mag_vis_inst_zp1);

	ImageGroup & image_group1 = *survey.add_image_group();
	survey.add_image_groups(2);

	// Check we have three image groups now
	BOOST_CHECK_EQUAL(int_t(survey.get_children().size()),3);

	ImageGroup & image_group2 = *static_cast<ImageGroup *>(survey.get_child(1));
	ImageGroup & image_group3 = *static_cast<ImageGroup *>(survey.get_child(2));

	Image & image11 = *image_group1.add_image();

	image_group1.set_p_param_params(exp_time_name,&exp_time1);
	image_group2.set_p_param_params(exp_time_name,&exp_time2);
	image_group3.set_p_param_params(exp_time_name,&exp_time3);

	Image & image12 = *image_group1.add_image();
	Image & image21 = *image_group2.add_image();

	// Check that each image group gets the zp from the survey values
	BOOST_CHECK_CLOSE(image_group1.get_param_value(mag_vis_zp_name),expected_mag_vis_zp10,1e-9);
	BOOST_CHECK_CLOSE(image_group2.get_param_value(mag_vis_zp_name),expected_mag_vis_zp10,1e-9);
	BOOST_CHECK_CLOSE(image_group3.get_param_value(mag_vis_zp_name),expected_mag_vis_zp10,1e-9);
	BOOST_CHECK_CLOSE(image11.get_param_value(mag_vis_zp_name),expected_mag_vis_zp10,1e-9);
	BOOST_CHECK_CLOSE(image12.get_param_value(mag_vis_zp_name),expected_mag_vis_zp10,1e-9);
	BOOST_CHECK_CLOSE(image21.get_param_value(mag_vis_zp_name),expected_mag_vis_zp10,1e-9);

	survey.set_generation_level(exp_time_name, dv::image_level);
	survey.set_generation_level(mag_vis_zp_name, dv::image_level);

	// Check that each image group gets the zp from its own values now
	BOOST_CHECK_CLOSE(image_group1.get_param_value(mag_vis_zp_name),expected_mag_vis_zp11,1e-9);
	BOOST_CHECK_CLOSE(image_group2.get_param_value(mag_vis_zp_name),expected_mag_vis_zp12,1e-9);
	BOOST_CHECK_CLOSE(image_group3.get_param_value(mag_vis_zp_name),expected_mag_vis_zp13,1e-9);
	BOOST_CHECK_CLOSE(image11.get_param_value(mag_vis_zp_name),expected_mag_vis_zp11,1e-9);
	BOOST_CHECK_CLOSE(image12.get_param_value(mag_vis_zp_name),expected_mag_vis_zp11,1e-9);
	BOOST_CHECK_CLOSE(image21.get_param_value(mag_vis_zp_name),expected_mag_vis_zp12,1e-9);

	// Change the survey's inst zp
	survey.set_p_param_params(mag_vis_inst_zp_name, &mag_vis_inst_zp2);

	// Check that each image group gets the correct new zp now
	BOOST_CHECK_CLOSE(image_group1.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	BOOST_CHECK_CLOSE(image_group2.get_param_value(mag_vis_zp_name),expected_mag_vis_zp22,1e-9);
	BOOST_CHECK_CLOSE(image_group3.get_param_value(mag_vis_zp_name),expected_mag_vis_zp23,1e-9);
	BOOST_CHECK_CLOSE(image11.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	BOOST_CHECK_CLOSE(image12.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	BOOST_CHECK_CLOSE(image21.get_param_value(mag_vis_zp_name),expected_mag_vis_zp22,1e-9);

	// Change the survey's generation levels
	survey.set_generation_level(exp_time_name, dv::survey_level);
	survey.set_generation_level(mag_vis_zp_name, dv::survey_level);

	// Check that each image group gets the zp from the survey values now
	BOOST_CHECK_CLOSE(image_group1.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image_group2.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image_group3.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image11.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image12.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image21.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);

	// Change image_group1's generation levels
	image_group1.set_generation_level(exp_time_name, dv::image_level);
	image_group1.set_generation_level(mag_vis_zp_name, dv::image_level);

	// Check
	BOOST_CHECK_CLOSE(image_group1.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	BOOST_CHECK_CLOSE(image_group2.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image_group3.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);
	BOOST_CHECK_CLOSE(image11.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	BOOST_CHECK_CLOSE(image12.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	BOOST_CHECK_CLOSE(image21.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);

	// Check if new parameters will inherit from their parents
	Image & image13 = *image_group1.add_image();
	BOOST_CHECK_CLOSE(image13.get_param_value(mag_vis_zp_name),expected_mag_vis_zp21,1e-9);
	Image & image22 = *image_group2.add_image();
	BOOST_CHECK_CLOSE(image22.get_param_value(mag_vis_zp_name),expected_mag_vis_zp20,1e-9);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
