/**********************************************************************\
 @file PHL_ID_test.cpp
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

#include <vector>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/param_declarations.hpp>
#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/levels/ImageGroup.hpp"
#include "SHE_SIM_gal_params/levels/Image.hpp"

namespace SHE_SIM
{

struct PHL_ID_fixture {

	Survey survey;

	int_t expected_survey_local_ID = 0;
	int_t expected_ig1_local_ID = 0;
	int_t expected_ig2_local_ID = 1;
	int_t expected_ig3_local_ID = 2;
	int_t expected_i11_local_ID = 0;
	int_t expected_i12_local_ID = 1;
	int_t expected_i21_local_ID = 0;

	std::vector<int_t> expected_survey_ID_seq = {expected_survey_local_ID};
	std::vector<int_t> expected_ig1_ID_seq = {expected_survey_local_ID,expected_ig1_local_ID};
	std::vector<int_t> expected_ig2_ID_seq = {expected_survey_local_ID,expected_ig2_local_ID};
	std::vector<int_t> expected_ig3_ID_seq = {expected_survey_local_ID,expected_ig3_local_ID};
	std::vector<int_t> expected_i11_ID_seq = {expected_survey_local_ID,expected_ig1_local_ID,expected_i11_local_ID};
	std::vector<int_t> expected_i12_ID_seq = {expected_survey_local_ID,expected_ig1_local_ID,expected_i12_local_ID};
	std::vector<int_t> expected_i21_ID_seq = {expected_survey_local_ID,expected_ig2_local_ID,expected_i21_local_ID};

};


BOOST_AUTO_TEST_SUITE (PHL_ID_Test)

BOOST_FIXTURE_TEST_CASE(test_PHL_ID, PHL_ID_fixture) {

	// Check the survey's ID
	BOOST_CHECK_EQUAL(survey.get_local_ID(), expected_survey_local_ID);

	// Add image groups to the survey
	ImageGroup & image_group1 = *survey.add_image_group();
	survey.add_image_groups(2);
	ImageGroup & image_group2 = *static_cast<ImageGroup *>(survey.get_child(1));
	ImageGroup & image_group3 = *static_cast<ImageGroup *>(survey.get_child(2));

	// Check that each image group has the proper ID
	BOOST_CHECK_EQUAL(image_group1.get_local_ID(), expected_ig1_local_ID);
	BOOST_CHECK_EQUAL(image_group2.get_local_ID(), expected_ig2_local_ID);
	BOOST_CHECK_EQUAL(image_group3.get_local_ID(), expected_ig3_local_ID);

	// Add images
	Image & image11 = *image_group1.add_image();
	Image & image12 = *image_group1.add_image();
	Image & image21 = *image_group2.add_image();

	// Check that each image has the proper ID
	BOOST_CHECK_EQUAL(image11.get_local_ID(), expected_i11_local_ID);
	BOOST_CHECK_EQUAL(image12.get_local_ID(), expected_i12_local_ID);
	BOOST_CHECK_EQUAL(image21.get_local_ID(), expected_i21_local_ID);

	// Check that everything has the proper ID sequence
	auto survey_ID_seq = survey.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(survey_ID_seq.begin(), survey_ID_seq.end(), expected_survey_ID_seq.begin(), expected_survey_ID_seq.end());
	auto ig1_ID_seq = image_group1.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(ig1_ID_seq.begin(), ig1_ID_seq.end(), expected_ig1_ID_seq.begin(), expected_ig1_ID_seq.end());
	auto ig2_ID_seq = image_group2.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(ig2_ID_seq.begin(), ig2_ID_seq.end(), expected_ig2_ID_seq.begin(), expected_ig2_ID_seq.end());
	auto ig3_ID_seq = image_group3.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(ig3_ID_seq.begin(), ig3_ID_seq.end(), expected_ig3_ID_seq.begin(), expected_ig3_ID_seq.end());
	auto i11_ID_seq = image11.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(i11_ID_seq.begin(), i11_ID_seq.end(), expected_i11_ID_seq.begin(), expected_i11_ID_seq.end());
	auto i12_ID_seq = image12.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(i12_ID_seq.begin(), i12_ID_seq.end(), expected_i12_ID_seq.begin(), expected_i12_ID_seq.end());
	auto i21_ID_seq = image21.get_ID_seq();
	BOOST_CHECK_EQUAL_COLLECTIONS(i21_ID_seq.begin(), i21_ID_seq.end(), expected_i21_ID_seq.begin(), expected_i21_ID_seq.end());

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
