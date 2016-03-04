/**********************************************************************\
 @file PHL_seed_test.cpp
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
#include "SHE_SIM_gal_params/param_params/IndGaussian.hpp"

namespace SHE_SIM
{

struct PHL_seed_fixture {

	Survey survey1;

	int_t test_seed1 = 15;

};


BOOST_AUTO_TEST_SUITE (PHL_seed_Test)

BOOST_FIXTURE_TEST_CASE(test_PHL_ID, PHL_seed_fixture) {

	// Setup
	survey1.set_param_params<IndGaussian>(exp_time_name,565.,10.);
	survey1.set_generation_level(exp_time_name,dv::image_group_level);

	Survey survey2 = survey1;

	survey1.set_seed(test_seed1);

	ImageGroup image_group11 = *survey1.add_image_group();
	ImageGroup image_group12 = *survey1.add_image_group();

	ImageGroup image_group21 = *survey2.add_image_group();
	ImageGroup image_group22 = *survey2.add_image_group();

	survey2.set_seed(test_seed1);

	// Check it's as expected
	BOOST_CHECK_NE(survey1.get_param_value(exp_time_name),image_group11.get_param_value(exp_time_name));
	BOOST_CHECK_NE(image_group11.get_param_value(exp_time_name),image_group12.get_param_value(exp_time_name));
	BOOST_CHECK_CLOSE(image_group11.get_param_value(exp_time_name),image_group21.get_param_value(exp_time_name),1e-9);
	BOOST_CHECK_CLOSE(image_group12.get_param_value(exp_time_name),image_group22.get_param_value(exp_time_name),1e-9);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
