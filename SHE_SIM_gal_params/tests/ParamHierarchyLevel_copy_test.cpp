/**********************************************************************\
 @file ParamHierarchyLevel_copy_test.cpp
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

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/levels/Survey.hpp"

namespace SHE_SIM
{

struct PHL_copy_fixture {

	Survey survey1;

	const IndFixed exp_time1 = IndFixed(1234.5);
	const IndFixed exp_time2 = IndFixed(2468.0);

};


BOOST_AUTO_TEST_SUITE (PHL_Copy_Test)

BOOST_FIXTURE_TEST_CASE(test_PHL_copy_construct, PHL_copy_fixture) {

	survey1.set_generation_level(exp_time_name,0);

	survey1.set_p_param_params(exp_time_name,&exp_time1);

	BOOST_CHECK_EQUAL(survey1.get_param_value(exp_time_name),exp_time1.get_independently());

	Survey survey2(survey1);

	survey1.set_p_param_params(exp_time_name,&exp_time2);

	BOOST_CHECK_EQUAL(survey1.get_param_value(exp_time_name),exp_time2.get_independently());

	BOOST_CHECK_EQUAL(survey2.get_param_value(exp_time_name),exp_time1.get_independently());

}

BOOST_FIXTURE_TEST_CASE(test_PHL_copy_assign, PHL_copy_fixture) {

	survey1.set_generation_level(exp_time_name,0);

	survey1.set_p_param_params(exp_time_name,&exp_time1);

	Survey survey2;

	survey2 = survey1;

	survey1.set_p_param_params(exp_time_name,&exp_time2);

	BOOST_CHECK_EQUAL(survey1.get_param_value(exp_time_name),exp_time2.get_independently());

	BOOST_CHECK_EQUAL(survey2.get_param_value(exp_time_name),exp_time1.get_independently());

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
