/**********************************************************************\
 @file ExposureTime_test.cpp
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
#include <SHE_SIM_gal_params/param_declarations.hpp>
#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"

namespace SHE_SIM
{

struct exp_time_fixture {

	Survey survey;

	const IndFixed exp_time1_param = IndFixed(1234.5);

	const flt_t exp_time2 = 2468.0;

	const flt_t exp_time3 = 34151.1;

};


BOOST_AUTO_TEST_SUITE (Exp_Time_Test)

BOOST_FIXTURE_TEST_CASE(test_exp_time, exp_time_fixture) {

	survey.set_generation_level(exp_time_name,0);

	survey.set_p_param_params(exp_time_name,&exp_time1_param);

	BOOST_CHECK_EQUAL(survey.get_param_value(exp_time_name),exp_time1_param.get_independently());

	survey.set_param_params<IndFixed>(exp_time_name,exp_time2);

	BOOST_CHECK_NE(survey.get_param_value(exp_time_name),exp_time1_param.get_independently());

	BOOST_CHECK_EQUAL(survey.get_param_value(exp_time_name),exp_time2);

	survey.set_param_params(exp_time_name,"fixed",exp_time3);

	BOOST_CHECK_NE(survey.get_param_value(exp_time_name),exp_time2);

	BOOST_CHECK_EQUAL(survey.get_param_value(exp_time_name),exp_time3);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
