/**********************************************************************\
 @file MagVis_test.cpp
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
#include "SHE_SIM_gal_params/levels/Survey.hpp"

namespace SHE_SIM
{

struct mag_vis_fixture {

	Survey survey;

	const flt_t exp_time1 = 1234.5;
	const flt_t exp_time2 = 2468.0;

	const flt_t mag_vis_inst_zp = 25.0;

	const flt_t expected_mag_vis_zp1 = 25 + 2.5 * std::log10(exp_time1);
	const flt_t expected_mag_vis_zp2 = 25 + 2.5 * std::log10(exp_time2);

};


BOOST_AUTO_TEST_SUITE (Mag_Vis_Test)

BOOST_FIXTURE_TEST_CASE(test_mag_vis, mag_vis_fixture) {

	survey.set_generation_level("exp_time",0);
	survey.set_generation_level("mag_vis_inst_zp",0);
	survey.set_generation_level("mag_vis_zp",0);

	survey.set_param_params("exp_time",std::vector<flt_t>({exp_time1}));
	survey.set_param_params("mag_vis_inst_zp",std::vector<flt_t>({mag_vis_inst_zp}));

	BOOST_CHECK_EQUAL(survey.get_param_value("mag_vis_zp"),expected_mag_vis_zp1);

	survey.set_param_params("exp_time",std::vector<flt_t>({exp_time2}));

	BOOST_CHECK_EQUAL(survey.get_param_value("mag_vis_zp"),expected_mag_vis_zp2);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
