/**********************************************************************\
 @file GenerateParameters_test.cpp
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
#include "SHE_SIM_gal_params/levels/Cluster.hpp"
#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"

namespace SHE_SIM
{

struct gen_params_fixture {

	Survey survey;

};


BOOST_AUTO_TEST_SUITE (Generate_Parameters_Test)

BOOST_FIXTURE_TEST_CASE(test_gen_params, gen_params_fixture) {

	BOOST_CHECK_NO_THROW(survey.generate_parameters());

	Cluster cluster;

	cluster.add_central_galaxy();

	BOOST_CHECK_NO_THROW(survey.autofill_children());

	BOOST_CHECK_NO_THROW(survey.generate_parameters());

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
