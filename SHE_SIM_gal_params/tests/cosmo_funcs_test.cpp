/**********************************************************************\
 @file cosmo_funcs_test.cpp
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

#include <Eigen/Core>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/dependency_functions/cosmology.hpp"

namespace SHE_SIM
{

struct cosmo_funcs_fixture {

	const flt_t tolerance = 5; // A bit loose, to allow different cosmologies

	const flt_t test_theta_0 = 0.;
	const flt_t test_theta_1 = 0.2;
	const flt_t test_theta_2 = 0.5;

	const flt_t test_z_0 = 0.;
	const flt_t test_z_1 = 0.1;
	const flt_t test_z_2 = 1.0;

	const flt_t exp_d_00 = 0.;
	const flt_t exp_d_01 = 0.;
	const flt_t exp_d_02 = 0.;
	const flt_t exp_d_10 = 0.;
	const flt_t exp_d_11 = 0.3692;
	const flt_t exp_d_12 = 1.6132;
	const flt_t exp_d_20 = 0.;
	const flt_t exp_d_21 = 0.923;
	const flt_t exp_d_22 = 4.033;

};


BOOST_AUTO_TEST_SUITE (Cosmo_Funcs_Test)

BOOST_FIXTURE_TEST_CASE(test_cosmo_funcs, cosmo_funcs_fixture) {

	flt_t d_00 = SHE_SIM::get_distance_from_angle(test_theta_0,test_z_0);
	BOOST_CHECK_EQUAL(d_00,exp_d_00);
	flt_t d_01 = SHE_SIM::get_distance_from_angle(test_theta_0,test_z_1);
	BOOST_CHECK_EQUAL(d_01,exp_d_01);
	flt_t d_02 = SHE_SIM::get_distance_from_angle(test_theta_0,test_z_2);
	BOOST_CHECK_EQUAL(d_02,exp_d_02);

	flt_t d_10 = SHE_SIM::get_distance_from_angle(test_theta_1,test_z_0);
	BOOST_CHECK_EQUAL(d_10,exp_d_10);
	flt_t d_11 = SHE_SIM::get_distance_from_angle(test_theta_1,test_z_1);
	BOOST_CHECK_CLOSE(d_11,exp_d_11,tolerance);
	flt_t d_12 = SHE_SIM::get_distance_from_angle(test_theta_1,test_z_2);
	BOOST_CHECK_CLOSE(d_12,exp_d_12,tolerance);

	flt_t d_20 = SHE_SIM::get_distance_from_angle(test_theta_2,test_z_0);
	BOOST_CHECK_EQUAL(d_20,exp_d_20);
	flt_t d_21 = SHE_SIM::get_distance_from_angle(test_theta_2,test_z_1);
	BOOST_CHECK_CLOSE(d_21,exp_d_21,tolerance);
	flt_t d_22 = SHE_SIM::get_distance_from_angle(test_theta_2,test_z_2);
	BOOST_CHECK_CLOSE(d_22,exp_d_22,tolerance);

	// Check inverse

	flt_t a_11 = SHE_SIM::get_angle_from_distance(exp_d_11,test_z_1);
	BOOST_CHECK_CLOSE(a_11,test_theta_1,tolerance);
	flt_t a_12 = SHE_SIM::get_angle_from_distance(exp_d_12,test_z_2);
	BOOST_CHECK_CLOSE(a_12,test_theta_1,tolerance);
	flt_t a_21 = SHE_SIM::get_angle_from_distance(exp_d_21,test_z_1);
	BOOST_CHECK_CLOSE(a_21,test_theta_2,tolerance);
	flt_t a_22 = SHE_SIM::get_angle_from_distance(exp_d_22,test_z_2);
	BOOST_CHECK_CLOSE(a_22,test_theta_2,tolerance);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
