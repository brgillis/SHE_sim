/**********************************************************************\
 @file ParamParam_test.cpp
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
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/param_params/IndFixed.hpp"
#include "SHE_SIM_gal_params/param_params/IndContRayleigh.hpp"
#include "SHE_SIM_gal_params/param_params/IndGaussian.hpp"
#include "SHE_SIM_gal_params/param_params/IndLogNormalMean.hpp"
#include "SHE_SIM_gal_params/param_params/IndLogNormalPeak.hpp"
#include "SHE_SIM_gal_params/param_params/IndPoisson.hpp"
#include "SHE_SIM_gal_params/param_params/IndRayleigh.hpp"
#include "SHE_SIM_gal_params/param_params/IndTruncGaussian.hpp"
#include "SHE_SIM_gal_params/param_params/IndTruncRayleigh.hpp"
#include "SHE_SIM_gal_params/param_params/IndUniform.hpp"

namespace SHE_SIM
{

struct param_param_fixture {

	const flt_t fixed_val = 3412.1;

	const flt_t uni_min = 21.;
	const flt_t uni_max = 124.;

	const flt_t mean = 101.;
	const flt_t stddev = 0.01;

	const flt_t min = 101.010;
	const flt_t max = 101.015;

	const flt_t l10_mean = -7;
	const flt_t l10_peak = 17;
	const flt_t l10_stddev = 0.001;

	const flt_t r_sigma = 27.;
	const flt_t r_max = 10;
	const flt_t r_p = 4;

	const flt_t lambda = 100000;

	const flt_t acceptable_sigma = 5.;

	const IndFixed fixed_param = IndFixed(fixed_val);
	const IndUniform uni_param = IndUniform(uni_min,uni_max);

	const IndGaussian gaus_param = IndGaussian(mean,stddev);
	const IndTruncGaussian trunc_gaus_param = IndTruncGaussian(mean,stddev,min,max);
	const IndLogNormalMean lognorm_mean_param = IndLogNormalMean(l10_mean,l10_stddev);
	const IndLogNormalPeak lognorm_peak_param = IndLogNormalPeak(l10_peak,l10_stddev);

	const IndRayleigh ray_param = IndRayleigh(r_sigma);
	const IndTruncRayleigh trunc_ray_param = IndTruncRayleigh(r_sigma,r_max);
	const IndContRayleigh cont_ray_param = IndContRayleigh(r_sigma,r_max,r_p);

	const IndPoisson pois_param = IndPoisson(lambda);

	const seed_t seed_1 = 0;
	const seed_t seed_2 = 15;
};


BOOST_AUTO_TEST_SUITE (Param_Param_Test)

BOOST_FIXTURE_TEST_CASE(test_param_param, param_param_fixture) {

	// Test the fixed value
	const ParamParam * p_param_param = &fixed_param;

	BOOST_CHECK_EQUAL( p_param_param->get_mode(), ParamParam::INDEPENDENT);

	BOOST_CHECK_EQUAL( p_param_param->get_independently(), fixed_val);

	// Test the uniform distribution
	p_param_param = &uni_param;

	flt_t uni_val_1 = p_param_param->get_independently();
	flt_t uni_val_2 = p_param_param->get_independently();

	BOOST_CHECK_NE( uni_val_1, uni_val_2  );

	BOOST_CHECK_GE( uni_val_1, uni_min );
	BOOST_CHECK_LE( uni_val_1, uni_max );
	BOOST_CHECK_GE( uni_val_2, uni_min );
	BOOST_CHECK_LE( uni_val_2, uni_max );

	// Test that we can pass our own generator
	gen_t gen;

	gen.seed(seed_1);
	flt_t seeded_uni_val_1 = p_param_param->get_independently(gen);
	gen.seed(seed_2);
	flt_t seeded_uni_val_2 = p_param_param->get_independently(gen);
	gen.seed(seed_1);
	flt_t seeded_uni_val_3 = p_param_param->get_independently(gen);

	BOOST_CHECK_NE( seeded_uni_val_1, seeded_uni_val_2  );
	BOOST_CHECK_EQUAL( seeded_uni_val_1, seeded_uni_val_3  );

	BOOST_CHECK_GE( seeded_uni_val_1, uni_min );
	BOOST_CHECK_LE( seeded_uni_val_1, uni_max );
	BOOST_CHECK_GE( seeded_uni_val_2, uni_min );
	BOOST_CHECK_LE( seeded_uni_val_2, uni_max );

}

BOOST_FIXTURE_TEST_CASE(test_gaussian, param_param_fixture) {

	flt_t val = gaus_param.get_independently();

	BOOST_CHECK_GE( val, mean-acceptable_sigma*stddev );
	BOOST_CHECK_LE( val, mean+acceptable_sigma*stddev );

}

BOOST_FIXTURE_TEST_CASE(test_trunc_gaussian, param_param_fixture) {

	flt_t val = trunc_gaus_param.get_independently();

	BOOST_CHECK_GE( val, min );
	BOOST_CHECK_LE( val, max );

}

BOOST_FIXTURE_TEST_CASE(test_lognormal_mean, param_param_fixture) {

	flt_t val = lognorm_mean_param.get_independently();

	BOOST_CHECK_GE( val, std::pow(10.,l10_mean-acceptable_sigma*l10_stddev) );
	BOOST_CHECK_LE( val, std::pow(10.,l10_mean+acceptable_sigma*l10_stddev) );

}

BOOST_FIXTURE_TEST_CASE(test_lognormal_peak, param_param_fixture) {

	flt_t val = lognorm_peak_param.get_independently();

	BOOST_CHECK_GE( val, std::pow(10.,l10_peak-acceptable_sigma*l10_stddev) );
	BOOST_CHECK_LE( val, std::pow(10.,l10_peak+acceptable_sigma*l10_stddev) );

}

BOOST_FIXTURE_TEST_CASE(test_rayleigh, param_param_fixture) {

	flt_t val = ray_param.get_independently();

	BOOST_CHECK_LE( val, acceptable_sigma*r_sigma );

}

BOOST_FIXTURE_TEST_CASE(test_trunc_rayleigh, param_param_fixture) {

	flt_t val = trunc_ray_param.get_independently();

	BOOST_CHECK_LE( val, r_max );

}

BOOST_FIXTURE_TEST_CASE(test_cont_rayleigh, param_param_fixture) {

	flt_t val = cont_ray_param.get_independently();

	BOOST_CHECK_LE( val, r_max );

}

BOOST_FIXTURE_TEST_CASE(test_poisson, param_param_fixture) {

	flt_t val = pois_param.get_independently();

	BOOST_CHECK_GE( val, lambda-acceptable_sigma*std::sqrt(lambda) );
	BOOST_CHECK_LE( val, lambda+acceptable_sigma*std::sqrt(lambda) );

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
