/**********************************************************************\
 @file cdf_rand_test.cpp
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
#include "IceBRG_main/math/random/random_functions.hpp"

namespace SHE_SIM
{

struct cdf_rand_fixture {

	typedef Eigen::Array<flt_t,Eigen::Dynamic,1> flt_array_t;

	const seed_t seed_1 = 1;
	const seed_t seed_2 = 5;

	flt_array_t xvals = flt_array_t(6);

	flt_array_t cdf_low = flt_array_t(6);
	flt_array_t cdf_high = flt_array_t(6);

	flt_array_t pdf_low = flt_array_t(6);
	flt_array_t pdf_high = flt_array_t(6);
};


BOOST_AUTO_TEST_SUITE (CDF_Rand_Test)

BOOST_FIXTURE_TEST_CASE(test_cdf_rand, cdf_rand_fixture) {

	xvals << 0, 1, 2, 3, 4, 5;

	cdf_low << 0, 1, 1, 1, 1, 1;
	cdf_high << 0, 0, 0, 0, 10, 10;

	pdf_low << 0, 1, 0, 0, 0, 0;
	pdf_high << 0, 0, 0, 0, 10, 0;

	std::function<flt_t(flt_t const &, flt_array_t const &)> const & get_from_array = [&] ( flt_t const & x, flt_array_t const & a)
	{
		int i = x;
		if(i >= a.size()) i = a.size()-1;
		if(i < 0) i = 0;
		return a[i];
	};

	std::function<flt_t(flt_t const &)> const & get_from_cdf_low_array = [&] ( flt_t const & x ) { return get_from_array(x,cdf_low); };
	std::function<flt_t(flt_t const &)> const & get_from_cdf_high_array = [&] ( flt_t const & x ) { return get_from_array(x,cdf_high); };

	std::function<flt_t(flt_t const &)> const & get_from_pdf_low_array = [&] ( flt_t const & x ) { return get_from_array(x,pdf_low); };
	std::function<flt_t(flt_t const &)> const & get_from_pdf_high_array = [&] ( flt_t const & x ) { return get_from_array(x,pdf_high); };

	flt_t r_cdf_low = IceBRG::rand_from_cdf(get_from_cdf_low_array,6.,0.,5.);
	BOOST_CHECK_GE(r_cdf_low,0.);
	BOOST_CHECK_LE(r_cdf_low,1.);

	flt_t r_cdf_high = IceBRG::rand_from_cdf(get_from_cdf_high_array,6.,0.,5.);
	BOOST_CHECK_GE(r_cdf_high,3.);
	BOOST_CHECK_LE(r_cdf_high,4.);

	flt_t r_pdf_low = IceBRG::rand_from_pdf(get_from_pdf_low_array,6.,0.,5.);
	BOOST_CHECK_GE(r_pdf_low,0.);
	BOOST_CHECK_LE(r_pdf_low,1.);

	flt_t r_pdf_high = IceBRG::rand_from_pdf(get_from_pdf_high_array,6.,0.,5.);
	BOOST_CHECK_GE(r_pdf_high,3.);
	BOOST_CHECK_LE(r_pdf_high,4.);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
