/**********************************************************************\
 @file PhysicalSize_test.cpp
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
#include <SHE_SIM_gal_params/default_values.hpp>
#include <SHE_SIM_gal_params/math.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include "SHE_SIM_gal_params/dependency_functions/galaxy_size_detail.hpp"
#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"

namespace SHE_SIM
{

struct physical_size_fixture {

	Survey survey;

	int_t const num_gals = 1000;

	flt_t const z_1 = 0.25;
	flt_t const z_2 = 1.25;

	flt_t const l10_ms_1 = 9.75;
	flt_t const l10_ms_2 = 11.25;

	flt_t const ex_mean_l10_r_bulge_11 = std::log10(1.9);
	flt_t const ex_mean_l10_r_bulge_12 = std::log10(8.7);
	flt_t const ex_mean_l10_r_bulge_21 = std::log10(1.8);
	flt_t const ex_mean_l10_r_bulge_22 = std::log10(4.0);

	flt_t const ex_scatter_l10_r_bulge_11 = 0.10;
	flt_t const ex_scatter_l10_r_bulge_12 = 0.12;
	flt_t const ex_scatter_l10_r_bulge_21 = 0.10;
	flt_t const ex_scatter_l10_r_bulge_22 = 0.12;

	flt_t const ex_mean_l10_r_disk_11 = std::log10(4.1);
	flt_t const ex_mean_l10_r_disk_12 = std::log10(9.8);
	flt_t const ex_mean_l10_r_disk_21 = std::log10(3.2);
	flt_t const ex_mean_l10_r_disk_22 = std::log10(6.3);

	flt_t const ex_scatter_l10_r_disk_11 = 0.16;
	flt_t const ex_scatter_l10_r_disk_12 = 0.17;
	flt_t const ex_scatter_l10_r_disk_21 = 0.16;
	flt_t const ex_scatter_l10_r_disk_22 = 0.17;

	flt_t const tol_mean = 10; // Since I'm guessing from a graph here
	flt_t const tol_scatter = 4;

};


BOOST_AUTO_TEST_SUITE (PhysicalSize_Test)

BOOST_FIXTURE_TEST_CASE(test_physical_size, physical_size_fixture) {

	// Now see if we can inherit all the way down to the galaxy level
	survey.set_param_params(num_clusters_name,"fixed",0.);
	survey.set_param_params(num_field_galaxies_name,"fixed",flt_t(num_gals));

	survey.set_generation_level(redshift_name,dv::survey_level);
	survey.set_generation_level(stellar_mass_name,dv::survey_level);

	survey.autofill_children();

	auto galaxies = survey.get_galaxy_descendants();

	gal_size_array_t l10_r_bulge_array = gal_size_array_t::Zero(num_gals);
	gal_size_array_t l10_r_disk_array = gal_size_array_t::Zero(num_gals);

	// Test for each set of redshift and stellar mass

	// 1-1

	survey.set_param_params(redshift_name,"fixed",z_1);
	survey.set_param_params(stellar_mass_name,"fixed",std::pow(10.,l10_ms_1));

	for( int_t i=0; i<num_gals; ++i )
	{
		l10_r_bulge_array[i] = std::log10(galaxies.at(i)->get_param_value(physical_size_bulge_name));
		l10_r_disk_array[i] = std::log10(galaxies.at(i)->get_param_value(physical_size_disk_name));
	}

	flt_t mean_r_bulge = l10_r_bulge_array.mean();
	flt_t mean_r_disk = l10_r_disk_array.mean();

	flt_t scatter_r_bulge = std::sqrt(square(l10_r_bulge_array).mean() - square(l10_r_bulge_array.mean()));
	flt_t scatter_r_disk = std::sqrt(square(l10_r_disk_array).mean() - square(l10_r_disk_array.mean()));

	BOOST_CHECK_CLOSE(mean_r_bulge,ex_mean_l10_r_bulge_11,tol_mean);
	BOOST_CHECK_CLOSE(mean_r_disk,ex_mean_l10_r_disk_11,tol_mean);

	BOOST_CHECK_CLOSE(scatter_r_bulge,ex_scatter_l10_r_bulge_11,tol_scatter);
	BOOST_CHECK_CLOSE(scatter_r_disk,ex_scatter_l10_r_disk_11,tol_scatter);

}

BOOST_AUTO_TEST_SUITE_END ()

} // namespace SHE_SIM
