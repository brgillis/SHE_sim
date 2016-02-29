/**********************************************************************\
 @file galaxy_physical_size.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2016 brg

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

#include "IceBRG_main/math/random/random_functions.hpp"
#include "IceBRG_physics/constants.hpp"
#include "IceBRG_physics/cosmology.hpp"

#include "SHE_SIM_gal_params/common.hpp"

#include "SHE_SIM_gal_params/dependency_functions/cosmology.hpp"
#include "SHE_SIM_gal_params/dependency_functions/galaxy_size_detail.hpp"
#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/math.hpp"

namespace SHE_SIM {

using namespace IceBRG;

// Implementation of functions

flt_t generate_physical_size( flt_t const & redshift,
		flt_t const & stellar_mass,
		med_gal_size_cache_t const & med_gal_size_cache,
		gal_size_scatter_cache_t const & gal_size_scatter_cache,
		gen_t & rng )
{
	// Simplify labels for arrays contained in the caches

	gal_size_array_t const & l10_ms_array = std::get<0>(med_gal_size_cache);
	gal_size_array_t const & l10_B_array = std::get<1>(med_gal_size_cache);
	gal_size_array_t const & beta_array = std::get<2>(med_gal_size_cache);

	gal_size_array_t const & z_array = std::get<0>(gal_size_scatter_cache);
	gal_size_array_t const & sigma_l10_R_array = std::get<1>(gal_size_scatter_cache);

	// Get median from interpolation of functional form

	flt_t l10_ms = std::log10(stellar_mass);

	flt_t l10_B = interpolate(l10_ms,l10_ms_array,l10_B_array);
	flt_t beta = interpolate(l10_ms,l10_ms_array,beta_array);

	flt_t l10_med_gal_size = l10_B + beta * std::log10(H(redshift)/H_0);

	// Get scatter, then generate a random value

	flt_t sigma_l10_gal_size = interpolate(redshift,z_array,sigma_l10_R_array);

	flt_t l10_gal_size = IceBRG::Gaus_rand(l10_med_gal_size,sigma_l10_gal_size,rng);

	return std::pow(10.,l10_gal_size);

}

flt_t generate_physical_size_bulge( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, gen_t & rng  )
{
	return generate_physical_size( redshift,
			stellar_mass,
			med_gal_size_bulge_cache,
			gal_size_scatter_bulge_cache,
			rng );
}

flt_t generate_physical_size_disk( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, gen_t & rng  )
{
	return generate_physical_size( redshift,
			stellar_mass,
			med_gal_size_disk_cache,
			gal_size_scatter_disk_cache,
			rng );
}

// Implementation of global constant caches

const med_gal_size_cache_t med_gal_size_bulge_cache = load_med_gal_size_bulge_cache();
const med_gal_size_cache_t med_gal_size_disk_cache = load_med_gal_size_disk_cache();
const gal_size_scatter_cache_t gal_size_scatter_bulge_cache = load_gal_size_scatter_bulge_cache();
const gal_size_scatter_cache_t gal_size_scatter_disk_cache = load_gal_size_scatter_disk_cache();

// Implementation of loading functions for caches

med_gal_size_cache_t load_med_gal_size_bulge_cache()
{
	gal_size_array_t l10_ms_array(4), l10_B_array(4), beta_array(4);

	l10_ms_array <<
		9.75,
		10.25,
		10.75,
		11.25;

	l10_B_array <<
		0.27,
		0.42,
		0.68,
		0.97;

	beta_array <<
		-0.19,
		-0.97,
		-1.13,
		-1.29;

	return med_gal_size_cache_t(std::move(l10_ms_array),
			std::move(l10_B_array),
			std::move(beta_array));
}

med_gal_size_cache_t load_med_gal_size_disk_cache()
{
	gal_size_array_t l10_ms_array(5), l10_B_array(5), beta_array(5);

	l10_ms_array <<
		9.25,
		9.75,
		10.25,
		10.75,
		11.25;

	l10_B_array <<
		0.52,
		0.65,
		0.71,
		0.86,
		1.01;

	beta_array <<
		-0.52,
		-0.58,
		-0.49,
		-0.65,
		-0.76;

	return med_gal_size_cache_t(std::move(l10_ms_array),
			std::move(l10_B_array),
			std::move(beta_array));
}

gal_size_scatter_cache_t load_gal_size_scatter_bulge_cache()
{
	gal_size_array_t z_array(6), sigma_l10_R_array(6);

	z_array <<
		0.25,
		0.75,
		1.25,
		1.75,
		2.25,
		2.75;

	sigma_l10_R_array <<
		0.10,
		0.11,
		0.12,
		0.14,
		0.14,
		0.14;

	return gal_size_scatter_cache_t(std::move(z_array),
			std::move(sigma_l10_R_array));
}

gal_size_scatter_cache_t load_gal_size_scatter_disk_cache()
{
	gal_size_array_t z_array(6), sigma_l10_R_array(6);

	z_array <<
		0.25,
		0.75,
		1.25,
		1.75,
		2.25,
		2.75;

	sigma_l10_R_array <<
		0.16,
		0.16,
		0.17,
		0.18,
		0.19,
		0.19;

	return gal_size_scatter_cache_t(std::move(z_array),
			std::move(sigma_l10_R_array));
}

} // namespace SHE_SIM


