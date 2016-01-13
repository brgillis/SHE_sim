/**********************************************************************\
 @file cosmology.cpp
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

#include <math.h>

#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/dependency_functions/cosmology.hpp"
#include "SHE_SIM_gal_params/dependency_functions/dfa_cache.hpp"
#include "SHE_SIM_gal_params/math.hpp"

namespace SHE_SIM {

flt_t H( flt_t const & z )
{
	if(z==0) return H_0;

	// Friedmann equation, assuming omega = -1

	flt_t zp1 = 1.+z;

	return H_0 * std::sqrt( Omega_r * square( zp1 )*square( zp1 )
							+ Omega_m * zp1 * square( zp1 )
							+ Omega_k * square( zp1 ) + Omega_l );
}

flt_t get_dfa( const flt_t & z )
{
	dfa_array_t z_array = dfa_cache.first;
	dfa_array_t dfa_array = dfa_cache.second;

	return interpolate(z,z_array,dfa_array);
}

flt_t get_distance_from_angle( const flt_t & theta_arcsec, const flt_t & z )
{
	flt_t res = theta_arcsec * get_dfa(z);

	return res;
}

flt_t get_angle_from_distance( const flt_t & d_kpc, const flt_t & z )
{
	flt_t res = d_kpc / get_dfa(z);

	return res;
}

} // namespace SHE_SIM


