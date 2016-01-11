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

#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/dependency_functions/cosmology.hpp"
#include "SHE_SIM_gal_params/dependency_functions/dfa_cache.hpp"

namespace SHE_SIM {

flt_t get_distance_from_angle( const flt_t & theta_arcsec, const flt_t & z )
{
	dfa_array_t z_array = dfa_cache.first;
	dfa_array_t dfa_array = dfa_cache.second;

	dfa_array_t diffs = (z_array-z).abs();
	int_t z_i,z_j;
	diffs.minCoeff(&z_i,&z_j);

	// We want the lower index around this redshift if possible
	if(z_array[z_i] > z) --z_i;
	if(z_i<0) z_i = 0;

	flt_t z_lo = z_array[z_i];
	flt_t z_hi = z_array[z_i+1];
	flt_t dfa_lo = dfa_array[z_i];
	flt_t dfa_hi = dfa_array[z_i+1];

	flt_t dfa = dfa_lo + (dfa_hi-dfa_lo)/(z_hi-z_lo) * (z-z_lo);

	flt_t res = theta_arcsec * dfa;

	return res;
}

} // namespace SHE_SIM


