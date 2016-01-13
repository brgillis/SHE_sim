/**********************************************************************\
 @file math.hpp
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

#ifndef MATH_HPP_
#define MATH_HPP_

#include <cassert>

#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

template< typename T >
decltype(T()*T()) square( T const & v1)
{
	return v1*v1;
}

template< typename T >
int_t round_int( T const & v )
{
	return static_cast<int_t>(v+0.5);
}

template< typename T_array >
flt_t interpolate( flt_t const & x, T_array const & x_array, T_array const & y_array )
{
	assert(x_array.size()==y_array.size());

	T_array diffs = (x_array-x).abs();
	int_t x_i,x_j;
	diffs.minCoeff(&x_i,&x_j);

	// We want the lower index around this redshift if possible
	if(x_array[x_i] > x) --x_i;
	if(x_i<0) x_i = 0;
	if(x_i>=x_array.size()-1) x_i = x_array.size()-2;

	flt_t x_lo = x_array[x_i];
	flt_t x_hi = x_array[x_i+1];
	flt_t y_lo = y_array[x_i];
	flt_t y_hi = y_array[x_i+1];

	flt_t y = y_lo + (y_hi-y_lo)/(x_hi-x_lo) * (x-x_lo);

	return y;
}

} // namespace SHE_SIM

#endif // MATH_HPP_
