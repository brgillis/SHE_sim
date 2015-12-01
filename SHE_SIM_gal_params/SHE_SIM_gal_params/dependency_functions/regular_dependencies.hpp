/**********************************************************************\
 @file regular_dependencies.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_REGULAR_DEPENDENCIES_HPP_
#define SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_REGULAR_DEPENDENCIES_HPP_

#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM {

// Inlined functions first

inline flt_t get_image_area( const flt_t & image_size_x_pix, const flt_t & image_size_y_pix, const flt_t & pixel_scale )
{
	return image_size_x_pix * image_size_y_pix * square( pixel_scale );
}

inline flt_t get_zp( const flt_t & inst_zp, const flt_t & exp_time)
{
	return inst_zp + 2.5 * std::log10(exp_time);
}

inline flt_t generate_count( const flt_t & lambda, gen_t & rng )
{
	return Pois_rand( lambda , rng);
}

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_REGULAR_DEPENDENCIES_HPP_
