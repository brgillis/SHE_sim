/**********************************************************************\
 @file cosmology.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_COSMOLOGY_HPP_
#define SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_COSMOLOGY_HPP_

#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

/**
 * Get a transverse distance in kpc from an angle in arcsec
 *
 * @param theta_arcsec angle in arcsec
 * @param z
 * @return distance in kpc
 */
flt_t get_distance_from_angle( const flt_t & theta_arcsec, const flt_t & z );

/**
 * Get an angle in kpc from a transverse distance in kpc
 *
 * @param d_kpc transverse distance in kpc
 * @param z
 * @return angle in arcsec
 */
flt_t get_angle_from_distance( const flt_t & d_kpc, const flt_t & z );

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_COSMOLOGY_HPP_
