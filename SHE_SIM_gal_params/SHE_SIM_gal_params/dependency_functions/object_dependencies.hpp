/**********************************************************************\
 @file object_dependencies.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_OBJECT_DEPENDENCIES_HPP_
#define SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_OBJECT_DEPENDENCIES_HPP_

#include "SHE_SIM_gal_params/common.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM {

// Image-level

array_2d_t get_background_psf( flt_t const & psf_params);

// Galaxy-level

array_1d_t get_psf_model( flt_t const & psf_params );

array_1d_t get_sed( flt_t const & morphology, flt_t const & redshift, flt_t const & stellar_mass );

} // namespace SHE_SIM



#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_OBJECT_DEPENDENCIES_HPP_
