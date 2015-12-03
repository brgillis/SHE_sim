/**********************************************************************\
 @file object_types.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_OBJECT_TYPES_HPP_
#define SHE_SIM_GAL_PARAMS_OBJECT_TYPES_HPP_

#include <memory>

#include <GalSim.h>

#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

// Typedefs for the type of each object

typedef std::shared_ptr<galsim::SBProfile> p_profile_t;

typedef p_profile_t background_psf_t;

typedef array_t<p_profile_t> binned_observed_flux_distribution_t;
typedef array_t<p_profile_t> binned_psf_t;
typedef p_profile_t core_observed_flux_distribution_t;
typedef array_1d_t core_sed_t;
typedef p_profile_t disk_observed_flux_distribution_t;
typedef array_1d_t disk_sed_t;
typedef p_profile_t observed_flux_distribution_t;
typedef array_2d_t psf_model_t;
typedef array_2d_t pix_galaxy_w_pois_noise_t;

}




#endif // SHE_SIM_GAL_PARAMS_OBJECT_TYPES_HPP_
