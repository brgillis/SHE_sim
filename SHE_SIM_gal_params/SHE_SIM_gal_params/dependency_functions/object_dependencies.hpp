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
#include "SHE_SIM_gal_params/object_types.hpp"
#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM {

// Image-level

background_psf_t get_background_psf( flt_t const & psf_params );

// Galaxy-level

binned_observed_flux_distribution_t get_binned_observed_flux_distribution(
		core_sed_t const & (core_sed),
		disk_sed_t const & (disk_sed),
		core_observed_flux_distribution_t const & core_observed_flux_distribution,
		disk_observed_flux_distribution_t const & disk_observed_flux_distribution);

binned_psf_t get_binned_psf(psf_model_t const & psf_model, flt_t const & xp, flt_t const & yp);

core_sed_t get_core_sed( flt_t const & morphology, flt_t const & redshift, flt_t const & stellar_mass );

core_observed_flux_distribution_t get_core_observed_flux_distribution(
		flt_t const & morphology, flt_t const & rotation, flt_t const & tilt);

disk_sed_t get_disk_sed( flt_t const & morphology, flt_t const & redshift, flt_t const & stellar_mass );

disk_observed_flux_distribution_t get_disk_observed_flux_distribution(
		flt_t const & morphology, flt_t const & rotation, flt_t const & tilt);

observed_flux_distribution_t get_observed_flux_distribution(
		binned_observed_flux_distribution_t const & binned_observed_flux_distribution,
		binned_psf_t const & binned_psf);

psf_model_t get_psf_model( flt_t const & psf_params );

pix_galaxy_w_pois_noise_t get_pix_galaxy_w_pois_noise(
		observed_flux_distribution_t observed_flux_distribution,
		flt_t const & xp, flt_t const & yp, flt_t const & pixel_scale, flt_t const & gain);

} // namespace SHE_SIM



#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_OBJECT_DEPENDENCIES_HPP_
