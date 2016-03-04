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

#ifndef SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_MISC_DEPENDENCIES_HPP_
#define SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_MISC_DEPENDENCIES_HPP_

#include "IceBRG_main/math/misc_math.hpp"
#include "IceBRG_main/math/random/random_functions.hpp"

#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

// Inlined functions first

/**
 * Get the image area in square arcminutes.
 *
 * @param image_size_x_pix
 * @param image_size_y_pix
 * @param pixel_scale Pixel scale in arcsec/pixel
 * @return
 */
inline flt_t get_image_area( flt_t const & image_size_x_pix, flt_t const & image_size_y_pix, flt_t const & pixel_scale )
{
	return image_size_x_pix * image_size_y_pix * IceBRG::square( pixel_scale/60. );
}

inline flt_t get_zp( flt_t const & inst_zp, flt_t const & exp_time)
{
	return inst_zp + 2.5 * std::log10(exp_time);
}

inline flt_t generate_count( flt_t const & lambda, gen_t & rng )
{
	return IceBRG::Pois_rand( lambda , rng);
}

// Non-inlined functions

// Image-level

flt_t get_background_noise( flt_t const & subtracted_background, flt_t const & unsubtracted_background,
		flt_t const & read_noise, flt_t const & gain, flt_t const & pixel_scale );

// Cluster-level

flt_t generate_cluster_mass( flt_t const & cluster_redshift, gen_t & _rng );

flt_t get_cluster_richness( flt_t const & cluster_mass, flt_t const & cluster_redshift );

// Galaxy-level

flt_t generate_abs_mag_vis(flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & cluster_mass, gen_t & rng);

flt_t get_apparent_mag_vis( flt_t const & absolute_mag_vis, flt_t const & redshift );

flt_t get_stellar_mass( flt_t const & absolute_mag_vis );

/**
 * Generates a random physical bulge size in kpc from galaxy type, redshift, and stellar mass in M_sun.
 *
 * @param galaxy_type Value indicating whether the galaxy is a field galaxy, central, or satellite.
 * @param redshift Redshift
 * @param stellar_mass Stellar mass of the galaxy in M_sun
 * @param rng Random-number generating engine
 *
 * @return Physical size of the bulge in kpc
 */
flt_t generate_physical_size_bulge( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, gen_t & rng );
/**
 * Generates a random physical disk size in kpc from galaxy type, redshift, and stellar mass in M_sun.
 *
 * @param galaxy_type Value indicating whether the galaxy is a field galaxy, central, or satellite.
 * @param redshift Redshift
 * @param stellar_mass Stellar mass of the galaxy in M_sun
 * @param rng Random-number generating engine
 *
 * @return Physical size of the disk in kpc
 */
flt_t generate_physical_size_disk( flt_t const & galaxy_type, flt_t const & redshift,
		flt_t const & stellar_mass, gen_t & rng );

flt_t generate_rotation( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & sersic_index, flt_t const & stellar_mass, gen_t & rng );

flt_t generate_rp( flt_t const & galaxy_type, flt_t const & cluster_mass, flt_t const & cluster_redshift,
		flt_t const & pixel_scale, gen_t & rng );

flt_t generate_tilt( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & sersic_index, flt_t const & stellar_mass, gen_t & rng );

flt_t generate_shear_angle( flt_t const & xp, flt_t const & yp, gen_t & rng );

flt_t generate_shear_magnitude( flt_t const & xp, flt_t const & yp, flt_t const & redshift, gen_t & rng );

flt_t generate_stellar_mass( flt_t const & galaxy_type, flt_t const & redshift, gen_t & rng );

flt_t generate_xp( flt_t const & rp, flt_t const & theta_sat, flt_t const & cluster_xp, gen_t & rng );

flt_t generate_yp( flt_t const & rp, flt_t const & theta_sat, flt_t const & cluster_yp, gen_t & rng );

flt_t generate_apparent_size_bulge( flt_t const & apparent_mag_vis, gen_t & rng );

flt_t generate_apparent_size_disk( flt_t const & apparent_mag_vis, gen_t & rng );

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_MISC_DEPENDENCIES_HPP_
