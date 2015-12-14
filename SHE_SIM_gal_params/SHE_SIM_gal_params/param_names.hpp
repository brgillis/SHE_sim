/**********************************************************************\
 @file param_names.h
 ------------------

 Names of the parameters, stored as variables to help avoid errors from
 refactoring.

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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_NAMES_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_NAMES_HPP_

#define DEF_NAME(param) constexpr const char * param##_name = #param;

// Survey level

DEF_NAME(gain);
DEF_NAME(mag_vis_inst_zp);
DEF_NAME(mag_i_inst_zp);
DEF_NAME(mag_vis_zp);
DEF_NAME(mag_i_zp);
DEF_NAME(num_images);
DEF_NAME(pixel_scale);
DEF_NAME(read_noise);

// Image level

DEF_NAME(background_galaxy_density);
DEF_NAME(background_noise);
DEF_NAME(cluster_density);
DEF_NAME(exp_time);
DEF_NAME(field_galaxy_density);
DEF_NAME(image_area);
DEF_NAME(image_size_xp);
DEF_NAME(image_size_yp);
DEF_NAME(num_background_galaxies);
DEF_NAME(num_clusters);
DEF_NAME(num_field_galaxies);
DEF_NAME(num_fields);
DEF_NAME(num_stars);
DEF_NAME(psf_params);
DEF_NAME(star_density);
DEF_NAME(subtracted_background);
DEF_NAME(unsubtracted_background);

// Cluster level

DEF_NAME(cluster_mass);
DEF_NAME(cluster_redshift);
DEF_NAME(cluster_num_satellites);
DEF_NAME(cluster_xp);
DEF_NAME(cluster_yp);

// Galaxy level

DEF_NAME(apparent_mag_vis);
DEF_NAME(apparent_size);
DEF_NAME(galaxy_type);
DEF_NAME(morphology);
DEF_NAME(physical_size);
DEF_NAME(psf_model);
DEF_NAME(redshift);
DEF_NAME(rotation);
DEF_NAME(rp);
DEF_NAME(shear_angle);
DEF_NAME(shear_magnitude);
DEF_NAME(stellar_mass);
DEF_NAME(theta_sat);
DEF_NAME(tilt);
DEF_NAME(xp);
DEF_NAME(yp);

#undef DEF_NAME

#endif // SHE_SIM_GAL_PARAMS_PARAM_NAMES_HPP_
