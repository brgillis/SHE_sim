/**********************************************************************\
 @file values.h
 ------------------

 Default values for galaxy simulation parameters.

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

#ifndef SHE_SIM_GAL_PARAMS_VALUES_H_
#define SHE_SIM_GAL_PARAMS_VALUES_H_

#include <SHE_SIM_gal_params/common.hpp>

namespace SHE_SIM { namespace dv {

// Survey-level

constexpr const int_t survey_level = 0;

constexpr const flt_t gain = 3.3; // e-/ADU
constexpr const flt_t num_images = 1;
constexpr const flt_t pixel_scale = 0.1; // arcsec/pixel
constexpr const flt_t read_noise = 5.4; // e-/pixel
constexpr const flt_t mag_vis_inst_zp = 25.6527; // Mag_vis instrumental zeropoint, from Lance's code
constexpr const flt_t mag_i_inst_zp = 25.3884; // Mag_i instrumental zeropoint, from Lance's code

constexpr const flt_t vis_filter_response = 1; // FIXME dummy value

// ImageGroup-level

constexpr const int_t image_group_level = 1;

// Image-level

constexpr const int_t image_level = 2;

constexpr const flt_t background_galaxy_density = 60.; // #/arcmin^2
constexpr const flt_t cluster_density = 1.; // #/arcmin^2
constexpr const flt_t field_galaxy_density = 30.; // #/arcmin^2
constexpr const flt_t star_density = 30.; // #/arcmin^2

constexpr const flt_t num_fields = 1.; // #/arcsec^2

constexpr const flt_t image_size_xp = 4096; // pixels
constexpr const flt_t image_size_yp = 2048; // pixels

constexpr const flt_t exp_time = 565.; // seconds

constexpr const flt_t sky_level = 32570.; // ADU/arcsec
constexpr const flt_t subtracted_background_l10_mean = std::log10(sky_level);
constexpr const flt_t subtracted_background_l10_stddev = 0.05;

constexpr const flt_t unsubtracted_background = 0.; // ADU/arcsec

constexpr const flt_t psf_params = 0.; // FIXME dummy value

// ClusterGroup-level

constexpr const int_t cluster_group_level = 3;

// FieldGroup-level

constexpr const int_t field_group_level = 3;

// Cluster-level

constexpr const int_t cluster_level = 4;

constexpr const flt_t cluster_mass_l10_mean = 14; // log10 Msun
constexpr const flt_t cluster_mass_l10_stddev = 0.5; // log10 Msun

constexpr const flt_t cluster_redshift_scale = 1.;
constexpr const flt_t cluster_redshift_median = 0.5;
constexpr const flt_t cluster_redshift_min = 0.05;
constexpr const flt_t cluster_redshift_max = 2.00;

constexpr const flt_t cluster_richness = 10.;

constexpr const flt_t cluster_xp_min = 0.; // pixels
constexpr const flt_t cluster_xp_max = image_size_xp; // pixels

constexpr const flt_t cluster_yp_min = 0.; // pixels
constexpr const flt_t cluster_yp_max = image_size_yp; // pixels

constexpr const flt_t cluster_num_satellites = 10;

// Field-level

constexpr const int_t field_level = 4;

// GalaxyGroup-level

constexpr const int_t galaxy_group_level = 5;

// GalaxyGroup-level

constexpr const int_t galaxy_pair_level = 6;

// Galaxy-level

constexpr const int_t galaxy_level = 7;

constexpr const flt_t apparent_mag_vis_min = 17;
constexpr const flt_t apparent_mag_vis_max = 24.5;

constexpr const flt_t apparent_size_bulge = std::pow(10.,-0.5); // arcsec

constexpr const flt_t apparent_size_bulge_l10_mean = -0.5; // log10 arcsec
constexpr const flt_t apparent_size_bulge_l10_stddev = 0.15; // log10 arcsec

constexpr const flt_t apparent_size_disk = std::pow(10.,-0.5); // arcsec

constexpr const flt_t apparent_size_disk_l10_mean = -0.5; // log10 arcsec
constexpr const flt_t apparent_size_disk_l10_stddev = 0.15; // log10 arcsec

constexpr const flt_t bulge_fraction_l10_mean = -0.5; // log10
constexpr const flt_t bulge_fraction_l10_stddev = 0.15; // log10
constexpr const flt_t bulge_fraction_l10_max = 1; // log10
constexpr const flt_t bulge_fraction_l10_min = -999; // log10

constexpr const int_t galaxy_type = 0.; // 0 = field, <0 = central, >0 = satellite

constexpr const flt_t morphology_min = 0.5; // Sersic index
constexpr const flt_t morphology_max = 6.5; // Sersic index

constexpr const flt_t physical_size_bulge = 1.; // kpc

constexpr const flt_t physical_size_disk = 1.; // kpc

constexpr const flt_t redshift_scale = 1.;
constexpr const flt_t redshift_median = 1.00;
constexpr const flt_t redshift_min = 0.05;
constexpr const flt_t redshift_max = 2.00;

constexpr const flt_t rotation_min = 0.; // degrees
constexpr const flt_t rotation_max = 180.; // degrees

constexpr const flt_t shear_angle_min = 0.; // degrees
constexpr const flt_t shear_angle_max = 180.; // degrees

constexpr const flt_t shear_magnitude_sigma = 0.25;
constexpr const flt_t shear_magnitude_max = 0.9;
constexpr const flt_t shear_magnitude_p = 4.;

constexpr const flt_t stellar_mass = 1e9; // kpc

constexpr const flt_t theta_sat_min = 0.;
constexpr const flt_t theta_sat_max = 360.;

constexpr const flt_t tilt_min = 0.;
constexpr const flt_t tilt_max = 90.;

constexpr const flt_t xp_min = 0.;
constexpr const flt_t xp_max = image_size_xp;

constexpr const flt_t yp_min = 0.;
constexpr const flt_t yp_max = image_size_yp;

} } // namespace SHE_SIM{ namespace dv{



#endif // SHE_SIM_GAL_PARAMS_VALUES_H_
