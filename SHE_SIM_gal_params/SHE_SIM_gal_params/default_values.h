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

#include "SHE_SIM_gal_params/common.h"

namespace SHE_SIM { namespace dv {

// Survey-level

constexpr const int_t survey_level = 0;

constexpr const flt_t gain = 3.3; // e-/ADU
constexpr const flt_t pixel_scale = 0.1; // arcsec/pixel
constexpr const flt_t read_noise = 5.4; // e-/pixel
constexpr const flt_t mag_vis_inst_zp = 25.6527; // Mag_vis instrumental zeropoint, from Lance's code
constexpr const flt_t mag_i_inst_zp = 25.3884; // Mag_i instrumental zeropoint, from Lance's code

// ImageGroup-level

constexpr const int_t image_group_level = 10;

// Image-level

constexpr const int_t image_level = 20;

constexpr const flt_t exp_time = 565.; // seconds
constexpr const flt_t sky_level = 32570.; // ADU/arcsec

// ClusterGroup-level

constexpr const int_t cluster_group_level = 30;

// FieldGroup-level

constexpr const int_t field_group_level = 30;

// Cluster-level

constexpr const int_t cluster_level = 40;

// Field-level

constexpr const int_t field_level = 40;

// GalaxyGroup-level

constexpr const int_t galaxy_group_level = 50;

// Galaxy-level

constexpr const int_t galaxy_level = 60;

// GalaxyDither-level

constexpr const int_t galaxy_dither_level = 70;

} } // namespace SHE_SIM{ namespace dv{



#endif // SHE_SIM_GAL_PARAMS_VALUES_H_
