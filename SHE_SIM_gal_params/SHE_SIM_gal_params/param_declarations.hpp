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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_DECLARATIONS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_DECLARATIONS_HPP_

#include "SHE_SIM_gal_params/common.hpp"

#include "SHE_SIM_gal_params/ParamGenerator.hpp"

namespace SHE_SIM {

#define DECLARE_PARAM(param) \
extern const name_t param##_name; \
 \
class param##_obj : public ParamGenerator \
{ \
private: \
	 \
	virtual void _generate() override; \
	 \
public: \
	 \
	param##_obj( owner_t * const & p_owner); \
	 \
	virtual ~param##_obj() \
	{ \
	} \
	 \
	virtual name_t name() const override \
	{ \
		return param##_name; \
	} \
	 \
	virtual ParamGenerator * clone() const override \
	{ \
		return new param##_obj(*this); \
	} \
}; \

// Survey level

DECLARE_PARAM(num_images);
DECLARE_PARAM(pixel_scale);

// Image level

DECLARE_PARAM(cluster_density);
DECLARE_PARAM(exp_time);
DECLARE_PARAM(galaxy_density);
DECLARE_PARAM(image_area);
DECLARE_PARAM(image_size_xp);
DECLARE_PARAM(image_size_yp);
DECLARE_PARAM(num_clusters);
DECLARE_PARAM(num_fields);
DECLARE_PARAM(subtracted_background);
DECLARE_PARAM(unsubtracted_background);

// Cluster level

DECLARE_PARAM(cluster_mass);
DECLARE_PARAM(cluster_redshift);
DECLARE_PARAM(cluster_num_satellites);
DECLARE_PARAM(cluster_xp);
DECLARE_PARAM(cluster_yp);

// Field level

DECLARE_PARAM(num_field_galaxies);

// Galaxy level

DECLARE_PARAM(absolute_mag_vis);
DECLARE_PARAM(apparent_mag_vis);
DECLARE_PARAM(apparent_size_bulge);
DECLARE_PARAM(apparent_size_disk);
DECLARE_PARAM(bulge_class);
DECLARE_PARAM(bulge_fraction);
DECLARE_PARAM(galaxy_type);
DECLARE_PARAM(physical_size_bulge);
DECLARE_PARAM(physical_size_disk);
DECLARE_PARAM(psf_model);
DECLARE_PARAM(redshift);
DECLARE_PARAM(rotation);
DECLARE_PARAM(rp);
DECLARE_PARAM(sersic_index);
DECLARE_PARAM(shear_angle);
DECLARE_PARAM(shear_magnitude);
DECLARE_PARAM(stellar_mass);
DECLARE_PARAM(theta_sat);
DECLARE_PARAM(tilt);
DECLARE_PARAM(xp);
DECLARE_PARAM(yp);

#undef DECLARE_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_DECLARATIONS_HPP_
