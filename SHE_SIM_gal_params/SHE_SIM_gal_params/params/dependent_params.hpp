/**********************************************************************\
 @file DEPENDENT_PARAMs.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/math.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

// Define a macro for each param

#define DEPENDENT_PARAM( class_name, param_name, dependent_generation ) \
class class_name : public ParamGenerator \
{ \
private: \
\
	virtual void _generate() override \
	{ \
		if(_params->get_mode()==ParamParam::DEPENDENT) \
		{ \
			_cached_value = dependent_generation; \
		} \
		else if(_params->get_mode()==ParamParam::INDEPENDENT) \
		{ \
			_cached_value = _params->get_independently(_rng); \
		} \
		else \
		{ \
			throw bad_mode_error(_params->get_mode_name()); \
		} \
	} \
\
public: \
	class_name( owner_t & owner) \
	: ParamGenerator(owner) \
	{ \
		_params = default_param_params_map.at(name()).get(); \
	} \
\
	virtual ~class_name() \
	{ \
	} \
\
	virtual name_t name() const override \
	{ \
		return param_name##_name; \
	} \
\
	virtual ParamGenerator * clone() const override \
	{ \
		return new class_name(*this); \
	} \
};

// Define each param

// Survey level

// Image level

DEPENDENT_PARAM(BackgroundNoise, background_noise,
		0; throw std::logic_error("Dependent calculation for background_noise not yet implemented."));

DEPENDENT_PARAM(ImageArea, image_area,
		_request_param_value(image_size_xp_name) * _request_param_value(image_size_yp_name)
		* square(_request_param_value(pixel_scale_name)));

DEPENDENT_PARAM(MagIZp, mag_i_zp, _request_param_value(mag_i_inst_zp_name)
				                  + 2.5 * std::log10(_request_param_value(exp_time_name)));

DEPENDENT_PARAM(MagVisZp, mag_vis_zp, _request_param_value(mag_vis_inst_zp_name)
								      + 2.5 * std::log10(_request_param_value(exp_time_name)));

DEPENDENT_PARAM(NumBackgroundGalaxies, num_background_galaxies,
		Pois_rand( _request_param_value(image_area_name) *
		_request_param_value(background_galaxy_density_name) , _rng) );

DEPENDENT_PARAM(NumClusters, num_clusters,
		Pois_rand( _request_param_value(image_area_name) *
		_request_param_value(cluster_density_name) , _rng) );

DEPENDENT_PARAM(NumFieldGalaxies, num_field_galaxies,
		Pois_rand( _request_param_value(image_area_name) *
		_request_param_value(field_galaxy_density_name) , _rng) );

DEPENDENT_PARAM(NumStars, num_stars,
		Pois_rand( _request_param_value(image_area_name) *
		_request_param_value(star_density_name) , _rng) );

// Cluster level

DEPENDENT_PARAM(ClusterMass, cluster_mass,
		0; throw std::logic_error("Dependent calculation for cluster_mass not yet implemented."));

DEPENDENT_PARAM(ClusterNumSatellites, cluster_num_satellites,
		0; throw std::logic_error("Dependent calculation for cluster_num_satellites not yet implemented."));

// Galaxy level

DEPENDENT_PARAM(ApparentMagVis, apparent_mag_vis,
		0; throw std::logic_error("Dependent calculation for apparent_mag_vis not yet implemented."));

DEPENDENT_PARAM(ApparentSize, apparent_size,
		0; throw std::logic_error("Dependent calculation for apparent_size not yet implemented."));

DEPENDENT_PARAM(GalaxyType, galaxy_type,
		0; throw std::logic_error("Dependent calculation for galaxy_type not yet implemented."));

DEPENDENT_PARAM(Morphology, morphology,
		0; throw std::logic_error("Dependent calculation for morphology not yet implemented."));

DEPENDENT_PARAM(PhysicalSize, physical_size,
		0; throw std::logic_error("Dependent calculation for physical_size not yet implemented."));

DEPENDENT_PARAM(Redshift, redshift,
		0; throw std::logic_error("Dependent calculation for redshift not yet implemented."));

DEPENDENT_PARAM(Rotation, rotation,
		0; throw std::logic_error("Dependent calculation for rotation not yet implemented."));

DEPENDENT_PARAM(ShearAngle, shear_angle,
		0; throw std::logic_error("Dependent calculation for shear_angle not yet implemented."));

DEPENDENT_PARAM(ShearMagnitude, shear_magnitude,
		0; throw std::logic_error("Dependent calculation for shear_magnitude not yet implemented."));

DEPENDENT_PARAM(StellarMass, stellar_mass,
		0; throw std::logic_error("Dependent calculation for stellar_mass not yet implemented."));

DEPENDENT_PARAM(Tilt, tilt,
		0; throw std::logic_error("Dependent calculation for tilt not yet implemented."));

DEPENDENT_PARAM(Xp, xp,
		0; throw std::logic_error("Dependent calculation for xp not yet implemented."));

DEPENDENT_PARAM(Yp, yp,
		0; throw std::logic_error("Dependent calculation for yp not yet implemented."));

// Undef the macro
#undef DEPENDENT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_PARAMS_HPP_
