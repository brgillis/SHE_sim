/**********************************************************************\
 @file dependent_object_params.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_OBJECT_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_OBJECT_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/object_types.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/dependency_functions/object_dependencies.hpp"
#include "SHE_SIM_gal_params/ObjectParamGenerator.hpp"
#include "SHE_SIM_gal_params/ObjectParamParam.hpp"

namespace SHE_SIM
{

// Define a macro for each param

#define DEPENDENT_OBJECT_PARAM( param_name, dependent_generation ) \
class param_name##_obj : public ObjectParamGenerator<param_name##_t> \
{ \
private: \
\
	param_name##_t _cached_object; \
\
	virtual void _generate() override \
	{ \
		if(_params->get_mode()==ParamParam::DEPENDENT) \
		{ \
			dependent_generation; \
		} \
		else if(_params->get_mode()==ParamParam::INDEPENDENT) \
		{ \
			auto _object_params = dynamic_cast<const ObjectParamParam<param_name##_t> *>(_params); \
			if(_object_params==nullptr) throw std::logic_error("This object requires an ObjectParamParam."); \
			_cached_object = _object_params->get_object_independently(_rng); \
			_cached_value = 0.; \
		} \
		else \
		{ \
			throw bad_mode_error(_params->get_mode_name()); \
		} \
	} \
\
private: \
\
public: \
	param_name##_obj( owner_t & owner) \
	: ObjectParamGenerator(owner) \
	{ \
		_params = default_param_params_map.at(name()).get(); \
	} \
\
	virtual ~param_name##_obj() \
	{ \
	} \
\
	virtual name_t name() const override \
	{ \
		return param_name##_name; \
	} \
\
	const param_name##_t & get_object() \
	{ \
		return _cached_object; \
	} \
\
	virtual ParamGenerator * clone() const override \
	{ \
		return new param_name##_obj(*this); \
	} \
};

// Define a macro to request a parameter
#define REQUEST(param) _request_param_value(param##_name)
#define REQUEST_OBJECT(param) static_cast<ObjectParamGenerator<param##_t> *>(_request_param(param##_name))->get_object()

// Define each param

// Survey level

// Image level

DEPENDENT_OBJECT_PARAM(background_psf,
		_cached_object = get_background_psf(REQUEST(psf_params)));

// Cluster level

// Galaxy level

DEPENDENT_OBJECT_PARAM(binned_observed_flux_distribution,
		_cached_object = get_binned_observed_flux_distribution(REQUEST_OBJECT(core_sed), REQUEST_OBJECT(disk_sed),
				REQUEST_OBJECT(core_observed_flux_distribution),
				REQUEST_OBJECT(disk_observed_flux_distribution)));

DEPENDENT_OBJECT_PARAM(binned_psf,
		_cached_object = get_binned_psf(REQUEST_OBJECT(psf_model), REQUEST(xp), REQUEST(yp)));

DEPENDENT_OBJECT_PARAM(core_observed_flux_distribution,
		_cached_object = get_core_observed_flux_distribution(REQUEST(morphology),
				REQUEST(rotation), REQUEST(tilt)));

DEPENDENT_OBJECT_PARAM(core_sed,
		_cached_object = get_core_sed(REQUEST(morphology), REQUEST(redshift), REQUEST(stellar_mass)));

DEPENDENT_OBJECT_PARAM(disk_observed_flux_distribution,
		_cached_object = get_disk_observed_flux_distribution(REQUEST(morphology),
				REQUEST(rotation), REQUEST(tilt)));

DEPENDENT_OBJECT_PARAM(disk_sed,
		_cached_object = get_disk_sed(REQUEST(morphology), REQUEST(redshift), REQUEST(stellar_mass)));

DEPENDENT_OBJECT_PARAM(observed_flux_distribution,
		_cached_object = get_observed_flux_distribution(REQUEST_OBJECT(binned_observed_flux_distribution),
				REQUEST_OBJECT(binned_psf)));

DEPENDENT_OBJECT_PARAM(psf_model,
		_cached_object = get_psf_model(REQUEST(psf_params)));

// GalaxyDither level

DEPENDENT_OBJECT_PARAM(pix_galaxy_w_pois_noise,
		_cached_object = get_pix_galaxy_w_pois_noise(REQUEST_OBJECT(observed_flux_distribution),
				REQUEST(xp), REQUEST(yp), REQUEST(pixel_scale), REQUEST(gain)));

// Undef the macro
#undef DEPENDENT_OBJECT_PARAM
#undef REQUEST
#undef REQUEST_OBJECT

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_DEPENDENT_OBJECT_PARAMS_HPP_
