/**********************************************************************\
 @file independent_object_params.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_OBJECT_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_OBJECT_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ObjectParamParam.hpp"

namespace SHE_SIM
{

// Define a macro for each param

#define INDEPENDENT_OBJECT_PARAM( class_name, param_name, object_type ) \
class class_name : public ParamGenerator \
{ \
private: \
\
	object_type _cached_object; \
\
	virtual void _generate() override \
	{ \
		if(_params->get_mode()==ParamParam::INDEPENDENT) \
		{ \
			auto _object_params = dynamic_cast<const ObjectParamParam<object_type> *>(_params); \
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
	const object_type & get_object() \
	{ \
		return _cached_object; \
	} \
\
	virtual ParamGenerator * clone() const override \
	{ \
		return new class_name(*this); \
	} \
};

// Define each param

// Survey level

INDEPENDENT_OBJECT_PARAM(VisFilterResponse, vis_filter_response, array_1d_t);

// Image level

INDEPENDENT_OBJECT_PARAM(BackgroundPSF, background_psf, array_2d_t);

// Cluster level

// Galaxy level

INDEPENDENT_OBJECT_PARAM(BinnedIntrinsicFluxDistribution, binned_intrinsic_flux_distribution, array_1d_array_t);
INDEPENDENT_OBJECT_PARAM(BinnedObservedFluxDistribution, binned_observed_flux_distribution, array_1d_array_t);
INDEPENDENT_OBJECT_PARAM(BinnedPSF, binned_psf, array_2d_array_t);
INDEPENDENT_OBJECT_PARAM(ObservedFluxDistribution, observed_flux_distribution, array_2d_t);
INDEPENDENT_OBJECT_PARAM(PSFModel, psf_model, array_2d_t);
INDEPENDENT_OBJECT_PARAM(SED, sed, array_1d_t);

// GalaxyDither level

INDEPENDENT_OBJECT_PARAM(PixGalaxyWPoisNoise, pix_galaxy_w_pois_noise, array_2d_t);

// Undef the macro
#undef INDEPENDENT_OBJECT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_INDEPENDENT_OBJECT_PARAMS_HPP_
