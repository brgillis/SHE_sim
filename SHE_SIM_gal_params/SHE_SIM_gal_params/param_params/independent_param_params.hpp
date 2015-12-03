/**********************************************************************\
 @file IndContRayleigh.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDEPENDENT_PARAM_PARAMS_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDEPENDENT_PARAM_PARAMS_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <stdexcept>
#include <initializer_list>

#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

// Define a macro for an independent ParamParam
#define INDEPENDENT_PARAM_PARAM( name, get_method ) \
class name##_pp: public ParamParam \
{ \
private: \
\
	\/\/ Private members \
	flt_t _sigma, _max, _p; \
\
	\/\/ Private methods \
	virtual bool is_equal( ParamParam const * const & other ) const override \
	{ \
		name##_pp const * other_derived = dynamic_cast<name##_pp const *>(other); \
		return other_derived != nullptr; \
	} \
\
public: \
\
	\/\/ Constructor and destructor \
	name##_pp( ) \
	: ParamParam(ParamParam::INDEPENDENT), \
	  _sigma(sigma), \
	  _max(max), \
	  _p(p) \
	{ \
	} \
	virtual ~name##_pp() {} \
\
	\/\/ Get the name of this \
	virtual name_t name() const override { return "independent_" #name; }; \
\
	\/\/ Get the value \
	virtual flt_t get_independently( gen_t & gen = rng ) const override \
	{ \
		get_method; \
	} \
\
	virtual ParamParam * clone() const override \
	{ \
		return new IndContRayleigh(*this); \
	} \
\
	virtual ParamParam * recreate(const std::initializer_list<flt_t> & params) const override \
	{ \
		return new IndContRayleigh(); \
	} \
}; \

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDEPENDENT_PARAM_PARAMS_HPP_
