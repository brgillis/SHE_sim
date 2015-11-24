/**********************************************************************\
 @file IndRayleigh.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDRAYLEIGH_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDRAYLEIGH_HPP_

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/random_functions.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndRayleigh: public ParamParam
{
private:

	flt_t _sigma;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndRayleigh const * other_derived = dynamic_cast<IndRayleigh const *>(other);
		if(other_derived==nullptr) return false;
		return (_sigma==other_derived->_sigma);
	}

public:

	// Constructor and destructor
	IndRayleigh( flt_t const & sigma = 1. )
	: ParamParam(ParamParam::INDEPENDENT),
	  _sigma(sigma)

	{
	}
	virtual ~IndRayleigh() {}

	// Get the name of this
	virtual name_t name() const override { return "rayleigh"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return Rayleigh_rand(_sigma, gen);
	}

	virtual ParamParam * clone() const override
	{
		return new IndRayleigh(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		return new IndRayleigh(params.at(0));
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDRAYLEIGH_HPP_
