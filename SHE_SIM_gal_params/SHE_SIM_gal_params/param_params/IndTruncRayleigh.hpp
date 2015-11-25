/**********************************************************************\
 @file IndTruncRayleigh.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTruncRayleigh_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTruncRayleigh_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndTruncRayleigh: public ParamParam
{
private:

	flt_t _sigma, _max;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndTruncRayleigh const * other_derived = dynamic_cast<IndTruncRayleigh const *>(other);
		if(other_derived==nullptr) return false;
		return (_sigma==other_derived->_sigma) and (_max==other_derived->_max);
	}

public:

	// Constructor and destructor
	IndTruncRayleigh( flt_t const & sigma = 1., flt_t const & max = 1. )
	: ParamParam(ParamParam::INDEPENDENT),
	  _sigma(sigma),
	  _max(max)
	{
	}
	virtual ~IndTruncRayleigh() {}

	// Get the name of this
	virtual name_t name() const override { return "truncated_Rayleigh"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return trunc_Rayleigh_rand(_sigma,_max,gen);
	}

	virtual ParamParam * clone() const override
	{
		return new IndTruncRayleigh(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		return new IndTruncRayleigh(params.at(0),params.at(1));
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTruncRayleigh_HPP_
