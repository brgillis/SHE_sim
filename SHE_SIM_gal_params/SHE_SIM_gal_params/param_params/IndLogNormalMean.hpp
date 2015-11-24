/**********************************************************************\
 @file IndLogNormalMean.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDLOGNORMALMEAN_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDLOGNORMALMEAN_HPP_

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndLogNormalMean: public ParamParam
{
private:

	flt_t _l10_mean, _l10_stddev;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndLogNormalMean const * other_derived = dynamic_cast<IndLogNormalMean const *>(other);
		if(other_derived==nullptr) return false;
		return (_l10_mean==other_derived->_l10_mean) and (_l10_stddev==other_derived->_l10_stddev);
	}

public:

	// Constructor and destructor
	IndLogNormalMean( flt_t const & l10_mean = 0., flt_t const & l10_stddev = 1. )
	: ParamParam(ParamParam::INDEPENDENT),
	  _l10_mean(l10_mean),
	  _l10_stddev(l10_stddev)
	{
	}
	virtual ~IndLogNormalMean() {}

	// Get the name of this
	virtual name_t name() const override { return "lognormal_peak"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return std::pow(10.,Gaus_rand(_l10_mean,_l10_stddev,gen));
	}

	virtual ParamParam * clone() const override
	{
		return new IndLogNormalMean(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		return new IndLogNormalMean(params.at(0),params.at(1));
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDLOGNORMALMEAN_HPP_
