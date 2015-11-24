/**********************************************************************\
 @file IndTruncGaussian.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTRUNCGAUSSIAN_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTRUNCGAUSSIAN_HPP_

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndTruncGaussian: public ParamParam
{
private:

	flt_t _mean, _stddev, _min, _max;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndTruncGaussian const * other_derived = dynamic_cast<IndTruncGaussian const *>(other);
		if(other_derived==nullptr) return false;
		return (_mean==other_derived->_mean) and (_stddev==other_derived->_stddev) and
				(_min==other_derived->_min) and (_max==other_derived->_max);
	}

public:

	// Constructor and destructor
	IndTruncGaussian( flt_t const & mean = 0., flt_t const & stddev = 1.,
			flt_t const & min = -1., flt_t const & max = 1.)
	: ParamParam(ParamParam::INDEPENDENT),
	  _mean(mean),
	  _stddev(stddev),
	  _min(min),
	  _max(max)
	{
	}
	virtual ~IndTruncGaussian() {}

	// Get the name of this
	virtual name_t name() const override { return "truncated_gaussian"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return trunc_Gaus_rand(_mean,_stddev,_min,_max,gen);
	}

	virtual ParamParam * clone() const override
	{
		return new IndTruncGaussian(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		return new IndTruncGaussian(params.at(0),params.at(1),params.at(2),params.at(3));
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTRUNCGAUSSIAN_HPP_
