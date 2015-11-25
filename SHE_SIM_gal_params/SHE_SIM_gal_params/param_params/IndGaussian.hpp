/**********************************************************************\
 @file IndGaussian.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDGAUSSIAN_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDGAUSSIAN_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndGaussian: public ParamParam
{
private:

	flt_t _mean, _stddev;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndGaussian const * other_derived = dynamic_cast<IndGaussian const *>(other);
		if(other_derived==nullptr) return false;
		return (_mean==other_derived->_mean) and (_stddev==other_derived->_stddev);
	}

public:

	// Constructor and destructor
	IndGaussian( flt_t const & mean = 0., flt_t const & stddev = 1. )
	: ParamParam(ParamParam::INDEPENDENT),
	  _mean(mean),
	  _stddev(stddev)
	{
	}
	virtual ~IndGaussian() {}

	// Get the name of this
	virtual name_t name() const override { return "gaussian"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return Gaus_rand(_mean,_stddev,gen);
	}

	virtual ParamParam * clone() const override
	{
		return new IndGaussian(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		return new IndGaussian(params.at(0),params.at(1));
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDGAUSSIAN_HPP_
