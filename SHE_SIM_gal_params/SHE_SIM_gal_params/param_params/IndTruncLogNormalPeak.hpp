/**********************************************************************\
 @file IndTruncLogNormalPeak.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTRUNCLOGNORMALPEAK_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTRUNCLOGNORMALPEAK_HPP_

#include <initializer_list>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "IceBRG_main/math/random/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndTruncLogNormalPeak: public ParamParam
{
private:

	flt_t _l10_peak, _l10_stddev, _l10_min, _l10_max;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndTruncLogNormalPeak const * other_derived = dynamic_cast<IndTruncLogNormalPeak const *>(other);
		if(other_derived==nullptr) return false;
		return (_l10_peak==other_derived->_l10_peak) and (_l10_stddev==other_derived->_l10_stddev) and
				(_l10_min==other_derived->_l10_min) and (_l10_max==other_derived->_l10_max);
	}

public:

	// Constructor and destructor
	IndTruncLogNormalPeak( flt_t const & l10_peak = 0., flt_t const & l10_stddev = 1.,
			flt_t const & l10_min = -5, flt_t const & l10_max = 5.)
	: ParamParam(ParamParam::INDEPENDENT),
	  _l10_peak(l10_peak),
	  _l10_stddev(l10_stddev),
	  _l10_min(l10_min),
	  _l10_max(l10_max)
	{
	}
	virtual ~IndTruncLogNormalPeak() {}

	// Get the name of this
	virtual name_t name() const override { return "trunc_lognormal_peak"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = IceBRG::rng ) const override
	{
		return std::pow(10.,IceBRG::trunc_Gaus_rand(_l10_peak,_l10_stddev,_l10_min,_l10_max,gen));
	}

	virtual ParamParam * clone() const override
	{
		return new IndTruncLogNormalPeak(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		if(params.size() != 4) throw std::runtime_error("Invalid number of arguments used for trunc_lognormal_peak param param.\n"
				"Exactly 4 arguments are required.");
		return new IndTruncLogNormalPeak(params[0],params[1],params[2],params[3]);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDTRUNCLOGNORMALPEAK_HPP_
