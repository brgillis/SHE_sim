/**********************************************************************\
 @file AltCalculated.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_ALTCALCULATED_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_ALTCALCULATED_HPP_

#include <initializer_list>

#include <SHE_SIM_gal_params/common.hpp>
#include "IceBRG_main/math/random/random_functions.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class AltCalculated: public ParamParam
{
public:

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		AltCalculated const * other_derived = dynamic_cast<AltCalculated const *>(other);
		if(other_derived==nullptr) return false;
		return true;
	}

public:

	// Constructor and destructor
	AltCalculated()
	: ParamParam(ParamParam::DEPENDENT)

	{
	}
	virtual ~AltCalculated() {}

	// Get the name of this
	virtual name_t name() const override { return "alt_calculated"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = IceBRG::rng ) const override
	{
		throw std::logic_error("AltCalculated parameters cannot use the 'get_independently' method.");
	}

	virtual ParamParam * clone() const override
	{
		return new AltCalculated(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		if(params.size() != 0) throw std::runtime_error("Invalid number of arguments used for alt_calculated param param.\n"
				"Exactly 0 arguments are required.");
		return new AltCalculated();
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_ALTCALCULATED_HPP_
