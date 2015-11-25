/**********************************************************************\
 @file MagIZp.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_MAGIZP_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_MAGIZP_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <cmath>
#include <vector>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class MagIZp : public ParamGenerator
{
private:

	virtual void _generate() override
	{
		if(_params->get_mode()==ParamParam::DEPENDENT)
		{
			_cached_value = _request_param_value(mag_i_inst_zp_name)
				+ 2.5* std::log10(_request_param_value(exp_time_name));
		}
		else if(_params->get_mode()==ParamParam::INDEPENDENT)
		{
			_cached_value = _params->get_independently(_rng);
		}
		else
		{
			throw bad_mode_error(_params->get_mode_name());
		}
	}

public:
	MagIZp( owner_t & owner)
	: ParamGenerator(owner)
	{
		_params = default_param_params_map.at(name()).get();
	}

	virtual ~MagIZp()
	{
	}

	virtual name_t name() const override
	{
		return mag_i_zp_name;
	}

	virtual ParamGenerator * clone() const override
	{
		return new MagIZp(*this);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_MAGIZP_HPP_