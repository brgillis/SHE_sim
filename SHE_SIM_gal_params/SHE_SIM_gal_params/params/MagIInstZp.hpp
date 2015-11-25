/**********************************************************************\
 @file MagIInstZp.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_MAGIINSTZP_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_MAGIINSTZP_HPP_

#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <SHE_SIM_gal_params/param_names.hpp>
#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class MagIInstZp : public ParamGenerator
{

public:
	MagIInstZp( owner_t & owner )
	: ParamGenerator(owner)
	{
		_params = default_param_params_map.at(name()).get();
	}

	virtual ~MagIInstZp()
	{
	}

	virtual name_t name() const override
	{
		return mag_i_inst_zp_name;
	}

	virtual ParamGenerator * clone() const override
	{
		return new MagIInstZp(*this);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_MAGIINSTZP_HPP_
