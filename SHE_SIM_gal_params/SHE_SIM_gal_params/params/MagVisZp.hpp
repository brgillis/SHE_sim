/**********************************************************************\
 @file MagVisZp.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_MAGVISZP_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_MAGVISZP_HPP_

#include <cassert>
#include <cmath>
#include <vector>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/default_values.h"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class MagVisZp : public ParamGenerator
{
private:

	virtual void _generate() override
	{
		ParamGenerator::_cached_value = _owner.get_param_value("mag_vis_inst_zp")
				+ 2.5* std::log10(_owner.get_param_value("exp_time"));
	}

	virtual void _set_params(const std::vector<flt_t> & v) override
	{
		assert(v.size()==0);
	}

public:
	MagVisZp( owner_t & owner)
	: ParamGenerator(owner)
	{
	}

	virtual ~MagVisZp()
	{
	}

	virtual ParamGenerator::name_t name() const override
	{
		return "mag_vis_zp";
	}

	virtual ParamGenerator * clone() const override
	{
		return new MagVisZp(*this);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_MAGVISZP_HPP_
