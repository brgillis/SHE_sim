/**********************************************************************\
 @file MagVisInstZp.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_MAGVISINSTZP_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_MAGVISINSTZP_HPP_

#include <cassert>
#include <vector>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/default_values.h"
#include "SHE_SIM_gal_params/param_names.h"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class MagVisInstZp : public ParamGenerator
{
private:

	flt_t _mag_vis_inst_zp;

	virtual void _generate() override
	{
		ParamGenerator::_cached_value = _mag_vis_inst_zp;
	}

	virtual void _set_params(const std::vector<flt_t> & v) override
	{
		assert(v.size()==1);
		assert(v[0]>0.);
		_mag_vis_inst_zp = v[0];
	}

public:
	MagVisInstZp( owner_t & owner, const flt_t & mag_vis_inst_zp = dv::mag_vis_inst_zp)
	: ParamGenerator(owner),
	  _mag_vis_inst_zp(mag_vis_inst_zp)
	{
	}

	virtual ~MagVisInstZp()
	{
	}

	virtual name_t name() const override
	{
		return mag_vis_inst_zp_name;
	}

	virtual ParamGenerator * clone() const override
	{
		return new MagVisInstZp(*this);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_MAGVISINSTZP_HPP_
