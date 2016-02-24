/**********************************************************************\
 @file IndRedshift.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDCLUSTERREDSHIFT_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDCLUSTERREDSHIFT_HPP_

#include <initializer_list>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/dependency_functions/galaxy_redshift.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndClusterRedshift: public ParamParam
{
private:

	flt_t _enhancement, _z_min, _z_max;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndClusterRedshift const * other_derived = dynamic_cast<IndClusterRedshift const *>(other);
		if(other_derived==nullptr) return false;
		return (_enhancement==other_derived->_enhancement)
				and (_z_min==other_derived->_z_min) and (_z_max==other_derived->_z_max);
	}

public:

	// Constructor and destructor
	IndClusterRedshift( flt_t const & enhancement = 1.,
			flt_t const & z_min = 0.2, flt_t const & z_max = 1.3
			)
	: ParamParam(ParamParam::INDEPENDENT),
	  _enhancement(enhancement),
	  _z_min(z_min),
	  _z_max(z_max)
	{
	}
	virtual ~IndClusterRedshift() {}

	// Get the name of this
	virtual name_t name() const override { return "cluster_redshift"; };

	// Get the value
	virtual flt_t get_independently( gen_t & gen = IceBRG::rng ) const override
	{
		return generate_cluster_z(_z_min,_z_max,gen);
	}

	virtual ParamParam * clone() const override
	{
		return new IndClusterRedshift(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		if(params.size() != 3) throw std::runtime_error("Invalid number of arguments used for redshift param param.\n"
				"Exactly 3 arguments are required.");
		return new IndClusterRedshift(params[0],params[1],params[2]);
	}

	// Get parameter values
	flt_t const & get_enhancement() const
	{
		return _enhancement;
	}
	flt_t const & get_z_min() const
	{
		return _z_min;
	}
	flt_t const & get_z_max() const
	{
		return _z_max;
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDCLUSTERREDSHIFT_HPP_
