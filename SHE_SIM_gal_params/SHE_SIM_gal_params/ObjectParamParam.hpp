/**********************************************************************\
 @file ObjectParamParam.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_ObjectParamParam_HPP_
#define SHE_SIM_GAL_PARAMS_ObjectParamParam_HPP_

#include <string>
#include <stdexcept>

#include <boost/algorithm/string.hpp>
#include <boost/bimap.hpp>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/random_functions.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
template< typename T_object >
class ObjectParamParam : public ParamParam
{

public:

	// Constructor and destructor

	ObjectParamParam( Mode const & mode = UNSPECIFIED )
	: ParamParam(mode)
	{
	};

	ObjectParamParam( name_t const & mode_name )
	: ParamParam(mode_name)
	{
	};

	virtual ~ObjectParamParam() {};

	// Comparisons

	virtual T_object get_object_independently( gen_t & gen=rng ) const = 0;

	virtual flt_t get_independently( gen_t & gen=rng ) const override
	{
		throw std::logic_error("get_independently should not be called for ObjectParamParam objects. Use "
				"get_object_independently instead.");
	}

};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_OBJECTPARAMPARAM_HPP_
