/**********************************************************************\
 @file ParamParam.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMPARAM_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMPARAM_HPP_

#include <string>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class ParamParam
{
public:

	/// Enum for mode. Deliberately not making it a scoped enum so that ints can be converted to it
	enum Mode
	{
		UNSPECIFIED = -1,
		INDEPENDENT = 0,
		DEPENDENT = 1,
		ALT_DEPENDENT = 2,
		OTHER = 3
	};

private:

	Mode _mode;

protected:

	Mode _get_mode_from_string( std::string mode_str )
	{
		boost::algorithm::to_lower(mode_str);

		if( mode_str == "independent" )
		{
			return INDEPENDENT;
		}
		else if( mode_str == "dependent" )
		{
			return DEPENDENT;
		}
		else if( mode_str == "alt_dependent" )
		{
			return ALT_DEPENDENT;
		}
		else if( mode_str == "other" )
		{
			return OTHER;
		}
		else if( mode_str == "unspecified" or mode_str == "none" )
		{
			return UNSPECIFIED;
		}
		else
		{
			throw std::runtime_error("Unrecognised mode: " + mode_str);
		}
	}

public:

	// Constructor and destructor

	ParamParam( Mode const & mode = UNSPECIFIED)
	: _mode(mode)
	{
	};
	virtual ~ParamParam() {};

	// Public methods

	Mode const & get_mode() const noexcept
	{
		return _mode;
	}

	virtual flt_t get_independently( gen_t & gen=rng ) const = 0;

	virtual name_t name() const = 0;

	virtual ParamParam * clone() const = 0;

};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMPARAM_HPP_
