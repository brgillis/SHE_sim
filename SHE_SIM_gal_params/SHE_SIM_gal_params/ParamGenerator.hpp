/**********************************************************************\
 @file ParamGenerator.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMGENERATOR_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMGENERATOR_HPP_

#include <limits>
#include <unordered_map>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"

#define UNCACHED_VALUE std::numeric_limits<flt_t>::infinity()

namespace SHE_SIM
{

/**
 * An abstract base class representing a generator for a certain type of parameter.
 */
template< int_t HierarchyLevel >
class ParamGenerator
{
public:
	// Public typedefs
	typedef ParamHierarchyLevel<HierarchyLevel> owner_t;
	typedef std::string name_t;
	typedef std::unordered_map<name_t,int_t> generation_level_map_t;

private:

	// Private members
	owner_t & _owner;
	const int_t _level_generated_at;
	flt_t _cached_value;

	// Private methods

	bool _is_cached() const
	{
		return _cached_value == UNCACHED_VALUE;
	}

	void _clear_cache()
	{
		_cached_value = UNCACHED_VALUE;
	}

	virtual void _generate() = 0;

	bool _generated_at_this_level() const
	{
		return HierarchyLevel <= _level_generated_at;
	}

	void _determine_value()
	{
		if(_generated_at_this_level())
		{
			_generate();
		}
		else
		{
			_cached_value = _parent_version()->get();
		}
	}

	void _determine_new_value()
	{
		_clear_cache();
		this->_determine_value();
	}

	ParamGenerator<HierarchyLevel-1> & _parent_version()
	{
		return _owner.parent()->param(name());
	}

	const ParamGenerator<HierarchyLevel-1> & _parent_version() const
	{
		return _owner.parent()->param(name());
	}

public:

	/**
	 * Constructor - initializes with value uncached.
	 *
	 * @param level_determined_at Which level of the hierarchy this will be generated at.
	 */
	ParamGenerator( owner_t & owner,
					const generation_level_map_t & generation_level_map = owner.get_generation_level_map() )
	: _owner(owner),
	  _cached_value(UNCACHED_VALUE)
	{
		_level_generated_at = generation_level_map.at(name());
	}

	/**
	 * Virtual destructor.
	 */
	virtual ~ParamGenerator() {}

	virtual const name_t & name() const = 0;

	const flt_t & get()
	{
		if(!_is_cached())
		{
			_generate();
		}
		return _cached_value;
	}

	const flt_t & get_new()
	{
		_determine_new_value();
		return _cached_value;
	}
};

} // namespace SHE_SIM

#undef UNCACHED_VALUE

#endif // SHE_SIM_GAL_PARAMS_PARAMGENERATOR_HPP_
