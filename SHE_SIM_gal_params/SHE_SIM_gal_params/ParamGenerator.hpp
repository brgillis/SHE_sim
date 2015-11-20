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

#include <unordered_map>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"

namespace SHE_SIM
{

// Forward-declare ParamHierarchyLevel
class ParamHierarchyLevel;

/**
 * An abstract base class representing a generator for a certain type of parameter.
 */
class ParamGenerator
{
public:
	// Public typedefs
	typedef ParamHierarchyLevel owner_t;
	typedef str_t name_t;
	typedef std::unordered_map<name_t,int_t> generation_level_map_t;

protected:

	// Protected members
	flt_t _cached_value;

private:

	// Private members
	owner_t & _owner;

	// Private methods

	bool _is_cached() const;

	void _clear_cache();

	virtual void _generate() = 0;

	bool _generated_at_this_level() const;

	void _determine_value();

	void _determine_new_value();

	virtual void _set_params(const std::vector<flt_t> & v) = 0;

	virtual void _set_params(std::vector<flt_t> && v) {set_params(v);};

	ParamGenerator & _parent_version();

	const ParamGenerator & _parent_version() const;

public:

	/**
	 * Constructor - initializes with value uncached.
	 *
	 * @param level_determined_at Which level of the hierarchy this will be generated at.
	 */
	ParamGenerator( owner_t & owner, const int_t & level_generated_at = -1 );

	/**
	 * Virtual destructor.
	 */
	virtual ~ParamGenerator() {}

	virtual name_t name() const = 0;

	void set_params(const std::vector<flt_t> & v);

	void set_params(std::vector<flt_t> && v);

	const flt_t & get();

	const flt_t & get_new();

	const int_t & level_generated_at() const;

	virtual ParamGenerator * clone() const = 0;

}; // ParamGenerator

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMGENERATOR_HPP_
