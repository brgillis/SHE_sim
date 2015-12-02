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

#include <SHE_SIM_gal_params/common.hpp>
#include <unordered_map>
#include <unordered_set>

#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

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

protected:

	// Protected members

	flt_t _cached_value;
	owner_t & _owner;
	std::unordered_set<name_t> _dependant_names;
	const ParamParam * _params;
	const level_t * _generation_level;
	gen_t & _rng;

	// Protected methods

	virtual void _generate();

	flt_t _request_param_value(name_t const & name);
	ParamGenerator * _request_param(name_t const & name);

	bool _generated_at_this_level() const;

	ParamGenerator & _parent_version();

private:

	// Private methods

	virtual bool _is_cached() const;

	virtual void _decache();
	void _clear_cache();

	void _add_dependant(name_t const & dependant_name);

	virtual void _determine_value();

	void _determine_new_value();

	const ParamGenerator & _parent_version() const;

	friend class ParamHierarchyLevel; // So it can access _clear_cache

public:

	/**
	 * Constructor - initializes with value uncached.
	 *
	 * @param level_determined_at Which level of the hierarchy this will be generated at.
	 */
	ParamGenerator( owner_t & owner, level_t const * const & p_generation_level = nullptr );

	/**
	 * Virtual destructor.
	 */
	virtual ~ParamGenerator() {}

	virtual name_t name() const = 0;

	ParamParam const & get_params() const;

	ParamParam const * const & get_p_params() const noexcept;

	void set_p_params(ParamParam const * const & p);

	level_t const & get_generation_level() const;

	level_t const * const & get_p_generation_level() const;

	void set_generation_level( level_t const & level );

	void set_p_generation_level( level_t const * const & p_level );

	const flt_t & get();

	const flt_t & get_new();

	const flt_t & request_value(name_t const & requester_name);

	const flt_t & request_new_value(name_t const & requester_name);

	ParamGenerator * request(name_t const & requester_name);

	ParamGenerator * request_new(name_t const & requester_name);

	const level_t & level_generated_at() const;

	virtual ParamGenerator * clone() const = 0;

}; // ParamGenerator

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMGENERATOR_HPP_
