/**********************************************************************\
 @file ParamGenerator.cpp
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

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>
#include <limits>

#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

#define UNCACHED_VALUE std::numeric_limits<flt_t>::infinity()

namespace SHE_SIM
{

// Protected methods

flt_t ParamGenerator::_request_param_value(name_t const & param_name)
{
	return _owner._request_param_value(param_name, name());
}

void ParamGenerator::_generate()
{
	if(_params->get_mode()==ParamParam::INDEPENDENT)
	{
		_cached_value = _params->get_independently(_rng);
	}
	else
	{
		throw bad_mode_error(_params->get_mode_name());
	}
}

// Private methods

bool ParamGenerator::_is_cached() const
{
	return _cached_value != UNCACHED_VALUE;
}

void ParamGenerator::_clear_cache()
{
	_cached_value = UNCACHED_VALUE;

	// Uncache for any children
	for( auto const & child : _owner._children )
	{
		child->_clear_param_cache(name());
	}

	// Uncache any dependants as well
	for( auto const & dependant_name : _dependant_names )
	{
		_owner._clear_param_cache(dependant_name);
	}

	_dependant_names.clear();
}

void ParamGenerator::_add_dependant(name_t const & dependant_name)
{
	_dependant_names.insert(dependant_name);
}

bool ParamGenerator::_generated_at_this_level() const
{
	return _owner.get_hierarchy_level() <= level_generated_at();
}

void ParamGenerator::_determine_value()
{
	if(_generated_at_this_level())
	{
		_generate();
	}
	else
	{
		_cached_value = _parent_version().get();
	}
}

void ParamGenerator::_determine_new_value()
{
	_clear_cache();
	this->_determine_value();
}

ParamGenerator & ParamGenerator::_parent_version()
{
	return *(_owner.get_parent()->get_param(name()));
}

const ParamGenerator & ParamGenerator::_parent_version() const
{
	return *(_owner.get_parent()->get_param(name()));
}

ParamGenerator::ParamGenerator( owner_t & owner, level_t const * const & p_generation_level )
: _cached_value(UNCACHED_VALUE),
  _owner(owner),
  _params(nullptr),
  _generation_level(p_generation_level),
  _rng(_owner._rng)
{
}

void ParamGenerator::set_p_params(ParamParam const * const & p)
{
	_clear_cache();
	_params = p;
}

ParamParam const & ParamGenerator::get_params() const
{
	return *_params;
}

ParamParam const * const & ParamGenerator::get_p_params() const noexcept
{
	return _params;
}

level_t const & ParamGenerator::get_generation_level() const
{
	return *_generation_level;
}

level_t const * const & ParamGenerator::get_p_generation_level() const
{
	return _generation_level;
}

void ParamGenerator::set_generation_level( level_t const & level )
{
	_owner.set_generation_level(name(),level);
}

void ParamGenerator::set_p_generation_level( level_t const * const & p_level )
{
	_clear_cache();
	_generation_level = p_level;
}

const flt_t & ParamGenerator::get()
{
	if(!_is_cached())
	{
		_determine_value();
	}
	return _cached_value;
}

const flt_t & ParamGenerator::get_new()
{
	_determine_new_value();
	return _cached_value;
}

const flt_t & ParamGenerator::request(name_t const & requester_name)
{
	_add_dependant(requester_name);
	return get();
}

const flt_t & ParamGenerator::request_new(name_t const & requester_name)
{
	_add_dependant(requester_name);
	return get_new();
}

const level_t & ParamGenerator::level_generated_at() const
{
	return _owner.get_generation_level(name());
}

} // namespace SHE_SIM
