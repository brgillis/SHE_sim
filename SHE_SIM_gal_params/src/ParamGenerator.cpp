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

#include <limits>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"

#define UNCACHED_VALUE std::numeric_limits<flt_t>::infinity()

namespace SHE_SIM
{

// Protected methods

flt_t ParamGenerator::_request_param_value(const name_t & param_name)
{
	return _owner._request_param_value(param_name, name());
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

void ParamGenerator::_add_dependant(const name_t & dependant_name)
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

ParamGenerator::ParamGenerator( owner_t & owner, const int_t & level_generated_at )
: _cached_value(UNCACHED_VALUE),
  _owner(owner)
{
}

void ParamGenerator::set_params(const std::vector<flt_t> & v)
{
	_clear_cache();
	_set_params(v);
}

void ParamGenerator::set_params(std::vector<flt_t> && v)
{
	_clear_cache();
	_set_params(std::move(v));
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

const flt_t & ParamGenerator::request(const name_t & requester_name)
{
	_add_dependant(requester_name);
	return get();
}

const flt_t & ParamGenerator::request_new(const name_t & requester_name)
{
	_add_dependant(requester_name);
	return get_new();
}

const int_t & ParamGenerator::level_generated_at() const
{
	return _owner.get_generation_level(name());
}

} // namespace SHE_SIM
