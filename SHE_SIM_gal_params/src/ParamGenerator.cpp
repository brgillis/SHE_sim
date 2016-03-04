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
#include <stdexcept>

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/default_param_params.hpp>

#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

#define UNCACHED_VALUE std::numeric_limits<flt_t>::infinity()

namespace SHE_SIM
{

// Protected methods

flt_t ParamGenerator::_request_param_value(name_t const & param_name)
{
	if(!_p_owner) throw std::logic_error("Cannot request another param value from a default ParamGenerator.");
	return _p_owner->_request_param_value(param_name, name());
}

ParamGenerator * ParamGenerator::_request_param(name_t const & param_name)
{
	if(!_p_owner) return nullptr;
	return _p_owner->_request_param(param_name, name());
}

void ParamGenerator::_generate()
{
	if(_p_params->get_mode()==ParamParam::INDEPENDENT)
	{
		_cached_value = _p_params->get_independently(get_rng());
	}
	else
	{
		throw bad_mode_error(_p_params->get_mode_name());
	}
}

// Private methods

bool ParamGenerator::_is_cached() const
{
	return _cached_value != UNCACHED_VALUE;
}

void ParamGenerator::_decache()
{
	_cached_value = UNCACHED_VALUE;
}

void ParamGenerator::_clear_cache()
{
	_decache();

	if(_p_owner)
	{
		// Uncache for any children
		for( auto const & child : _p_owner->_children )
		{
			child->_clear_param_cache(name());
		}

		// Uncache any dependants as well
		for( auto const & dependant_name : _dependant_names )
		{
			_p_owner->_clear_param_cache(dependant_name);
		}
	}

	_dependant_names.clear();
}

void ParamGenerator::_add_dependant(name_t const & dependant_name)
{
	_dependant_names.insert(dependant_name);
}

bool ParamGenerator::_generated_at_this_level() const
{
	if(!_p_owner) return true;
	return _p_owner->get_hierarchy_level() <= level_generated_at();
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

ParamGenerator * ParamGenerator::_p_parent_version()
{
	if(!_p_owner) return nullptr;
	auto p_parent = _p_owner->get_parent();
	if(!p_parent) return nullptr;
	return p_parent->get_param(name());
}

ParamGenerator const * ParamGenerator::_p_parent_version() const
{
	if(!_p_owner) return nullptr;
	auto p_parent = _p_owner->get_parent();
	if(!p_parent) return nullptr;
	return p_parent->get_param(name());
}

ParamGenerator & ParamGenerator::_parent_version()
{
	assert(_p_parent_version());
	return *_p_parent_version();
}

ParamGenerator const & ParamGenerator::_parent_version() const
{
	assert(_p_parent_version());
	return *_p_parent_version();
}

ParamGenerator::ParamGenerator( owner_t * const & p_owner )
: _cached_value(UNCACHED_VALUE),
  _p_owner(p_owner),
  _p_params(nullptr),
  _p_generation_level(nullptr)
{
}


ParamGenerator::owner_t * ParamGenerator::get_p_owner()
{
	return _p_owner;
}
ParamGenerator::owner_t const * ParamGenerator::get_p_owner() const
{
	return _p_owner;
}
ParamGenerator::owner_t & ParamGenerator::get_owner()
{
	if(!_p_owner) throw std::logic_error("Owner of ParamGenerator requested for default generator.");
	return *_p_owner;
}
ParamGenerator::owner_t const & ParamGenerator::get_owner() const
{
	if(!_p_owner) throw std::logic_error("Owner of ParamGenerator requested for default generator.");
	return *_p_owner;
}
void ParamGenerator::set_p_owner(ParamGenerator::owner_t * const & p_owner)
{
	_p_owner = p_owner;
}
void ParamGenerator::set_owner(ParamGenerator::owner_t & owner)
{
	_p_owner = &owner;
}

gen_t * ParamGenerator::get_p_rng()
{
	if(!get_p_owner()) return nullptr;
	return &(get_p_owner()->_rng);
}
gen_t const * ParamGenerator::get_p_rng() const
{
	if(!get_p_owner()) return nullptr;
	return &(get_p_owner()->_rng);
}
gen_t & ParamGenerator::get_rng()
{
	if(!get_p_owner()) throw std::logic_error("RNG of ParamGenerator requested for default generator.");
	return get_p_owner()->_rng;
}
gen_t const & ParamGenerator::get_rng() const
{
	if(!get_p_owner()) throw std::logic_error("RNG of ParamGenerator requested for default generator.");
	return get_p_owner()->_rng;
}

void ParamGenerator::set_p_params(ParamParam const * const & p)
{
	_clear_cache();
	_p_params = p;
}

ParamParam const & ParamGenerator::get_params() const
{
	return *_p_params;
}

ParamParam const * const & ParamGenerator::get_p_params() const noexcept
{
	return _p_params;
}

level_t const & ParamGenerator::get_generation_level() const
{
	return *_p_generation_level;
}

level_t const * const & ParamGenerator::get_p_generation_level() const
{
	return _p_generation_level;
}

void ParamGenerator::set_generation_level( level_t const & level )
{
	_p_owner->set_generation_level(name(),level);
}

void ParamGenerator::set_p_generation_level( level_t const * const & p_level )
{
	_clear_cache();
	_p_generation_level = p_level;
}

flt_t const & ParamGenerator::get()
{
	if(!_is_cached())
	{
		_determine_value();
	}
	return _cached_value;
}

flt_t const & ParamGenerator::get_new()
{
	_determine_new_value();
	return _cached_value;
}

ParamGenerator * ParamGenerator::request(name_t const & requester_name)
{
	_add_dependant(requester_name);
	if(!_is_cached())
	{
		_determine_value();
	}
	return this;
}

ParamGenerator * ParamGenerator::request_new(name_t const & requester_name)
{
	_add_dependant(requester_name);
	_determine_new_value();
	return this;
}

flt_t const & ParamGenerator::request_value(name_t const & requester_name)
{
	_add_dependant(requester_name);
	return get();
}

flt_t const & ParamGenerator::request_new_value(name_t const & requester_name)
{
	_add_dependant(requester_name);
	return get_new();
}

const level_t & ParamGenerator::level_generated_at() const
{
	return _p_owner->get_generation_level(name());
}

} // namespace SHE_SIM
