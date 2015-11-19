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

bool ParamGenerator::_is_cached() const
{
	return _cached_value != UNCACHED_VALUE;
}

void ParamGenerator::_clear_cache()
{
	_cached_value = UNCACHED_VALUE;
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
		_generate();
	}
	return _cached_value;
}

const flt_t & ParamGenerator::get_new()
{
	_determine_new_value();
	return _cached_value;
}

const int_t & ParamGenerator::level_generated_at() const
{
	return _owner.get_generation_level(name());
}

} // namespace SHE_SIM
