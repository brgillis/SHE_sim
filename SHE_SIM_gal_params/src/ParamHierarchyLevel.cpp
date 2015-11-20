/**********************************************************************\
 @file ParamHierarchyLevel.cpp
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

#include <memory>
#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/params_list.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"

namespace SHE_SIM
{

// Private methods
void ParamHierarchyLevel::_update_parent(parent_ptr_t const & new_p_parent)
{
	_p_parent = new_p_parent;
}
void ParamHierarchyLevel::_update_child(child_t * const & old_p_child, child_t * const & new_p_child)
{
	for( auto & child : _children )
	{
		if(child.get()==old_p_child)
		{
			child.reset(new_p_child);
		}
	}
}

// Public methods

ParamHierarchyLevel::ParamHierarchyLevel(parent_ptr_t const & p_parent,
		const generation_level_map_t * p_generation_level_map,
		params_t && params)
: _p_parent(p_parent),
  _params(std::move(params))
{
	if(p_generation_level_map!=nullptr)
	{
		_generation_level_map = p_generation_level_map;
	}
	else if(p_parent != nullptr)
	{
		_generation_level_map = p_parent->get_generation_level_map();
	}
}

ParamHierarchyLevel::ParamHierarchyLevel(const ParamHierarchyLevel & other)
: _p_parent(other._p_parent),
  _generation_level_map(other._generation_level_map)
{
	// Deep-copy _params and _children
	for( auto const & param_name_and_ptr : other._params )
	{
		_params.insert( std::make_pair( param_name_and_ptr.first , param_ptr_t( param_name_and_ptr.second->clone() ) ) );
	}
	for( auto const & child_ptr : other._children )
	{
		_children.push_back( child_ptr_t( child_ptr->clone() ) );
		_children.back()->_update_parent(this);
	}
}


/**
 * Move constructor.
 *
 * @param other
 */
ParamHierarchyLevel::ParamHierarchyLevel(ParamHierarchyLevel && other)
: _p_parent(std::move(other._p_parent)),
  _children(std::move(other._children)),
  _params(std::move(other._params)),
  _generation_level_map(std::move(other._generation_level_map))
{
	// Update parent's pointer to this
	if(_p_parent != nullptr)
	{
		_p_parent->_update_child(&other,this);
	}

	// Update children's pointers to this
	for( auto & child : _children )
	{
		child->_update_parent(this);
	}
}

ParamHierarchyLevel & ParamHierarchyLevel::operator=(const ParamHierarchyLevel & other)
{
	_p_parent = other._p_parent;
	_generation_level_map = other._generation_level_map;

	// Deep-copy _params and _children
	_params.clear();
	for( auto const & param_name_and_ptr : other._params )
	{
		_params.insert( std::make_pair( param_name_and_ptr.first , param_ptr_t( param_name_and_ptr.second->clone() ) ) );
	}
	_children.clear();
	for( auto const & child_ptr : other._children )
	{
		_children.push_back( child_ptr_t( child_ptr->clone() ) );
		_children.back()->_update_parent(this);
	}

	return *this;
}

ParamHierarchyLevel & ParamHierarchyLevel::operator=(ParamHierarchyLevel && other)
{
	_p_parent = std::move(other._p_parent);
	_params = std::move(other._params);
	_children = std::move(other._children);
	_generation_level_map = std::move(other._generation_level_map);

	// Update parent's pointer to this
	if(_p_parent != nullptr)
	{
		_p_parent->_update_child(&other,this);
	}

	// Update children's pointers to this
	for( auto & child : _children )
	{
		child->_update_parent(this);
	}

	return *this;
}

// Public methods

int_t ParamHierarchyLevel::num_children() const
{
	return _children.size();
}

ParamHierarchyLevel::parent_t * const & ParamHierarchyLevel::get_parent() noexcept
{
	return _p_parent;
}

const ParamHierarchyLevel::parent_t * ParamHierarchyLevel::get_parent() const
{
	return _p_parent;
}

ParamHierarchyLevel::children_t const & ParamHierarchyLevel::get_children() noexcept
{
	return _children;
}

ParamHierarchyLevel::children_t const & ParamHierarchyLevel::get_children() const noexcept
{
	return _children;
}

ParamHierarchyLevel::child_t * ParamHierarchyLevel::get_child(const int & i)
{
	return _children.at(i).get();
}

const ParamHierarchyLevel::child_t * ParamHierarchyLevel::get_child(const int & i) const
{
	return _children.at(i).get();
}

ParamHierarchyLevel::param_t * ParamHierarchyLevel::get_param(const param_name_t & name)
{
	return _params.at(name).get();
}

const ParamHierarchyLevel::param_t * ParamHierarchyLevel::get_param(const param_name_t & name) const
{
	return _params.at(name).get();
}

flt_t const & ParamHierarchyLevel::get_param_value(const param_name_t & name)
{
	return _params.at(name).get()->get();
}

const generation_level_map_t * ParamHierarchyLevel::get_generation_level_map() const noexcept
{
	return _generation_level_map;
}

const int & ParamHierarchyLevel::get_generation_level( const str_t & name) const
{
	return _generation_level_map->at(name);
}

void ParamHierarchyLevel::set_param_params(const param_name_t & name, const std::vector<flt_t> & params)
{
	get_param(name)->set_params(params);
}

void ParamHierarchyLevel::set_param_params(const param_name_t & name, std::vector<flt_t> && params)
{
	get_param(name)->set_params(std::move(params));
}

} // namespace SHE_SIM
