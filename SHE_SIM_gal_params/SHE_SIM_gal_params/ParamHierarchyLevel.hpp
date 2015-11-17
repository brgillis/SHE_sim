/**********************************************************************\
 @file ParamHierarchyLevel.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMHIERARCHYLEVEL_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMHIERARCHYLEVEL_HPP_

#include <unordered_map>
#include <memory>
#include <set>
#include <utility>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"

namespace SHE_SIM
{

/**
 * An abstract base class template for a level in the hierarchy of parameter generation (eg. per-image, per-galaxy, etc.)
 */
template< int_t HierachyLevel >
class ParamHierarchyLevel
{
public:

	// Public typedefs

	typedef ParamHierarchyLevel<HierachyLevel-1> * parent_ptr_t;

	typedef std::unique_ptr<ParamHierarchyLevel<HierachyLevel-1>> child_ptr_t;
	typedef std::vector<child_ptr_t> children_t;

	typedef ParamGenerator<HierachyLevel> param_t;
	typedef typename param_t::name_t param_name_t;
	typedef std::unique_ptr<param_t> param_ptr_t;
	typedef std::unordered_map<param_name_t, param_ptr_t> params_t;

	typedef typename param_t::generation_level_map_t generation_level_map_t;

private:

	// Private members
	parent_ptr_t _p_parent;
	children_t _children;
	params_t _params;
	const generation_level_map_t & _generation_level_map;

public:

	/**
	 * Constructor, creates with no children and (by default) no parent.
	 *
	 * @param p_parent Pointer to the parent of this object, defaulting to nullptr
	 */
	ParamHierarchyLevel(parent_ptr_t const & p_parent = nullptr,
			const generation_level_map_t & generation_level_map = p_parent->get_generation_level_map())
	: _p_parent(p_parent),
	  _generation_level_map(generation_level_map)
	{
	}

	/**
	 * Copy constructor is deleted. There's no logical way to handle the question of what
	 * to do if a child is copied. Do you copy the parent as well?
	 *
	 * @param other
	 */
	ParamHierarchyLevel(const ParamHierarchyLevel<HierachyLevel> & other) = delete;

	/**
	 * Move constructor.
	 *
	 * @param other
	 */
	ParamHierarchyLevel(ParamHierarchyLevel<HierachyLevel> && other)
	: _p_parent(std::move(other._p_parent)),
	  _children(std::move(other._children)),
	  _generation_level_map(other._generation_level_map)
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

	/**
	 * Copy assigmnet is deleted. There's no logical way to handle the question of what
	 * to do if a child is copied. Do you copy the parent as well?
	 *
	 * @param other
	 */
	ParamHierarchyLevel & operator=(const ParamHierarchyLevel<HierachyLevel> & other) = delete;

	/**
	 * Move assignment.
	 *
	 * @param other
	 */
	ParamHierarchyLevel & operator=(ParamHierarchyLevel<HierachyLevel> && other)
	: _p_parent(std::move(other._p_parent)),
	  _children(std::move(other._children))
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

		return *this;
	}

	/**
	 * Deconstructor - virtual to ensure this isn't instantiated itself.
	 */
	virtual ~ParamHierarchyLevel() = 0;

	// Public class methods

	/**
	 * Get the hierarchy level for this class.
	 * @return The hierachy level. 0 = highest, 1 = just below 0, etc.
	 */
	static int_t get_hierarchy_level() const
	{
		return HierachyLevel;
	}

	// Public methods

	/**
	 * Get the number of children of this object.
	 *
	 * @return Number of children.
	 */
	int_t num_children() const
	{
		return _children.size();
	}

	/**
	 * Get a pointer to this object's parent.
	 *
	 * @return A pointer to this object's parent.
	 */
	parent_ptr_t const & get_parent() noexcept
	{
		return _p_parent;
	}

	/**
	 * Get a pointer to this object's parent.
	 *
	 * @return A pointer to this object's parent.
	 */
	const parent_ptr_t const & get_parent() const noexcept
	{
		return _p_parent;
	}

	/**
	 * Get a vector of this object's children.
	 *
	 * @return A vector of this object's children.
	 */
	children_t const & get_children() noexcept
	{
		return _children;
	}

	/**
	 * Get a vector of this object's children.
	 *
	 * @return A vector of this object's children.
	 */
	const children_t const & get_children() const noexcept
	{
		return _children;
	}

	/**
	 * Get a pointer to a specific child. Will throw an exception if no child with that index exists.
	 *
	 * @param i Index of the desired child.
	 *
	 * @return A pointer to the desired child.
	 */
	child_ptr_t const & get_child(const int & i)
	{
		return _children.at(i);
	}

	/**
	 * Get a pointer to a specific child. Will throw an exception if no child with that index exists.
	 *
	 * @param i Index of the desired child.
	 *
	 * @return A pointer to the desired child.
	 */
	const child_ptr_t const & get_child(const int & i) const
	{
		return _children.at(i);
	}

	/**
	 * Get a pointer to the parameter generator with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter generator.
	 * @return Pointer to the the desired parameter generator.
	 */
	param_ptr_t const & get_param(const param_name_t & name)
	{
		return _params.at(name);
	}

	/**
	 * Get a pointer to the parameter generator with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter generator.
	 * @return Pointer to the the desired parameter generator.
	 */
	const param_ptr_t const & get_param(const param_name_t & name) const
	{
		return _params.at(name);
	}

	/**
	 * Get the value for a parameter with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter.
	 * @return The value of the desired parameter.
	 */
	param_t const & get_param_value(const param_name_t & name)
	{
		return _params.at(name)->get();
	}

	/**
	 * Get the generation level map used by this object.
	 *
	 * @return The generation level map.
	 */
	generation_level_map_t const & get_generation_level_map() const noexcept
	{
		return _generation_level_map;
	}


};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMHIERARCHYLEVEL_HPP_
