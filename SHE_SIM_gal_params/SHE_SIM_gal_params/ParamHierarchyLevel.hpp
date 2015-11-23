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

#include <cassert>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

namespace SHE_SIM
{

// Forward declare ParamGenerator
class ParamGenerator;

/**
 * An abstract base class template for a level in the hierarchy of parameter generation (eg. per-image, per-galaxy, etc.)
 */
class ParamHierarchyLevel
{
public:

	// Public typedefs

	typedef ParamHierarchyLevel parent_t;
	typedef parent_t * parent_ptr_t;

	typedef ParamHierarchyLevel child_t;
	typedef std::unique_ptr<child_t> child_ptr_t;
	typedef std::vector<child_ptr_t> children_t;

	typedef ParamGenerator param_t;
	typedef std::unique_ptr<param_t> param_ptr_t;
	typedef std::unordered_map<name_t,param_ptr_t> params_t;

	typedef ParamParam param_param_t;
	typedef std::unique_ptr<param_param_t> param_param_ptr_t;
	typedef std::unordered_map<name_t,param_param_ptr_t> param_params_t;

private:

	// Private members
	parent_ptr_t _p_parent;
	children_t _children;
	param_params_t _local_param_params;

	// Private methods
	void _update_parent(parent_ptr_t const & new_p_parent);

	// Private methods
	void _update_child(child_t * const & old_p_child, child_t * const & new_p_child);

	/**
	 * Get the value for a parameter with a given name. Will throw an exception if none
	 * by that name exists. Will record the name of the requesting parameter so it can
	 * be updated later if need be.
	 *
	 * @param name Name of the desired parameter.
	 * @param name Name of the requesting parameter.
	 *
	 * @return The value of the desired parameter.
	 */
	flt_t const & _request_param_value(const name_t & name, const name_t & requester_name);

	void _drop_local_param_param(const name_t & name);

	friend class ParamGenerator; // So ParamGenerators can access _request_param_value and _clear_param_cache

protected:

	// Protected members
	params_t _params;
	const generation_level_map_t * _generation_level_map;

	// Protected methods
	/**
	 * Clears the cache of the parameter with the specified name, for both this and all its
	 * children.
	 *
	 * @param name The name of the parameters whose cache is to be cleared.
	 */
	void _clear_param_cache(const name_t & name);

public:

	/**
	 * Constructor, creates with no children and (by default) no parent.
	 *
	 * @param p_parent Pointer to the parent of this object, defaulting to nullptr
	 * @param p_generation_level_map Pointer to the generation level map this will use.
	 * @param params Parameters map.
	 */
	ParamHierarchyLevel(parent_ptr_t const & p_parent = nullptr,
			const generation_level_map_t * p_generation_level_map = nullptr,
			params_t && params = params_t());

	/**
	 * Copy constructor. Note that the copied object will maintain a pointer
	 * to the parent (by necessity), but the parent won't automatically be updated to
	 * manage the copied object.
	 *
	 * @param other
	 */
	ParamHierarchyLevel(const ParamHierarchyLevel & other);

	/**
	 * Move constructor.
	 *
	 * @param other
	 */
	ParamHierarchyLevel(ParamHierarchyLevel && other);

	/**
	 * Copy assignment. Note that the copied object will maintain a pointer
	 * to the parent (by necessity), but the parent won't automatically be updated to
	 * manage the copied object.
	 *
	 * @param other
	 */
	ParamHierarchyLevel & operator=(const ParamHierarchyLevel & other);

	/**
	 * Move assignment.
	 *
	 * @param other
	 */
	ParamHierarchyLevel & operator=(ParamHierarchyLevel && other);

	/**
	 * Deconstructor - virtual to ensure this isn't instantiated itself.
	 */
	virtual ~ParamHierarchyLevel() {};

	// Public methods

	// Get details on this object
#if(1)

	/**
	 * Get the hierarchy level for this class.
	 * @return The hierachy level. 0 = highest, 1 = just below 0, etc.
	 */
	virtual int_t get_hierarchy_level() const = 0;

	/**
	 * Get the number of children of this object.
	 *
	 * @return Number of children.
	 */
	int_t num_children() const;

	/**
	 * Get the generation level map used by this object.
	 *
	 * @return The generation level map.
	 */
	const generation_level_map_t * get_generation_level_map() const noexcept;

#endif // Get details on this object

	// Parent-related methods
#if(1)

	/**
	 * Get a pointer to this object's parent.
	 *
	 * @return A pointer to this object's parent.
	 */
	parent_t * get_parent();

	/**
	 * Get a pointer to this object's parent.
	 *
	 * @return A pointer to this object's parent.
	 */
	parent_t const * get_parent() const;

#endif

	// Child-related methods
#if(1)

	/**
	 * Get a vector of this object's children.
	 *
	 * @return A vector of this object's children.
	 */
	children_t const & get_children() noexcept;

	/**
	 * Get a vector of this object's children.
	 *
	 * @return A vector of this object's children.
	 */
	children_t const & get_children() const noexcept;

	/**
	 * Get a pointer to a specific child. Will throw an exception if no child with that index exists.
	 *
	 * @param i Index of the desired child.
	 *
	 * @return A pointer to the desired child.
	 */
	child_t * get_child(const int & i);

	/**
	 * Get a pointer to a specific child. Will throw an exception if no child with that index exists.
	 *
	 * @param i Index of the desired child.
	 *
	 * @return A pointer to the desired child.
	 */
	const child_t * get_child(const int & i) const;

	/**
	 * Create a new child with the specified arguments passed to its constructor.
	 *
	 * @param args Arguments to be passed to the child's constructor after the pointer to this.
	 *
	 * @return Pointer to the new child.
	 */
	template< typename T_child, typename... Args >
	child_t * spawn_child(Args... args)
	{
		_children.push_back(child_ptr_t( new T_child(this, args...) ));
		return _children.back().get();
	}

	/**
	 * Create multiple new children, each with the specified arguments pass to its constructor.
	 *
	 * @param N The number of children to be created.
	 * @param args The arguments to be passed to each child's constructor after the pointer to this.
	 */
	template< typename T_child, typename... Args >
	void spawn_children( int_t const & N, Args... args)
	{
		assert(N>=0);

		for( int_t i=0; i<N; ++i )
			_children.push_back(child_ptr_t( new T_child(this, args...) ));

		return;
	}

	/**
	 * Take ownership of a pre-existing child.
	 *
	 * @param p_child Pointer to the child to take ownership of.
	 */
	void adopt_child(child_t * const & p_child);

#endif

	// Parameter-related methods
#if(1)

	/**
	 * Get a pointer to the parameter generator with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter generator.
	 * @return Pointer to the the desired parameter generator.
	 */
	const param_t * get_param(const name_t & name) const;

	/**
	 * Get a pointer to the parameter generator with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter generator.
	 * @return Pointer to the the desired parameter generator.
	 */
	param_t * get_param(const name_t & name);

	/**
	 * Get the value for a parameter with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter.
	 * @return The value of the desired parameter.
	 */
	flt_t const & get_param_value(const name_t & name);

	/**
	 * Get the level at which a parameter should be generated
	 *
	 * @param name The name of the parameter
	 *
	 * @return The level it's generated at
	 */
	const int & get_generation_level( const str_t & name) const;

	void set_p_param_params(const name_t & name, ParamParam const * const & params);

	template< typename T_pp, typename... Args >
	void set_param_params(const name_t & name, Args... args)
	{
		_local_param_params[name] = param_param_ptr_t(new T_pp(args...));
		set_p_param_params( name, _local_param_params.at(name).get() );
	}

	ParamParam const * const & get_param_params(const name_t & name) const;

#endif

	virtual ParamHierarchyLevel * clone() const = 0;

}; // ParamHierarchyLevel

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMHIERARCHYLEVEL_HPP_
