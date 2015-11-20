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
#include <vector>

#include "SHE_SIM_gal_params/common.h"
#include "SHE_SIM_gal_params/params_list.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"

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
	typedef str_t param_name_t;
	typedef std::unique_ptr<param_t> param_ptr_t;
	typedef std::unordered_map<param_name_t,param_ptr_t> params_t;

	typedef std::unordered_map<param_name_t,int_t> generation_level_map_t;

private:

	// Private members
	parent_ptr_t _p_parent;
	children_t _children;

	// Private methods
	void _update_parent(parent_ptr_t const & new_p_parent);

	// Private methods
	void _update_child(child_t * const & old_p_child, child_t * const & new_p_child);

protected:

	// Protected members
	params_t _params;
	const generation_level_map_t * _generation_level_map;

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
	 * Get a pointer to this object's parent.
	 *
	 * @return A pointer to this object's parent.
	 */
	parent_t * const & get_parent() noexcept;

	/**
	 * Get a pointer to this object's parent.
	 *
	 * @return A pointer to this object's parent.
	 */
	const parent_t * get_parent() const;

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
	 * Get a pointer to the parameter generator with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter generator.
	 * @return Pointer to the the desired parameter generator.
	 */
	const param_t * get_param(const param_name_t & name) const;

	/**
	 * Get a pointer to the parameter generator with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter generator.
	 * @return Pointer to the the desired parameter generator.
	 */
	param_t * get_param(const param_name_t & name);

	/**
	 * Get the value for a parameter with a given name. Will throw an exception if none
	 * by that name exists.
	 *
	 * @param name Name of the desired parameter.
	 * @return The value of the desired parameter.
	 */
	flt_t const & get_param_value(const param_name_t & name);

	/**
	 * Get the generation level map used by this object.
	 *
	 * @return The generation level map.
	 */
	const generation_level_map_t * get_generation_level_map() const noexcept;

	/**
	 * Get the level at which a parameter should be generated
	 *
	 * @param name The name of the parameter
	 *
	 * @return The level it's generated at
	 */
	const int & get_generation_level( const str_t & name) const;

	void set_param_params(const param_name_t & name, const std::vector<flt_t> & params);

	void set_param_params(const param_name_t & name, std::vector<flt_t> && params);

	virtual ParamHierarchyLevel * clone() const = 0;

}; // ParamHierarchyLevel

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMHIERARCHYLEVEL_HPP_
