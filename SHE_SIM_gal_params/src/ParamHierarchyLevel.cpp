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

#include <ctime>
#include <memory>
#include <random>
#include <utility>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/params_list.hpp"
#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
#include "SHE_SIM_gal_params/ParamParam.hpp"

#include "SHE_SIM_gal_params/levels/Survey.hpp"
#include "SHE_SIM_gal_params/levels/ImageGroup.hpp"
#include "SHE_SIM_gal_params/levels/Image.hpp"
#include "SHE_SIM_gal_params/levels/ClusterGroup.hpp"
#include "SHE_SIM_gal_params/levels/Cluster.hpp"
#include "SHE_SIM_gal_params/levels/FieldGroup.hpp"
#include "SHE_SIM_gal_params/levels/Field.hpp"
#include "SHE_SIM_gal_params/levels/GalaxyGroup.hpp"
#include "SHE_SIM_gal_params/levels/Galaxy.hpp"

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

flt_t const & ParamHierarchyLevel::_request_param_value(name_t const & name, name_t const & requester_name)
{
	return _params.at(name)->request_value(requester_name);
}

ParamGenerator * ParamHierarchyLevel::_request_param(name_t const & name, name_t const & requester_name)
{
	return _params.at(name)->request(requester_name);
}

void ParamHierarchyLevel::_drop_local_param_param(name_t const & name)
{
	if ( _local_param_params.find(name) == _local_param_params.end() )
	{
		// Key isn't in the map
		return;
	}
	else
	{
		if(_local_param_params.at(name).get() == get_p_param_params(name))
		{
			// Don't drop if we're using it
			return;
		}
		else
		{
			// Unused, so drop it
			_local_param_params.erase(name);
		}
	}
}

void ParamHierarchyLevel::_drop_local_generation_level(name_t const & name)
{
	if ( _local_generation_levels.find(name) == _local_generation_levels.end() )
	{
		// Key isn't in the map
		return;
	}
	else
	{
		if(_local_generation_levels.at(name).get() == get_p_generation_level(name))
		{
			// Don't drop if we're using it
			return;
		}
		else
		{
			// Unused, so drop it
			_local_generation_levels.erase(name);
		}
	}
}

void ParamHierarchyLevel::_clear_param_cache(name_t const & name)
{
	// Clear for this
	_clear_own_param_cache(name);

	// Clear for all children
	for( auto & child : _children )
	{
		child->_clear_param_cache(name);
	}
}

void ParamHierarchyLevel::_clear_own_param_cache(name_t const & name)
{
	_params.at(name)->_clear_cache();
}

void ParamHierarchyLevel::_clear_param_cache()
{
	// Clear for this
	_clear_own_param_cache();

	// Clear for all children
	for( auto & child : _children )
	{
		child->_clear_param_cache();
	}
}

void ParamHierarchyLevel::_clear_own_param_cache()
{
	for( auto & param : _params )
	{
		param.second->_clear_cache();
	}
}

// Public methods

ParamHierarchyLevel::ParamHierarchyLevel(parent_ptr_t const & p_parent)
: _p_parent(p_parent),
  _params(get_full_params_map(*this))
{
	// Inherit parameters and generation_levels from parent if it exists
	if(p_parent)
	{
		// Get ID from the parent's number of children
		_local_ID = p_parent->num_children();
		set_seed(p_parent->get_seed());
	}
	else
	{
		// Set ID to zero
		_local_ID = 0;

		// Use default seed
		set_seed();
	}
}

ParamHierarchyLevel::ParamHierarchyLevel(const ParamHierarchyLevel & other)
: _p_parent(other._p_parent),
  _local_ID(other._local_ID),
  _seed(other._seed),
  _seed_seq(other._seed_seq),
  _rng(other._rng)
{
	// Deep-copy maps

	for( auto const & child_ptr : other._children )
	{
		_children.push_back( child_ptr_t( child_ptr->clone() ) );
		_children.back()->_update_parent(this);
	}

	for( auto const & param_name_and_ptr : other._params )
	{
		_params.insert( std::make_pair( param_name_and_ptr.first , param_ptr_t( param_name_and_ptr.second->clone() ) ) );
	}

	for( auto const & param_param_name_and_ptr : other._local_param_params )
	{
		_local_param_params.insert( std::make_pair( param_param_name_and_ptr.first,
				param_param_ptr_t( param_param_name_and_ptr.second->clone() ) ) );
	}

	for( auto const & gen_level_name_and_ptr : other._local_generation_levels )
	{
		_local_generation_levels.insert( std::make_pair( gen_level_name_and_ptr.first,
				level_ptr_t( new level_t( *(gen_level_name_and_ptr.second.get()) ) ) ) );
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
  _local_param_params(std::move(other._local_param_params)),
  _local_generation_levels(std::move(other._local_generation_levels)),
  _local_ID(std::move(other._local_ID)),
  _seed(std::move(other._seed)),
  _seed_seq(std::move(other._seed_seq)),
  _rng(std::move(other._rng)),
  _params(std::move(other._params))
{
	// Update parent's pointer to this
	if(_p_parent)
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
	_local_ID = other._local_ID;
	_seed = other._seed;
	_seed_seq = other._seed_seq;
	_rng = other._rng;

	// Deep-copy maps

	_children.clear();
	for( auto const & child_ptr : other._children )
	{
		_children.push_back( child_ptr_t( child_ptr->clone() ) );
		_children.back()->_update_parent(this);
	}

	_params.clear();
	for( auto const & param_name_and_ptr : other._params )
	{
		_params.insert( std::make_pair( param_name_and_ptr.first , param_ptr_t( param_name_and_ptr.second->clone() ) ) );
	}

	_local_param_params.clear();
	for( auto const & param_param_name_and_ptr : other._local_param_params )
	{
		_local_param_params.insert( std::make_pair( param_param_name_and_ptr.first,
				param_param_ptr_t( param_param_name_and_ptr.second->clone() ) ) );
	}

	_local_generation_levels.clear();
	for( auto const & gen_level_name_and_ptr : other._local_generation_levels )
	{
		_local_generation_levels.insert( std::make_pair( gen_level_name_and_ptr.first,
				level_ptr_t( new level_t( *(gen_level_name_and_ptr.second.get()) ) ) ) );
	}

	return *this;
}

ParamHierarchyLevel & ParamHierarchyLevel::operator=(ParamHierarchyLevel && other)
{
	_p_parent = std::move(other._p_parent);
	_params = std::move(other._params);
	_children = std::move(other._children);
	_local_param_params = std::move(other._local_param_params);
	_local_generation_levels = std::move(other._local_generation_levels);
	_local_ID = std::move(other._local_ID);
	_seed = std::move(other._seed);
	_seed_seq = std::move(other._seed_seq);
	_rng = std::move(other._rng);

	// Update parent's pointer to this
	if(_p_parent)
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

std::vector<ParamHierarchyLevel::child_t *> ParamHierarchyLevel::get_children( name_t const & type_name )
{
	std::vector<ParamHierarchyLevel::child_t *> res;

	for( auto & child : _children )
	{
		if( ( child->get_name()==type_name ) or ( type_name == "" ) ) res.push_back( child.get() );
	}

	return res;
}

std::vector<const ParamHierarchyLevel::child_t *> ParamHierarchyLevel::get_children( name_t const & type_name ) const
{
	std::vector<const ParamHierarchyLevel::child_t *> res;

	for( const auto & child : _children )
	{
		if( ( child->get_name()==type_name ) or ( type_name == "" ) ) res.push_back( child.get() );
	}

	return res;
}

std::vector<ParamHierarchyLevel::child_t *> ParamHierarchyLevel::get_descendants( name_t const & type_name )
{
	std::vector<child_t *> res;

	for( auto & child : _children )
	{
		if( ( child->get_name()==type_name ) or ( type_name == "" ) )
			res.push_back( child.get() );

		// Check if this one has any descendants of the desired type
		std::vector<child_t *> childs_descendants = child->get_descendants(type_name);
		for( auto & descendant : childs_descendants )
		{
			res.push_back(descendant);
		}
	}

	return res;
}

std::vector<const ParamHierarchyLevel::child_t *> ParamHierarchyLevel::get_descendants( name_t const & type_name ) const
{
	std::vector<const child_t *> res;

	for( const auto & child : _children )
	{
		if( ( child->get_name()==type_name ) or ( type_name == "" ) )
			res.push_back( child.get() );

		// Check if this one has any descendants of the desired type
		std::vector<child_t *> childs_descendants = child->get_descendants(type_name);
		for( auto & descendant : childs_descendants )
		{
			res.push_back(descendant);
		}
	}

	return res;
}

ParamHierarchyLevel::child_t * ParamHierarchyLevel::get_child(const int & i)
{
	return _children.at(i).get();
}

ParamHierarchyLevel::child_t const * ParamHierarchyLevel::get_child(const int & i) const
{
	return _children.at(i).get();
}

void ParamHierarchyLevel::adopt_child(child_t * const & p_child)
{
	_children.push_back( child_ptr_t(p_child) );
}

void ParamHierarchyLevel::autofill_children()
{
	fill_children();
	for( auto & child : _children )
	{
		child->autofill_children();
	}
}

// Methods to get children of specific types
#if(1)

std::vector<ImageGroup *> ParamHierarchyLevel::get_image_groups()
{
	return get_children<ImageGroup>();
}
std::vector<Image *> ParamHierarchyLevel::get_images()
{
	return get_children<Image>();
}
std::vector<ClusterGroup *> ParamHierarchyLevel::get_cluster_groups()
{
	return get_children<ClusterGroup>();
}
std::vector<Cluster *> ParamHierarchyLevel::get_clusters()
{
	return get_children<Cluster>();
}
std::vector<FieldGroup *> ParamHierarchyLevel::get_field_groups()
{
	return get_children<FieldGroup>();
}
std::vector<Field *> ParamHierarchyLevel::get_fields()
{
	return get_children<Field>();
}
std::vector<GalaxyGroup *> ParamHierarchyLevel::ParamHierarchyLevel::get_galaxy_groups()
{
	return get_children<GalaxyGroup>();
}
std::vector<Galaxy *> ParamHierarchyLevel::get_galaxies()
{
	return get_children<Galaxy>();
}
std::vector<Galaxy *> ParamHierarchyLevel::get_background_galaxies()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_background_galaxy())
				res.push_back(casted_child);
		}
	}

	return res;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_foreground_galaxies()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_foreground_galaxy())
				res.push_back(casted_child);
		}
	}

	return res;
}
Galaxy * ParamHierarchyLevel::get_central_galaxy()
{

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_central_galaxy())
				return casted_child;
		}
	}

	return nullptr;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_field_galaxies()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_field_galaxy())
				res.push_back(casted_child);
		}
	}

	return res;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_satellite_galaxies()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_satellite_galaxy())
				res.push_back(casted_child);
		}
	}

	return res;
}

#endif // Methods to get children of specific types

// Methods to get descendants of specific types
#if(1)

std::vector<ImageGroup *> ParamHierarchyLevel::get_image_group_descendants()
{
	return get_descendants<ImageGroup>();
}
std::vector<Image *> ParamHierarchyLevel::get_image_descendants()
{
	return get_descendants<Image>();
}
std::vector<ClusterGroup *> ParamHierarchyLevel::get_cluster_group_descendants()
{
	return get_descendants<ClusterGroup>();
}
std::vector<Cluster *> ParamHierarchyLevel::get_cluster_descendants()
{
	return get_descendants<Cluster>();
}
std::vector<FieldGroup *> ParamHierarchyLevel::get_field_group_descendants()
{
	return get_descendants<FieldGroup>();
}
std::vector<Field *> ParamHierarchyLevel::get_field_descendants()
{
	return get_descendants<Field>();
}
std::vector<GalaxyGroup *> ParamHierarchyLevel::ParamHierarchyLevel::get_galaxy_group_descendants()
{
	return get_descendants<GalaxyGroup>();
}
std::vector<Galaxy *> ParamHierarchyLevel::get_galaxy_descendants()
{
	return get_descendants<Galaxy>();
}
std::vector<Galaxy *> ParamHierarchyLevel::get_background_galaxy_descendants()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_background_galaxy())
				res.push_back(casted_child);
		}
		else
		{
			// Otherwise check if this one has any descendants of the desired type
			std::vector<Galaxy *> childs_descendants = child->get_descendants<Galaxy>();
			for( auto & descendant : childs_descendants )
			{
				res.push_back(descendant);
			}
		}
	}

	return res;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_foreground_galaxy_descendants()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_foreground_galaxy())
				res.push_back(casted_child);
		}
		else
		{
			// Otherwise check if this one has any descendants of the desired type
			std::vector<Galaxy *> childs_descendants = child->get_descendants<Galaxy>();
			for( auto & descendant : childs_descendants )
			{
				res.push_back(descendant);
			}
		}
	}

	return res;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_central_galaxy_descendants()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_central_galaxy())
				res.push_back(casted_child);
		}
		else
		{
			// Otherwise check if this one has any descendants of the desired type
			std::vector<Galaxy *> childs_descendants = child->get_descendants<Galaxy>();
			for( auto & descendant : childs_descendants )
			{
				res.push_back(descendant);
			}
		}
	}

	return res;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_field_galaxy_descendants()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_field_galaxy())
				res.push_back(casted_child);
		}
		else
		{
			// Otherwise check if this one has any descendants of the desired type
			std::vector<Galaxy *> childs_descendants = child->get_descendants<Galaxy>();
			for( auto & descendant : childs_descendants )
			{
				res.push_back(descendant);
			}
		}
	}

	return res;
}
std::vector<Galaxy *> ParamHierarchyLevel::get_satellite_galaxy_descendants()
{
	std::vector<Galaxy *> res;

	for( auto & child : get_children() )
	{
		Galaxy * casted_child = dynamic_cast<Galaxy *>(child.get());
		if( casted_child != nullptr )
		{
			if( casted_child->is_satellite_galaxy())
				res.push_back(casted_child);
		}
		else
		{
			// Otherwise check if this one has any descendants of the desired type
			std::vector<Galaxy *> childs_descendants = child->get_descendants<Galaxy>();
			for( auto & descendant : childs_descendants )
			{
				res.push_back(descendant);
			}
		}
	}

	return res;
}

#endif // Methods to get descendants of specific types

param_t * ParamHierarchyLevel::get_param( name_t const & name )
{
	return _params.at(name).get();
}

const param_t * ParamHierarchyLevel::get_param( name_t const & name) const
{
	return _params.at(name).get();
}

flt_t const & ParamHierarchyLevel::get_param_value( name_t name )
{
	boost::algorithm::to_lower(name);
	return _params.at(name).get()->get();
}

level_t const & ParamHierarchyLevel::get_generation_level( name_t name ) const
{
	boost::algorithm::to_lower(name);
	return get_param(name)->get_generation_level();
}

level_t const * const & ParamHierarchyLevel::get_p_generation_level( name_t const & name ) const
{
	return get_param(name)->get_p_generation_level();
}

void ParamHierarchyLevel::set_p_generation_level( name_t const & name, level_t const * const & p_level )
{
	get_param(name)->set_p_generation_level( p_level );

	// Pass this along to all children
	for( auto & child : _children )
	{
		child->set_p_generation_level( name, p_level );
	}

	_drop_local_generation_level(name);
}

void ParamHierarchyLevel::set_generation_level( name_t name, level_t const & level )
{
	boost::algorithm::to_lower(name);
	_local_generation_levels[name] = level_ptr_t( new level_t(level) );
	set_p_generation_level( name, _local_generation_levels.at(name).get() );
}

ParamParam const * const & ParamHierarchyLevel::get_p_param_params( name_t const & name ) const
{
	return get_param(name)->get_p_params();
}

void ParamHierarchyLevel::set_p_param_params( name_t const & name, ParamParam const * const & params )
{

	get_param(name)->set_p_params(params);

	// Pass this along to all children
	for( auto & child : _children )
	{
		child->set_p_param_params(name,params);
	}

	_drop_local_param_param(name);
}

void ParamHierarchyLevel::generate_parameters()
{
	// Get all parameters at this level
	for( auto & param_name_and_ptr : _params )
	{
		param_name_and_ptr.second->get();
	}

	// Generate for all children as well
	for( auto & child : _children )
	{
		child->generate_parameters();
	}

}

long_int_t ParamHierarchyLevel::get_full_ID() const
{
	long_int_t ID = get_local_ID();

	// Add the proper parent ID
	if(_p_parent)
	{
		long_int_t parent_ID = _p_parent->get_local_ID();

		// Multiply the parent ID by 256^(num levels above this)
		for( int i = _p_parent->get_hierarchy_level(); i<get_hierarchy_level(); ++i)
		{
			parent_ID *= 256;
		}
	}

	return ID;
}

std::vector<int_t> ParamHierarchyLevel::get_ID_seq() const
{
	// Append this to the parent's sequence if the parent exists
	if(_p_parent)
	{
		auto res = _p_parent->get_ID_seq();

		res.push_back(get_local_ID());

		return res;
	}
	else // Just use this one's ID
	{
		std::vector<int_t> res({get_local_ID()});

		return res;
	}
}

long_int_t ParamHierarchyLevel::get_full_seed() const
{
	// Start with the actual seed value
	long_int_t seed = get_seed();

	// Multiply it by 256^(depth of this level)
	for( int i = 0; i<get_hierarchy_level(); ++i)
	{
		seed *= 256;
	}

	// Add the full ID, so that each object will get a unique value for this function
	seed += get_full_ID();

	return seed;
}

void ParamHierarchyLevel::set_seed()
{
	set_seed(time(nullptr));
}


void ParamHierarchyLevel::set_seed( int_t const & seed )
{
	// Clear the cache
	_clear_own_param_cache();

	// Get a seed sequence
	auto seed_vec = get_ID_seq();
	seed_vec.push_back(seed);

	_seed_seq = std::seed_seq(seed_vec.begin(),seed_vec.end());
	_seed = seed;

	// Seed the generator
	_rng.seed(_seed_seq);

	// Seed all children with this
	for( auto & child : _children )
	{
		child->set_seed(seed);
	}
}

} // namespace SHE_SIM
