/**********************************************************************\
 @file ObjectParamGenerator.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_OBJECTPARAMGENERATOR_HPP_
#define SHE_SIM_GAL_PARAMS_OBJECTPARAMGENERATOR_HPP_

#include <boost/optional.hpp>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/ParamGenerator.hpp"

namespace SHE_SIM
{

/**
 * An abstract base class representing a generator for a certain type of parameter.
 */
template< typename T_object >
class ObjectParamGenerator : public ParamGenerator
{
public:

	typedef T_object object_type;

private:

	boost::optional<object_type> _cached_object;

	virtual bool _is_cached() const override
	{
		return static_cast<bool>(_cached_object);
	}

	virtual void _decache() override
	{
		_cached_object = boost::none;
	}

	virtual void _determine_value() override
	{
		if(_generated_at_this_level())
		{
			_generate();
		}
		else
		{
			_cached_object = static_cast<ObjectParamGenerator<object_type> &>(_parent_version()).get_object();
		}
	}

public:

	ObjectParamGenerator( owner_t & owner, level_t const * const & p_generation_level = nullptr )
	: ParamGenerator(owner, p_generation_level)
	{
	}

	/**
	 * Virtual destructor.
	 */
	virtual ~ObjectParamGenerator() {}

	object_type const & get_object()
	{
		if(!_is_cached())
		{
			_determine_value();
		}
		return *_cached_object;
	}

	object_type const & get_new_object()
	{
		_determine_new_value();
		return *_cached_object;
	}

	object_type const & request_object(name_t const & requester_name)
	{
		_add_dependant(requester_name);
		return get_object();
	}

	object_type const & request_new_object(name_t const & requester_name)
	{
		_add_dependant(requester_name);
		return get_new_object();
	}

}; // ObjectParamGenerator

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_OBJECTPARAMGENERATOR_HPP_
