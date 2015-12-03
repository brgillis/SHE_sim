/**********************************************************************\
 @file galaxy_type.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_GALAXY_TYPE_HPP_
#define SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_GALAXY_TYPE_HPP_

namespace SHE_SIM {

// Constant definitions for each type
constexpr const flt_t central_galaxy_type = -1.;
constexpr const flt_t field_galaxy_type = 0;
constexpr const flt_t satellite_galaxy_type = 1.;

inline bool is_central_galaxy( flt_t const & galaxy_type )
{
	return galaxy_type<0;
}

inline bool is_field_galaxy( flt_t const & galaxy_type )
{
	return galaxy_type==0;
}

inline bool is_satellite_galaxy( flt_t const & galaxy_type )
{
	return galaxy_type>0;
}

} // namespace SHE_SIM



#endif // SHE_SIM_GAL_PARAMS_DEPENDENCY_FUNCTIONS_GALAXY_TYPE_HPP_
