/**********************************************************************\
 @file galaxy_observed_properties.cpp
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdexcept>

#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

flt_t get_apparent_mag_vis( flt_t const & stellar_mass, flt_t const & redshift )
{
	throw std::logic_error("get_apparent_mag_vis NYI");
}

flt_t get_apparent_size( flt_t const & physical_size, flt_t const & redshift )
{
	throw std::logic_error("get_apparent_size NYI");
}

flt_t generate_shear_angle( flt_t const & xp, flt_t const & yp, gen_t & rng  )
{
	throw std::logic_error("get_apparent_size NYI");
}

flt_t generate_shear_magnitude( flt_t const & xp, flt_t const & yp, flt_t const & redshift, gen_t & rng  )
{
	throw std::logic_error("get_apparent_size NYI");
}


} // namespace SHE_SIM
