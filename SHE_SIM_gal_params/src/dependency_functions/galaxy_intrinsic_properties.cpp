/**********************************************************************\
 @file galaxy_intrinsic_properties.cpp
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

flt_t generate_morphology( flt_t const & galaxy_type, flt_t const & redshift, flt_t const & stellar_mass, gen_t & rng  )
{
	throw std::logic_error("generate_morphology NYI");
}

flt_t generate_physical_size( flt_t const & galaxy_type, flt_t const & redshift, flt_t const & stellar_mass, gen_t & rng  )
{
	throw std::logic_error("generate_physical_size NYI");
}

flt_t generate_redshift( flt_t const & galaxy_type, flt_t const & cluster_redshift, gen_t & rng  )
{
	throw std::logic_error("generate_redshift NYI");
}

flt_t generate_rotation( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & morphology, flt_t const & stellar_mass, gen_t & rng  )
{
	throw std::logic_error("generate_rotation NYI");
}

flt_t generate_tilt( flt_t const & xp, flt_t const & yp, flt_t const & cluster_xp, flt_t const & cluster_yp,
				         flt_t const & morphology, flt_t const & stellar_mass, gen_t & rng  )
{
	throw std::logic_error("generate_tilt NYI");
}

flt_t generate_stellar_mass( flt_t const & galaxy_type, flt_t const & redshift, gen_t & rng  )
{
	throw std::logic_error("generate_stellar_mass NYI");
}

} // namespace SHE_SIM


