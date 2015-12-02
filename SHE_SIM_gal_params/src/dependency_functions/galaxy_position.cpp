/**********************************************************************\
 @file galaxy_position.cpp
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

#include <cmath>
#include <stdexcept>

#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"

namespace SHE_SIM {

flt_t generate_rp( const flt_t & galaxy_type, const flt_t & cluster_mass, const flt_t & cluster_redshift, const gen_t & rng  )
{
	throw std::logic_error("generate_rp NYI");
}

flt_t generate_xp( const flt_t & rp, const flt_t & theta_sat, const flt_t & cluster_xp, gen_t & rng  )
{
	return cluster_xp + rp * std::cos(theta_sat*M_PI/180);
}

flt_t generate_yp( const flt_t & rp, const flt_t & theta_sat, const flt_t & cluster_yp, gen_t & rng  )
{
	return generate_xp( rp, 45-theta_sat, cluster_yp, rng);
}

}


