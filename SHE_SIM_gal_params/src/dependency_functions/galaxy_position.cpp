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

#include <boost/log/trivial.hpp>
#include <cmath>

#include "SHE_SIM_gal_params/dependency_functions/regular_dependencies.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"
#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

flt_t generate_rp( flt_t const & galaxy_type, flt_t const & cluster_mass, flt_t const & cluster_redshift, gen_t & rng  )
{
	BOOST_LOG_TRIVIAL(warning) << "Dummy function generate_rp used.";

	return drand(0.,10.,rng);
}

flt_t generate_xp( flt_t const & rp, flt_t const & theta_sat, flt_t const & cluster_xp, gen_t & rng  )
{
	return cluster_xp + rp * std::cos(theta_sat*M_PI/180);
}

flt_t generate_yp( flt_t const & rp, flt_t const & theta_sat, flt_t const & cluster_yp, gen_t & rng  )
{
	return generate_xp( rp, 45-theta_sat, cluster_yp, rng);
}

}


