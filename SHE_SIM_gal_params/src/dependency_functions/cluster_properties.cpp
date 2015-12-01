/**********************************************************************\
 @file cluster_properties.cpp
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

flt_t generate_cluster_mass( flt_t const & cluster_redshift, gen_t & _rng )
{
	throw std::logic_error("generate_cluster_mass NYI");
}

flt_t generate_cluster_redshift( gen_t & _rng )
{
	throw std::logic_error("generate_cluster_redshift NYI");
}

flt_t get_cluster_richness( flt_t const & cluster_mass, flt_t const & cluster_redshift )
{
	throw std::logic_error("get_cluster_richness NYI");
}

} // namespace SHE_SIM
