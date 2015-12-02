/**********************************************************************\
 @file galaxy_model.cpp
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

#include "SHE_SIM_gal_params/dependency_functions/object_dependencies.hpp"
#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

core_sed_t get_core_sed( flt_t const & morphology, flt_t const & redshift, flt_t const & stellar_mass )
{
	throw std::logic_error("get_core_sed NYI");
}

core_observed_flux_distribution_t get_core_observed_flux_distribution(
		flt_t const & morphology, flt_t const & rotation, flt_t const & tilt)
{
	throw std::logic_error("get_core_observed_flux_distribution NYI");
}

disk_sed_t get_disk_sed( flt_t const & morphology, flt_t const & redshift, flt_t const & stellar_mass )
{
	throw std::logic_error("get_disk_sed NYI");
}

disk_observed_flux_distribution_t get_disk_observed_flux_distribution(
		flt_t const & morphology, flt_t const & rotation, flt_t const & tilt)
{
	throw std::logic_error("get_disk_observed_flux_distribution NYI");
}

} // namespace SHE_SIM

