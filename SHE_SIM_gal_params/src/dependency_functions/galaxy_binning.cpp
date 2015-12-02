/**********************************************************************\
 @file galaxy_binning.cpp
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

binned_observed_flux_distribution_t get_binned_observed_flux_distribution(
		core_sed_t const & (core_sed),
		disk_sed_t const & (disk_sed),
		core_observed_flux_distribution_t const & core_observed_flux_distribution,
		disk_observed_flux_distribution_t const & disk_observed_flux_distribution)
{
	throw std::logic_error("get_binned_observed_flux_distribution NYI");
}

observed_flux_distribution_t get_observed_flux_distribution(
		binned_observed_flux_distribution_t const & binned_observed_flux_distribution,
		binned_psf_t const & binned_psf)
{
	throw std::logic_error("get_observed_flux_distribution NYI");
}

} // namespace SHE_SIM
