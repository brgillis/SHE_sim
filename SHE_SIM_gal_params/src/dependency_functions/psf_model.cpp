/**********************************************************************\
 @file psf_model.cpp
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

// Image-level

background_psf_t get_background_psf( flt_t const & psf_params )
{
	throw std::logic_error("get_background_psf NYI");
}

// Galaxy-level

binned_psf_t get_binned_psf(psf_model_t const & psf_model, flt_t const & xp, flt_t const & yp)
{
	throw std::logic_error("get_binned_psf NYI");
}

psf_model_t get_psf_model( flt_t const & psf_params )
{
	throw std::logic_error("get_psf_model NYI");
}

} // namespace SHE_SIM


