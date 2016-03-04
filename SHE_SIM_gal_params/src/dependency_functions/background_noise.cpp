/**********************************************************************\
 @file background_noise.cpp
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

#include "IceBRG_main/math/misc_math.hpp"

#include "SHE_SIM_gal_params/common.hpp"

namespace SHE_SIM {

flt_t get_background_noise( flt_t const & subtracted_background, flt_t const & unsubtracted_background,
		flt_t const & read_noise, flt_t const & gain, flt_t const & pixel_scale )
{
	flt_t background_ADU_per_arcsec = subtracted_background + unsubtracted_background;
	flt_t background_e_per_pixel = background_ADU_per_arcsec * IceBRG::square(pixel_scale) * gain;

	flt_t background_noise_e = std::sqrt( background_e_per_pixel + read_noise );
	flt_t background_noise_ADU = background_noise_e / gain;

	return background_noise_ADU;
}


} // namespace SHE_SIM


