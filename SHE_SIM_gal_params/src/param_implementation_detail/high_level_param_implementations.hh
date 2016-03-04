/**********************************************************************\
 @file high_level_param_implementations.hh
 ------------------

 Implementations for Survey and Image-level parameters, separated out
 for readability.

 **********************************************************************

 Copyright (C) 2016 brg

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

#ifndef SRC_PARAM_IMPLEMENTATION_DETAIL_HIGH_LEVEL_PARAM_IMPLEMENTATIONS_HH_
#define SRC_PARAM_IMPLEMENTATION_DETAIL_HIGH_LEVEL_PARAM_IMPLEMENTATIONS_HH_

// Survey level
IMPLEMENT_PARAM(num_images, dv::survey_level, IndFixed(dv::num_images)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(pixel_scale, dv::survey_level, IndFixed(dv::pixel_scale)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);

// Image level


IMPLEMENT_PARAM(exp_time, dv::image_level, IndFixed(dv::exp_time)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(cluster_density, dv::image_level, IndFixed(dv::cluster_density)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(galaxy_density, dv::image_level, IndFixed(dv::galaxy_density)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(image_area, dv::image_level, Calculated
	,
	 	 _cached_value = get_image_area(REQUEST(image_size_xp),REQUEST(image_size_yp),
			 REQUEST(pixel_scale));
	,
	 	 _cached_value = get_image_area(REQUEST(image_size_xp),REQUEST(image_size_yp),
			 REQUEST(pixel_scale));
	);
IMPLEMENT_PARAM(image_size_xp, dv::image_level, IndFixed(dv::image_size_xp)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(image_size_yp, dv::image_level, IndFixed(dv::image_size_xp)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(num_clusters, dv::image_level, Calculated
	,
		_cached_value = generate_count( REQUEST(image_area) *
					 REQUEST(cluster_density) , get_rng());
	,
		_cached_value = generate_count( REQUEST(image_area) *
					 REQUEST(cluster_density) , get_rng());
	);
IMPLEMENT_PARAM(num_fields, dv::image_level, IndFixed(dv::num_fields)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(subtracted_background, dv::image_level,
		IndLogNormalMean(dv::subtracted_background_l10_mean,dv::subtracted_background_l10_stddev)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(unsubtracted_background, dv::image_level, IndFixed(dv::unsubtracted_background)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);

#endif // SRC_PARAM_IMPLEMENTATION_DETAIL_HIGH_LEVEL_PARAM_IMPLEMENTATIONS_HH_
