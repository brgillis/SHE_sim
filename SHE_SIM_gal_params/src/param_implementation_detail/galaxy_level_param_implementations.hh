/**********************************************************************\
 @file galaxy_level_param_implementations.hh
 ------------------

 Implementations for Galaxy-level parameters, separated out
 for readability.

 TODO <Insert file description here>

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

#ifndef SRC_PARAM_IMPLEMENTATION_DETAIL_GALAXY_LEVEL_PARAM_IMPLEMENTATIONS_HH_
#define SRC_PARAM_IMPLEMENTATION_DETAIL_GALAXY_LEVEL_PARAM_IMPLEMENTATIONS_HH_

IMPLEMENT_PARAM(absolute_mag_vis, dv::galaxy_level, Calculated
	,
		_cached_value = generate_abs_mag_vis(REQUEST(galaxy_type), REQUEST(redshift),
					REQUEST(cluster_mass),get_rng());
	,
		_cached_value = generate_abs_mag_vis(REQUEST(galaxy_type), REQUEST(redshift),
					REQUEST(cluster_mass),get_rng());
	);
IMPLEMENT_PARAM(apparent_mag_vis, dv::galaxy_level, Calculated
	,
		_cached_value = get_apparent_mag_vis(REQUEST(absolute_mag_vis), REQUEST(redshift));
	,
		_cached_value = get_apparent_mag_vis(REQUEST(absolute_mag_vis), REQUEST(redshift));
	);
IMPLEMENT_PARAM(apparent_size_bulge, dv::galaxy_level, Calculated
	,
		_cached_value = get_angle_from_distance(REQUEST(physical_size_bulge), REQUEST(redshift))
	,
		_cached_value = generate_apparent_size_bulge(REQUEST(apparent_mag_vis), get_rng());
	);
IMPLEMENT_PARAM(apparent_size_disk, dv::galaxy_level, Calculated
	,
		_cached_value = get_angle_from_distance(REQUEST(physical_size_disk), REQUEST(redshift));
	,
		_cached_value = generate_apparent_size_disk(REQUEST(apparent_mag_vis), get_rng());
	);
IMPLEMENT_PARAM(bulge_class, dv::galaxy_level, Calculated
	,
		_cached_value = generate_bulge_class(REQUEST(stellar_mass), REQUEST(redshift), get_rng());
	,
		_cached_value = generate_bulge_class(REQUEST(stellar_mass), REQUEST(redshift), get_rng());
	);
IMPLEMENT_PARAM(bulge_fraction, dv::galaxy_level, Calculated
	,
		_cached_value = get_bulge_fraction_from_class(REQUEST(bulge_class));
	,
		_cached_value = generate_bulge_fraction(REQUEST(apparent_mag_vis), REQUEST(sersic_index),
					get_rng());
	);
IMPLEMENT_PARAM(bulge_axis_ratio, dv::galaxy_level, Calculated
	,
		_cached_value = get_bulge_axis_ratio(REQUEST(sersic_index));
	,
		_cached_value = get_bulge_axis_ratio(REQUEST(sersic_index));
	);
IMPLEMENT_PARAM(bulge_ellipticity, dv::galaxy_level, Calculated
	,
		_cached_value = get_bulge_ellipticity(REQUEST(bulge_axis_ratio),REQUEST(tilt));
	,
		_cached_value = get_bulge_ellipticity(REQUEST(bulge_axis_ratio),REQUEST(tilt));
	);
IMPLEMENT_PARAM(galaxy_type, dv::galaxy_level, IndFixed(dv::galaxy_type)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(physical_size_bulge, dv::galaxy_level, Calculated
	,
		_cached_value = generate_physical_size_bulge(REQUEST(galaxy_type), REQUEST(redshift),
					REQUEST(stellar_mass), get_rng());
	,
		_cached_value = generate_physical_size_bulge(REQUEST(galaxy_type), REQUEST(redshift),
					REQUEST(stellar_mass), get_rng());
	);
IMPLEMENT_PARAM(physical_size_disk, dv::galaxy_level, Calculated
	,
		_cached_value = generate_physical_size_disk(REQUEST(galaxy_type), REQUEST(redshift),
					REQUEST(stellar_mass), get_rng());
	,
		_cached_value = generate_physical_size_disk(REQUEST(galaxy_type), REQUEST(redshift),
					REQUEST(stellar_mass), get_rng());
	);
IMPLEMENT_PARAM(redshift, dv::galaxy_level, DepFieldRedshift(dv::galaxy_redshift_enhancement,
		dv::galaxy_redshift_min, dv::galaxy_redshift_max)
	,
		if(is_field_galaxy(REQUEST(galaxy_type)))
		{
			const DepFieldRedshift * p_redshift_pp = dynamic_cast<const DepFieldRedshift *>(_p_params);
			if(p_redshift_pp==nullptr)
			{
				_cached_value = _p_params->get_independently(get_rng());
			}
			else
			{
				_cached_value = p_redshift_pp->get_dependently(
						_request_param(cluster_redshift_name)->get_p_params(),
						REQUEST(galaxy_density), REQUEST(cluster_density),
						get_rng());
			}
		}
		else
		{
			_cached_value = REQUEST(cluster_redshift);
		}
	,
		if(is_field_galaxy(REQUEST(galaxy_type)))
		{
			const DepFieldRedshift * p_redshift_pp = dynamic_cast<const DepFieldRedshift *>(_p_params);
			if(p_redshift_pp==nullptr)
			{
				_cached_value = _p_params->get_independently(get_rng());
			}
			else
			{
				_cached_value = p_redshift_pp->get_dependently(
						_request_param(cluster_redshift_name)->get_p_params(),
						REQUEST(galaxy_density), REQUEST(cluster_density),
						get_rng());
			}
		}
		else
		{
			_cached_value = REQUEST(cluster_redshift);
		}
	);
IMPLEMENT_PARAM(rotation, dv::galaxy_level, IndUniform(dv::rotation_min, dv::rotation_max)
	,
		if(is_satellite_galaxy(REQUEST(galaxy_type)))
			_cached_value = generate_rotation( REQUEST(xp), REQUEST(yp), REQUEST(cluster_xp), REQUEST(cluster_yp),
				 REQUEST(sersic_index), REQUEST(stellar_mass), get_rng()  );
		else
			_cached_value = _p_params->get_independently(get_rng());
	,
		if(is_satellite_galaxy(REQUEST(galaxy_type)))
			_cached_value = generate_rotation( REQUEST(xp), REQUEST(yp), REQUEST(cluster_xp), REQUEST(cluster_yp),
				 REQUEST(sersic_index), REQUEST(stellar_mass), get_rng()  );
		else
			_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(rp, dv::galaxy_level, Calculated
	,
		_cached_value = generate_rp(REQUEST(galaxy_type), REQUEST(cluster_mass),
				REQUEST(cluster_redshift), REQUEST(pixel_scale), get_rng());
	,
		_cached_value = generate_rp(REQUEST(galaxy_type), REQUEST(cluster_mass),
				REQUEST(cluster_redshift), REQUEST(pixel_scale), get_rng());
	);
IMPLEMENT_PARAM(sersic_index, dv::galaxy_level, Calculated
	,
		_cached_value = generate_sersic_index_from_bulge_class(REQUEST(bulge_class), get_rng());
	,
		_cached_value = generate_sersic_index_from_apparent_mag_vis(REQUEST(apparent_mag_vis),
				get_rng());
	);
IMPLEMENT_PARAM(shear_angle, dv::galaxy_level, IndUniform(dv::shear_angle_min,
															dv::shear_angle_max)
	,
		_cached_value = generate_shear_angle(REQUEST(xp), REQUEST(yp), get_rng());
	,
		_cached_value = generate_shear_angle(REQUEST(xp), REQUEST(yp), get_rng());
	);
IMPLEMENT_PARAM(shear_magnitude, dv::galaxy_level, IndContRayleigh(dv::shear_magnitude_sigma,
	dv::shear_magnitude_max,
	dv::shear_magnitude_p)
	,
		_cached_value = generate_shear_magnitude(REQUEST(xp), REQUEST(yp),
				REQUEST(redshift), get_rng());
	,
		_cached_value = generate_shear_magnitude(REQUEST(xp), REQUEST(yp),
				REQUEST(redshift), get_rng());
	);
IMPLEMENT_PARAM(spin, dv::galaxy_level, IndUniform(dv::spin_min, dv::spin_max)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(stellar_mass, dv::galaxy_level, Calculated
	,
		_cached_value = get_stellar_mass(REQUEST(absolute_mag_vis));
	,
		_cached_value = get_stellar_mass(REQUEST(absolute_mag_vis));
	);
IMPLEMENT_PARAM(theta_sat, dv::galaxy_level, IndUniform(dv::theta_sat_min, dv::theta_sat_max)
	,
		_cached_value = _p_params->get_independently(get_rng());
	,
		_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(tilt, dv::galaxy_level, IndUniform(dv::tilt_min, dv::tilt_max)
	,
		if(is_satellite_galaxy(REQUEST(galaxy_type)))
			_cached_value = generate_tilt( REQUEST(xp), REQUEST(yp), REQUEST(cluster_xp), REQUEST(cluster_yp),
				 REQUEST(sersic_index), REQUEST(stellar_mass), get_rng()  );
		else
			_cached_value = _p_params->get_independently(get_rng());
	,
		if(is_satellite_galaxy(REQUEST(galaxy_type)))
			_cached_value = generate_tilt( REQUEST(xp), REQUEST(yp), REQUEST(cluster_xp), REQUEST(cluster_yp),
				 REQUEST(sersic_index), REQUEST(stellar_mass), get_rng()  );
		else
			_cached_value = _p_params->get_independently(get_rng());
	);
IMPLEMENT_PARAM(xp, dv::galaxy_level, Calculated
	,
		if(is_field_galaxy(REQUEST(galaxy_type)))
			_cached_value = IceBRG::drand(0.,REQUEST(image_size_xp));
		else if(is_central_galaxy(REQUEST(galaxy_type)))
			_cached_value = REQUEST(cluster_xp);
		else
			_cached_value = get_xp(REQUEST(rp), REQUEST(theta_sat),
					REQUEST(cluster_xp));
	,
		if(is_field_galaxy(REQUEST(galaxy_type)))
			_cached_value = IceBRG::drand(0.,REQUEST(image_size_xp));
		else if(is_central_galaxy(REQUEST(galaxy_type)))
			_cached_value = REQUEST(cluster_xp);
		else
			_cached_value = get_xp(REQUEST(rp), REQUEST(theta_sat),
					REQUEST(cluster_xp));
	);
IMPLEMENT_PARAM(yp, dv::galaxy_level, Calculated
	,
		if(is_field_galaxy(REQUEST(galaxy_type)))
			_cached_value = IceBRG::drand(0.,REQUEST(image_size_yp));
		else if(is_central_galaxy(REQUEST(galaxy_type)))
			_cached_value = REQUEST(cluster_yp);
		else
			_cached_value = get_yp(REQUEST(rp), REQUEST(theta_sat),
				REQUEST(cluster_yp));
	,
		if(is_field_galaxy(REQUEST(galaxy_type)))
			_cached_value = IceBRG::drand(0.,REQUEST(image_size_yp));
		else if(is_central_galaxy(REQUEST(galaxy_type)))
			_cached_value = REQUEST(cluster_yp);
		else
			_cached_value = get_yp(REQUEST(rp), REQUEST(theta_sat),
				REQUEST(cluster_yp));
	);

#endif // SRC_PARAM_IMPLEMENTATION_DETAIL_GALAXY_LEVEL_PARAM_IMPLEMENTATIONS_HH_
