/**********************************************************************\
 @file level_names.hpp
 ------------------

 Names of the PHLs.

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

#ifndef SHE_SIM_GAL_PARAMS_LEVEL_NAMES_HPP_
#define SHE_SIM_GAL_PARAMS_LEVEL_NAMES_HPP_

#define DEF_NAME(level) constexpr const char * level##_name = #level;

DEF_NAME(survey);
DEF_NAME(image_group);
DEF_NAME(image);
DEF_NAME(cluster_group);
DEF_NAME(cluster);
DEF_NAME(field_group);
DEF_NAME(field);
DEF_NAME(galaxy_group);
DEF_NAME(galaxy);
DEF_NAME(galaxy_dither);

#undef DEF_NAME

#endif // SHE_SIM_GAL_PARAMS_LEVEL_NAMES_HPP_
