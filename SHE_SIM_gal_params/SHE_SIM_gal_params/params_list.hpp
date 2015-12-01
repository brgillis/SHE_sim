/**********************************************************************\
 @file params_list.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAMS_LIST_HPP_
#define SHE_SIM_GAL_PARAMS_PARAMS_LIST_HPP_

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/params/dependent_params.hpp>
#include <SHE_SIM_gal_params/params/dependent_object_params.hpp>
#include <SHE_SIM_gal_params/params/independent_object_params.hpp>
#include <SHE_SIM_gal_params/params/independent_params.hpp>
#include <memory>
#include <unordered_map>
#include <utility>

#include "SHE_SIM_gal_params/ParamGenerator.hpp"
#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"


namespace SHE_SIM {

class ParamGenerator;
class ParamHierarchyLevel;

// Function to get a list of all params

template<typename T_in, typename T_map>
void insert_param(T_map & res, ParamHierarchyLevel & owner)
{
	typename T_map::mapped_type new_ptr(new T_in(owner));
	auto name(new_ptr->name());

	res.insert(std::make_pair(std::move(name),std::move(new_ptr)));
}

#define INSERT_PARAM(param) insert_param<param>(res,owner);

inline params_t get_full_params_map(ParamHierarchyLevel & owner)
{
	params_t res;

	// Survey level
	INSERT_PARAM(ExposureTime);
	INSERT_PARAM(Gain);
	INSERT_PARAM(MagIInstZp);
	INSERT_PARAM(MagVisInstZp);
	INSERT_PARAM(PixelScale);
	INSERT_PARAM(ReadNoise);
	INSERT_PARAM(VisFilterResponse);

	INSERT_PARAM(MagIZp);
	INSERT_PARAM(MagVisZp);

	// Image level
	INSERT_PARAM(BackgroundGalaxyDensity);
	INSERT_PARAM(BackgroundNoise);
	INSERT_PARAM(BackgroundPSF);
	INSERT_PARAM(ClusterDensity);
	INSERT_PARAM(FieldGalaxyDensity);
	INSERT_PARAM(ImageSizeXp);
	INSERT_PARAM(ImageSizeYp);
	INSERT_PARAM(PSFParams);
	INSERT_PARAM(StarDensity);
	INSERT_PARAM(SubtractedBackground);
	INSERT_PARAM(UnsubtractedBackground);

	INSERT_PARAM(ImageArea);
	INSERT_PARAM(NumBackgroundGalaxies);
	INSERT_PARAM(NumClusters);
	INSERT_PARAM(NumFieldGalaxies);
	INSERT_PARAM(NumStars);

	// Cluster level
	INSERT_PARAM(ClusterMass);
	INSERT_PARAM(ClusterRedshift);
	INSERT_PARAM(ClusterNumSatellites);
	INSERT_PARAM(ClusterXp);
	INSERT_PARAM(ClusterYp);

	// Galaxy-level params here

	INSERT_PARAM(ApparentMagVis);
	INSERT_PARAM(ApparentSize);
	INSERT_PARAM(BinnedIntrinsicFluxDistribution);
	INSERT_PARAM(BinnedObservedFluxDistribution);
	INSERT_PARAM(BinnedPSF);
	INSERT_PARAM(GalaxyType);
	INSERT_PARAM(Morphology);
	INSERT_PARAM(ObservedFluxDistribution);
	INSERT_PARAM(PhysicalSize);
	INSERT_PARAM(PSFModel);
	INSERT_PARAM(Redshift);
	INSERT_PARAM(Rotation);
	INSERT_PARAM(SED);
	INSERT_PARAM(ShearAngle);
	INSERT_PARAM(ShearMagnitude);
	INSERT_PARAM(StellarMass);
	INSERT_PARAM(Tilt);
	INSERT_PARAM(Xp);
	INSERT_PARAM(Yp);

	// Galaxy-dither-level

	INSERT_PARAM(DitherXpShift);
	INSERT_PARAM(DitherYpShift);

	INSERT_PARAM(PixGalaxyWPoisNoise);

	return res;

} // params_t get_full_params_map()

#undef INSERT_PARAM

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAMS_PARAMS_LIST_HPP_
