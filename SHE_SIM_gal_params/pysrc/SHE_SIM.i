/**********************************************************************\
 @file ParamHierarchyLevel.i
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

 %module SHE_SIM
 %{
	 
	/* Includes the header in the wrapper code */
	#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
	#include "SHE_SIM_gal_params/levels/Galaxy.hpp"
	#include "SHE_SIM_gal_params/levels/Survey.hpp"
	 
 %}
 
/* Parse the header file to generate wrappers */
%include "../SHE_SIM_gal_params/common.hpp"
%include "../SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
%include "../SHE_SIM_gal_params/levels/Galaxy.hpp"
%include "../SHE_SIM_gal_params/levels/Survey.hpp"

// SWIG includes
%include "typemaps.i"
%include "std_vector.i"

namespace std {

%template(GalVector) vector<SHE_SIM::Galaxy *>;

}