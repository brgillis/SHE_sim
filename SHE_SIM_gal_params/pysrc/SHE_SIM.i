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

// SWIG includes
%include "typemaps.i"
%include "std_string.i"
%include "std_vector.i"

%module SHE_SIM

%{
	 
	/* Includes the header in the wrapper code */
	#include "SHE_SIM_gal_params/common.hpp"
	#include "SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
	#include "SHE_SIM_gal_params/levels/Cluster.hpp"
	#include "SHE_SIM_gal_params/levels/ClusterGroup.hpp"
	#include "SHE_SIM_gal_params/levels/Field.hpp"
	#include "SHE_SIM_gal_params/levels/FieldGroup.hpp"
	#include "SHE_SIM_gal_params/levels/Galaxy.hpp"
	#include "SHE_SIM_gal_params/levels/GalaxyGroup.hpp"
	#include "SHE_SIM_gal_params/levels/Image.hpp"
	#include "SHE_SIM_gal_params/levels/ImageGroup.hpp"
	#include "SHE_SIM_gal_params/levels/Survey.hpp"
	 
%}
 
// Parse the header files to generate wrappers
%include "../SHE_SIM_gal_params/common.hpp"
%include "../SHE_SIM_gal_params/ParamHierarchyLevel.hpp"
%include "../SHE_SIM_gal_params/levels/Cluster.hpp"
%include "../SHE_SIM_gal_params/levels/ClusterGroup.hpp"
%include "../SHE_SIM_gal_params/levels/Field.hpp"
%include "../SHE_SIM_gal_params/levels/FieldGroup.hpp"
%include "../SHE_SIM_gal_params/levels/Galaxy.hpp"
%include "../SHE_SIM_gal_params/levels/GalaxyGroup.hpp"
%include "../SHE_SIM_gal_params/levels/Image.hpp"
%include "../SHE_SIM_gal_params/levels/ImageGroup.hpp"
%include "../SHE_SIM_gal_params/levels/Survey.hpp"

%apply int& {SHE_SIM::level_t&};

namespace SHE_SIM {

%extend ParamHierarchyLevel { 
	void set_param_params(name_t name, name_t param_type )
	{
		return $self->set_param_params( name, param_type );
	}
	void set_param_params(name_t name, name_t param_type, flt_t arg1 )
	{
		return $self->set_param_params( name, param_type, arg1 );
	}
	void set_param_params(name_t name, name_t param_type, flt_t arg1, flt_t arg2 )
	{
		return $self->set_param_params( name, param_type, arg1, arg2 );
	}
	void set_param_params(name_t name, name_t param_type, flt_t arg1, flt_t arg2, flt_t arg3 )
	{
		return $self->set_param_params( name, param_type, arg1, arg2, arg3 );
	}
	void set_param_params(name_t name, name_t param_type, flt_t arg1, flt_t arg2, flt_t arg3, flt_t arg4 )
	{
		return $self->set_param_params( name, param_type, arg1, arg2, arg3, arg4 );
	}
	void set_param_params(name_t name, name_t param_type, flt_t arg1, flt_t arg2, flt_t arg3, flt_t arg4, flt_t arg5 )
	{
		return $self->set_param_params( name, param_type, arg1, arg2, arg3, arg4, arg5 );
	}
	void set_param_params(name_t name, name_t param_type, flt_t arg1, flt_t arg2, flt_t arg3, flt_t arg4, flt_t arg5, flt_t arg6 )
	{
		return $self->set_param_params( name, param_type, arg1, arg2, arg3, arg4, arg5, arg6 );
	}
}

}

// Tell SWIG to implement vectors of the PHL types as lists
namespace std {

%template(PHLVector) vector<SHE_SIM::ParamHierarchyLevel *>;

%template(ClusterVector) vector<SHE_SIM::Cluster *>;
%template(ClusterGroupVector) vector<SHE_SIM::ClusterGroup *>;
%template(FieldVector) vector<SHE_SIM::Field *>;
%template(FieldGroupVector) vector<SHE_SIM::FieldGroup *>;
%template(GalaxyVector) vector<SHE_SIM::Galaxy *>;
%template(GalaxyGroupVector) vector<SHE_SIM::GalaxyGroup *>;
%template(ImageVector) vector<SHE_SIM::Image *>;
%template(ImageGroupVector) vector<SHE_SIM::ImageGroup *>;

}