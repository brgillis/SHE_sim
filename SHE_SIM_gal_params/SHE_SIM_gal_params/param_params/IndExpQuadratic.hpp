/**********************************************************************\
 @file IndExpQuadratic.hpp
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

#ifndef SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDEXPQUADRATIC_HPP_
#define SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDEXPQUADRATIC_HPP_

#include <initializer_list>

#include <SHE_SIM_gal_params/common.hpp>
#include "SHE_SIM_gal_params/ParamParam.hpp"
#include "SHE_SIM_gal_params/random_functions.hpp"

namespace SHE_SIM
{

/**
 * TODO Auto-generated comment stub
 */
class IndExpQuadratic: public ParamParam
{
private:

	flt_t _N_scale, _hinge_mag, _beta_0, _d_beta, _mag_min, _mag_max;

	// Private methods
	virtual bool is_equal( ParamParam const * const & other ) const override
	{
		IndExpQuadratic const * other_derived = dynamic_cast<IndExpQuadratic const *>(other);
		if(other_derived==nullptr) return false;
		return (_N_scale==other_derived->_N_scale) and (_hinge_mag==other_derived->_hinge_mag) and
				(_beta_0==other_derived->_beta_0) and (_d_beta==other_derived->_d_beta) and
				(_mag_min==other_derived->_mag_min) and (_mag_max==other_derived->_mag_max);
	}

public:

	// Constructor and destructor
	IndExpQuadratic( flt_t const & N_scale = 0., flt_t const & hinge_mag = 1.,
			flt_t const & beta_0 = -1., flt_t const & d_beta = 1.,
			flt_t const & mag_min = -1., flt_t const & mag_max = 1.)
	: ParamParam(ParamParam::INDEPENDENT),
	  _N_scale(N_scale),
	  _hinge_mag(hinge_mag),
	  _beta_0(beta_0),
	  _d_beta(d_beta),
	  _mag_min(mag_min),
	  _mag_max(mag_max)
	{
	}
	virtual ~IndExpQuadratic() {}

	// Get the name of this
	virtual name_t name() const override { return "exp_quadratic"; };

	// PDF generation function
	flt_t get_pdf(flt_t const & mag)
	{
		flt_t m = mag-_hinge_mag;
		flt_t p = _N_scale*std::exp(_beta_0*m + _d_beta*m*m);
		return p;
	}

	// Get the value
	virtual flt_t get_independently( gen_t & gen = rng ) const override
	{
		return rand_from_pdf(get_pdf,1000,_mag_min,_mag_max,gen);
	}

	virtual ParamParam * clone() const override
	{
		return new IndExpQuadratic(*this);
	}

	virtual ParamParam * recreate(const std::vector<flt_t> & params) const override
	{
		if(params.size() != 6) throw std::runtime_error("Invalid number of arguments used for exp_quadratic param param.\n"
				"Exactly 6 arguments are required.");
		return new IndExpQuadratic(params[0],params[1],params[2],params[3],params[4],params[5]);
	}
};

} // namespace SHE_SIM

#endif // SHE_SIM_GAL_PARAMS_PARAM_PARAMS_INDEXPQUADRATIC_HPP_
