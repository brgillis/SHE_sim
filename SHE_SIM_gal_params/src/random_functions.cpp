/**********************************************************************\
 @file random_functions.cpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015  Bryan R. Gillis

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

#include <cassert>
#include <Eigen/Core>
#include <numeric>

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/random_functions.hpp>

namespace SHE_SIM {

// Implement the random number generator engine
gen_t rng;

typedef Eigen::Array<flt_t,Eigen::Dynamic,1> flt_array_t;

flt_t rand_from_cdf_arrays( flt_array_t const & xvals, flt_array_t cvals, gen_t & gen = rng )
{
	assert(xvals.size()>=2 and xvals.size()==cvals.size());

	// Quietly normalize the cvals
	cvals -= cvals[0];
	cvals /= cvals[-1];

	// Generate a random value
	flt_t const r = drand(0.,1.,rng);

	// Get the index on the cdf where this lies
	int_t i = (cvals-r).abs().minCoeff();

	// If the value at the index is below r, move up one index
	if((cvals[i]>r) and (i<cvals.size())) ++i;

	// Due to the way random generation works, we can safely ignore the pathological edge cases here
	// Interpolate to estimate the value
	flt_t const clow = cvals[i - 1];
	flt_t const chi = cvals[i];
	flt_t const xlow = xvals[i - 1];
	flt_t const xhi = xvals[i];

	flt_t const res = xlow + (xhi - xlow) / (chi - clow) * (r - clow);

	return res;

}

flt_t rand_from_cdf( std::function<flt_t(flt_t)> const & f, int_t const & N_samples,
		flt_t const & xlow, flt_t const & xhigh, gen_t & gen )
{
	// Get an array of x points and cdf values at those points
	flt_array_t xvals = flt_array_t::LinSpaced(N_samples, xlow, xhigh);
	flt_array_t cvals = xvals.unaryExpr(f);

	flt_t res = rand_from_cdf_arrays( xvals, cvals, gen );

	return res;
}

flt_t rand_from_pdf( std::function<flt_t(flt_t)> const & f, int_t const & N_samples,
		flt_t const & xlow, flt_t const & xhigh, gen_t & gen )
{
	// Get an array of x points and pdf values at those points
	flt_array_t xvals = flt_array_t::LinSpaced(N_samples, xlow, xhigh);
	flt_array_t pvals = xvals.unaryExpr(f);

	// Get (unnormalized) cdf values
	flt_array_t cvals(pvals.size());

	std::partial_sum(pvals.data(), pvals.data()+pvals.size(), cvals.data(), std::plus<flt_t>());

	flt_t res = rand_from_cdf_arrays( xvals, cvals, gen );

	return res;
}

} // namespace SHE_SIM
