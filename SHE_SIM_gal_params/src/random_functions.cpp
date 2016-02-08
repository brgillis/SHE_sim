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
#include <functional>
#include <numeric>
#include <stdexcept>

#include <Eigen/Core>

#include <SHE_SIM_gal_params/common.hpp>
#include <SHE_SIM_gal_params/random_functions.hpp>

namespace SHE_SIM {

// Implement the random number generator engine
gen_t rng;

#include <functional>
typedef Eigen::Array<flt_t,Eigen::Dynamic,1> flt_array_t;

flt_t rand_from_cdf_arrays( flt_array_t const & xvals, flt_array_t cvals, gen_t & gen = rng )
{
	assert(xvals.size()>=2 and xvals.size()==cvals.size());

	// Quietly normalize the cvals
	cvals -= cvals[0];

	flt_t cmax = cvals[cvals.size()-1];

	if(cmax<=0)
	{
		throw std::runtime_error("Invalid values used for generating random value: Final CDF value is <= 0.");
	}

	cvals /= cmax;

	// Generate a random value
	flt_t const r = drand(0.,1.,rng);

	// Get the index on the cdf where this lies
	flt_array_t diffs = (cvals-r).abs();
	int_t i,j;
	diffs.minCoeff(&i,&j);

	// If the value at the index is below r, or the index is zero, move up one index
	while(((cvals[i]<r) and (i<cvals.size())) or (i==0))
	{
		++i;
	}

	// Interpolate to estimate the value
	flt_t const clow = cvals[i - 1];
	flt_t const chi = cvals[i];
	flt_t const xlow = xvals[i - 1];
	flt_t const xhi = xvals[i];

	flt_t res = xlow + (xhi - xlow) / (chi - clow) * (r - clow);

	// Check for edge cases
	if( res < xvals[0] ) res = xvals[0];
	if( res > xvals[xvals.size()-1] ) res = xvals[xvals.size()-1];

	return res;

}

flt_t rand_from_cdf( std::function<flt_t(flt_t const &)> const & f, int_t const & N_samples,
		flt_t const & xlow, flt_t const & xhigh, gen_t & gen )
{
	// Get an array of x points and cdf values at those points
	flt_array_t xvals = flt_array_t::LinSpaced(N_samples, xlow, xhigh);
	flt_array_t cvals = xvals.unaryExpr(f);

	flt_t res = rand_from_cdf_arrays( xvals, cvals, gen );

	return res;
}

flt_t rand_from_pdf( std::function<flt_t(flt_t const &)> const & f, int_t const & N_samples,
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
