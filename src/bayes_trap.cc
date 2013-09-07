#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/rvanal.h"

//! Function object to support Bayesian integration

/** This function object allows one to provide a f(x) type interface to 
 * Rvanal::bayes_prob(const Subs::Buffer1D<Subs::rv>&, double, double, double, int, long int&)
 * as needed by 'trapzd'
 */

class Bayes_Func {
public:

  //! Constructor to store necessary data
  Bayes_Func(const Subs::Buffer1D<Subs::rv>& data, double gmax, double mmax, int nmonte, Subs::INT4 seed) : 
    data_(data), gmax_(gmax), mmax_(mmax), nmonte_(nmonte), seed_(seed) {}

  //! Function operator that just needs the frequency
  double operator()(double f){
    double kmax = pow(Constants::TWOPI*Constants::G*Constants::MSUN*mmax_*f/Constants::DAY,1./3.)/1000.;
    return Rvanal::bayes_prob(data_, f, gmax_, kmax, nmonte_, seed_);
  }
private:
  const Subs::Buffer1D<Subs::rv>& data_;
  double gmax_, mmax_; 
  int nmonte_;
  Subs::INT4 seed_;
};

/**
 * Given a set of data that can be fitted with a constant plus a sinusoid,
 * this function carries out trapezoidal integration of the probability over a specified
 * frequency range. Call it successively to get finer and finer sampling of the integrand
 * (which can be VERY peaky).
 *
 *  \sa double Rvanal::bayes_prob(const Subs::Buffer1D<Subs::rv>&, double, double, double, int, long int)
 * for gory details of what is going one here.
 *
 * \param data   the set of times, velocities and errors
 * \param flo    the lower frequency of interest in cycles per unit time
 * \param fhi    the lower frequency of interest in cycles per unit time
 * \param gmax   the maximum gamma velocity (establishes range of integration)
 * \param mmax   the maximum companion mass, solar masses if data has units of days amnd km/s
 * \param nmonte the number of samples to take, 0 will just attempt the analytic result only.
 * \param seed   the seed integer for the monte carlo integration for nmonte > 0
 * \param first  true for the first call of a sequence, false for subsequent ones (the program keeps track internally of
 * the number of calls but you must not alter any other argument during such a sequence).
 * \return The function return the log10 of the probability integral. This only has meaning in a relative
 * sense when compared against another value for a different frequency range.
 * \exception The routine may throw a Rvanal:Rvanal_Error
 */

double Rvanal::bayes_trap(const Subs::Buffer1D<Subs::rv>& data, double flo, double fhi, double gmax, double mmax, int nmonte, Subs::INT4 seed, bool first){
  
  // The routine takes care to avoid overflow by extracting whatever value is largest from the logarithm at all times
  // This is called cfactor. The remaining part is called sum.
  static int ncall;
  if(first)
    ncall = 1;
  else
    ncall++;
  static double cfactor, sum;

  // Create function object
  Bayes_Func func(data, gmax, mmax, nmonte, seed);

  if(ncall == 1){

    // First point just evaluate at each end of the range
    double fval1 = func(flo);
    double fval2 = func(fhi);

    cfactor = std::max(flo, fhi);
    sum     = 0.5*(fhi-flo)*(exp(fval1-cfactor)+exp(fval2-cfactor));
  }else{

    // Subsequent calls add points halway between previous points.

    // Compute number of points to add.
    long int npts;
    int j;
    for(npts=1, j=1; j<ncall-1; j++) npts <<= 1;

    // Add new points, changing exponential factor as necessary
    long int n;
    double f, fval, nsum = 0.;
    for(n=0; n<npts; n++){
      f = flo + (fhi-flo)*(n+0.5)/double(npts);
      fval = func(f);
      if(fval > cfactor){
	double scale = exp(cfactor-fval);
	sum    *= scale;
	nsum   *= scale;
	cfactor = fval;
	nsum   += 1.;
      }else{
	nsum   += exp(fval-cfactor);
      }
    }
    sum = 0.5*(sum+(fhi-flo)*nsum/double(npts));
  }
  // Return log10 of result
  return (cfactor + log(sum))/log(10.);
}










