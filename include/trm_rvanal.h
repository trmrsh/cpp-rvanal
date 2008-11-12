#ifndef TRM_RVANAL
#define TRM_RVANAL

#include <string>
#include "trm_subs.h"

//! Rvanal namespace

/**
 * Namespace for the rvanal suite of programs and functions.
 */

namespace Rvanal {

  //! Default directory for command defaults
  const char RVANAL_DIR[] = ".rvanal";

  //! Environment variable for switching directory for command defaults
  const char RVANAL_ENV[] = "RVANAL_ENV";

  //! An exception class.

  /** Rvanal::Rvanal_Error is the error class for the Rvanal programs.
   * It is inherited from the standard string class.
   */
  class Rvanal_Error : public std::string {
  public:
    //! Default constructor
    Rvanal_Error() : std::string() {}
    //! Constructor from a string (e.g. an error message).
    Rvanal_Error(const std::string& err) : std::string(err) {} 
  };

  //! Structure returning covariances on a sinfit from sinfit_errors
  /*
   * This structure stores the co-variances on the parameters \c c, \c a, \c p and \c t0 
   * of the fit "c+a*sin(TWOPI*(t-t0)/p)
   * \var vcc variance of constant term var(c,c)
   * \var vca covariance of constant and semi-amplitude var(c,a)
   * \var vcp covariance of constant and period var(c,p)
   * \var vca covariance of constant and tzero var(c,t0)
   * \var vaa variance of semi-amplitude var(a,a)
   * \var vap covariance of semi-amplitude and period var(a,p)
   * \var vap covariance of semi-amplitude and tzero var(a,t0)
   * \var vpp variance of period var(p,p)
   * \var vpt covariance of period and tzero var(p,t0)
   * \var vtt variance of tzero var(t0,t0)
   */
  struct sinerr{
    double vcc, vca, vcp, vct, vaa, vap, vat, vpp, vpt, vtt;
  };

  //! Calculates probability of a given frequency
  double bayes_prob(const Subs::Buffer1D<Subs::rv>& data, double f, double gmax, double kmax, int nmonte, Subs::INT4& seed);
  
  //! Integrates probability over a range of frequency
  double bayes_trap(const Subs::Buffer1D<Subs::rv>& data, double flo, double fhi, double gmax, double kmax, int nmonte, Subs::INT4 seed, bool first);

  //! Calculates optimised best sine fit 
  void bestsin(const Subs::Buffer1D<Subs::rv>& data, double flo, double fhi, 
	       double over, double &fbest, double &cbest);

  //! Calculates best sine fit over 'natural' frequencies
  void bestsin(const Subs::Buffer1D<Subs::rv>& data, double &fbest, double &cbest);

  //! Calculates best sine fit over oversampled 'natural' frequencies
  void bestsin(const Subs::Buffer1D<Subs::rv>& data, double over, double &fbest, 
	       double &cbest);
  
  float dpb(const Subs::Buffer1D<Subs::rv>& data, double period, double phase, float k, size_t ntrial, double thresh, Subs::INT4& seed);

  //! Computes mass function
  float fm(double freq, double amp);
  
  //! Computes minimum mass of secondary star
  float m2min(double freq, double amp, double m1);
  
  //! Computes Chi**2 over an array of frequencies
  void sinfit_chisq(const Subs::Buffer1D<Subs::rv>& data, double f1, double f2, 
		    double freq[], double chisq[], size_t nfreq);
  
  //! Computes Chi**2 at a single frequency
  double sinfit_chisq(const Subs::Buffer1D<Subs::rv>& data, double f, 
		      float& con, float& amp, double& t0);
  
  //! Computes variances on best fit
  sinerr sinfit_errors(const Subs::Buffer1D<Subs::rv>& data, double freq); 
  
  //! Computes Scargle periodogram
  void scargle(double x[], float y[], float e[], size_t n, 
	       double f1, double f2, double freq[], float pgram[], 
	       size_t nfreq, int &status);

  //! Function object for call by brent
  
  /**
   * Stores data array so that an f(x) type call can return
   * the Chi**2 for a given frequency to allow minimisation
   */
  class Chisq : public Subs::Sfunc {
    const Subs::Buffer1D<Subs::rv>& dat;
  public:
    Chisq(const Subs::Buffer1D<Subs::rv>& data) : dat(data) {}
    
    double operator()(double f) {
      float con, amp;
      double t0;
      return sinfit_chisq(dat, f, con, amp, t0);
    }
  };
}

#endif






