#include <cmath>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/buffer2d.h"
#include "trm/rvanal.h"

/**
 * Computes the chi**2 after fitting a sine wave plus
 * a constant to some data. It does this over an array of equally spaced
 * frequencies (some time svaings can be made versus calling a single frequency
 * routine many times). It is very similar to the 
 * Lomb-Scargle method except that it accounts properly for the constant term.
 * \param data times, velocities and uncertainties
 * \param f1 the first frequency in cycles/(unit of x).
 * \param f2 the last  frequency in cycles/(unit of x).
 * \param freq   array of frequencies (returned)
 * \param chisq array of Chi**2 (returned)
 * \param nfreq The number of frequencies
 */

void Rvanal::sinfit_chisq(const Subs::Buffer1D<Subs::rv>& data, double f1, double f2, double freq[], double chisq[], size_t nfreq){

  size_t i, j, n = data.size();

  double tmin = 0, tmax = 0, tmid;
  for(int i=0; i<data.size(); i++){
    if(data[i].z > 0.) {
      tmin = tmax = data[i].x;
      break;
    }
  }
  for(int i=1; i<data.size(); i++){
    if(data[i].z > 0.){
      if(data[i].x < tmin) tmin = data[i].x;
      if(data[i].x > tmax) tmax = data[i].x;
    }
  }   
  tmid = 0.5*(tmin+tmax);

  // Initialise recurrences
  double delta, theta, t, df;
  Subs::Buffer1D<double> alpha(n), beta(n), czero(n), szero(n);

  df    = (f2-f1)/(nfreq-1);
  for(j=0; j<n; j++){
    if(data[j].z > 0.){
      delta     = df*(t = Constants::TWOPI*(data[j].x-tmid));
      alpha[j]  = -2.*Subs::sqr(sin(0.5*delta));
      beta[j]   = sin(delta);
      theta     = f1*t;
      czero[j]  = cos(theta);
      szero[j]  = sin(theta);
    }
  }

  // Workspace arrays
  Subs::Buffer1D<float>  a(3), w(3);
  Subs::Buffer2D<float>  u(3,n), v(3,3);

  for(i=0; i<nfreq; i++){
    freq[i]  = f1+df*i;
    chisq[i] = Subs::svdfit(data,a,czero,szero,u,v,w);

    // update cos and sin arrays (NR, page 178/179)
    for(j=0;j<n;j++){
      if(data[j].z > 0.){
	t         = czero[j];
	czero[j] += (alpha[j]*t-beta[j]*szero[j]);
	szero[j] += (alpha[j]*szero[j]+beta[j]*t);
      }
    }
  }
}

/**
 * Computes the chi**2 after fitting a sine wave plus
 * a constant to some data. It also returns the best fit parameters.
 * Uses singular value decomposition.
 * \param data times, velocities and uncertainties
 * \param f the frequency in cycles/(unit of x).
 * \param con the constant term of the best fit
 * \param amp the semi-amplitude of the best fit
 * \param t0 the positive-going zero crossing time of the best fit
 * \return The chi**2 of the best fit at the particular frequency is returned
 */

double Rvanal::sinfit_chisq(const Subs::Buffer1D<Subs::rv>& data, double f, float& con, float& amp, double& t0){

  size_t i, j, n;
  double tmin, tmax, tmid;

  n = data.size();

  tmin = tmax = data[0].x;
  for(i=1;i<n;i++){
    if(data[i].z > 0.){
      if(data[i].x < tmin) tmin = data[i].x;
      if(data[i].x > tmax) tmax = data[i].x;
    }
  }
  tmid = 0.5*(tmin+tmax);
  
  // Workspace arrays
  Subs::Buffer1D<float>   a(3), w(3);
  Subs::Buffer1D<double>  cosine(n), sine(n);
  Subs::Buffer2D<float>   u(3,n), v(3,3);

  // Compute cos and sin  
  double theta;
  for(j=0;j<n;j++){
    if(data[j].z > 0.){
      theta     = f*Constants::TWOPI*(data[j].x-tmid);
      cosine[j] = cos(theta);
      sine[j]   = sin(theta);
    }
  }
  double chisq = Subs::svdfit(data,a,cosine,sine,u,v,w); 
  con = a[0];
  amp = sqrt(Subs::sqr(a[1])+Subs::sqr(a[2]));
  t0  = tmid + atan2(-a[1],a[2])/Constants::TWOPI/f;
  return chisq;
}









