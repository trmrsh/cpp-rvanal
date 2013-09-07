#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/rvanal.h"

/**
 * bestsin is a callable function for performing 
 * a chi**2 minimisation over a specified frequency range
 * based upon a  constant + sinusoid model. It comes in two
 * version. 
 * \param data array of times, rv and errors
 * \param flo  lower frequency limt, cycles/unit of t
 * \param fhi  upper frequency limt, cycles/unit of t
 * \param over over-sampling factor to define the number of frequencies
 * to compute. If set to 1, the number of frequencies is such that from one
 * to the next there is a change of 0.1 cycles over the whole baseline of the 
 * input data. Although the routine is faster if this number is reduced,
 * if it is too low, the frequencies become too coarsely sampled and
 * the global minimum may be missed.
 * \param fbest the best (minimum chi**2 frequency.
 * \param cbest the best (minimum) chi**2.
 */

void Rvanal::bestsin(const Subs::Buffer1D<Subs::rv>& data, double flo, double fhi, double over, double &fbest, double &cbest){

  size_t ndat = data.size();
  if(ndat < 4) throw Rvanal_Error("Too few data points in bestsin");

  // Compute number of frequencies needed
  double tmin = 0, tmax = 0;
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

  size_t nfreq = (size_t)(10.*over*(tmax-tmin)*(fhi-flo));
  nfreq = nfreq > 2 ? nfreq : 3;
  
  // Do the search in chunks to avoid grabbing ludicrous amounts of memory
  
  const size_t CHUNK = 100000;
  int nchunk = (nfreq-1)/CHUNK+1;
  size_t ngrab = nfreq/nchunk+1;
  
  double freq[ngrab], chisq[ngrab];
  
  // Construct function object once

  Chisq chi(data);
  
  double f1, f2, fmin, cmin;
  double tol = 1.e-5/(tmax-tmin);

  // Loop through, chunk by chunk with a little overlap
  
  int nmin = 0;
  for(int n=0; n<nchunk; n++){
    
    f1 = flo + 0.99*(fhi-flo)*n/nchunk;
    f2 = flo + (fhi-flo)*(n+1)/nchunk;
    
    // Compute Chi**2/frequency arrays
    
    sinfit_chisq(data,f1,f2,freq,chisq,ngrab);
    if(n == 0){
      fbest = freq[0];
      cbest = chisq[0];
    }

    for(size_t i = 1; i<ngrab-1; i++){
      
      if(chisq[i] < chisq[i-1] && chisq[i] < chisq[i+1]){
	
	// Refine the minimum value
	
	nmin++;
	cmin = Subs::brent(freq[i], freq[i-1], freq[i+1], chi, tol, fmin);	
	if(cmin < cbest){
	  cbest = cmin;
	  fbest = fmin;
	}
      }
    }
  }
}    

/**
 * bestsin is a callable function for performing  a chi**2 minimisation over a specified frequency range
 * based upon a  constant + sinusoid model. This version finds the minimum only over the set of "natural"
 * frequencies. These are N/2 frequencies, where N is the number of valid data points, which are regular multiples
 * of (N-1)/N/T where T is the time-base of the data.
 * \param data array of times, rv and errors
 * \param fbest the best (minimum chi**2 frequency.
 * \param cbest the best (minimum) chi**2.
 */

void Rvanal::bestsin(const Subs::Buffer1D<Subs::rv>& data, double &fbest, double &cbest){

  size_t ndat = data.size();
  if(ndat < 4) throw Rvanal_Error("Too few data points in bestsin");

  // Compute Nyquist frequency (only really for equal spacing)
  
  double tmin, tmax;
  tmin = tmax = data[0].x;
  for(long unsigned int i=1; i<ndat; i++){
    if(data[i].x < tmin) tmin = data[i].x;
    if(data[i].x > tmax) tmax = data[i].x;
  }

  double tunit = (tmax-tmin)/(ndat-1);
  double flo   = 1./tunit/ndat;
  double fhi   = (ndat/2)/tunit/ndat;
  size_t nfreq = ndat/2;
  
  double freq[nfreq], chisq[nfreq];
  
  // Construct function object once

  Chisq chi(data);
  
  sinfit_chisq(data,flo,fhi,freq,chisq,nfreq);
  fbest = freq[0];
  cbest = chisq[0];

  for(long unsigned int i=1; i<nfreq; i++){
    if(chisq[i] < cbest){
      cbest = chisq[i];
      fbest = freq[i];
    }
  }
}    

// calculates over natural frequency range but at a 
// controllable number of frequencies.

/**
 * bestsin is a callable function for performing  a chi**2 minimisation over a specified frequency range
 * based upon a  constant + sinusoid model. This version finds the minimum only over the set of "natural"
 * frequencies but with a defined over-sampling factor allowed on top of this. The natural frequencies are N/2 frequencies, 
 * where N is the number of valid data points, which are regular multiples of (N-1)/N/T where T is the time-base of the data.
 * \param data array of times, rv and errors
 * \param over over-sampling factor to define the number of frequencies
 * to compute. If set to 1, the number of frequencies is such that from one
 * to the next there is a change of 0.1 cycles over the whole baseline of the 
 * input data. Although the routine is faster if this number is reduced,
 * if it is too low, the frequencies become too coarsely sampled and
 * the global minimum may be missed.
 * \param fbest the best (minimum chi**2 frequency.
 * \param cbest the best (minimum) chi**2.
 */

void Rvanal::bestsin(const Subs::Buffer1D<Subs::rv>& data, double over, double &fbest, double &cbest){

  size_t ndat = data.size();
  if(ndat < 4) throw Rvanal_Error("Too few data points in bestsin");

  // Compute Nyquist frequency (only really for equal spacing)
  
  double tmin, tmax;
  tmin = tmax = data[0].x;
  for(long unsigned int i=1; i<ndat; i++){
    if(data[i].x < tmin) tmin = data[i].x;
    if(data[i].x > tmax) tmax = data[i].x;
  }

  double tunit = (tmax-tmin)/(ndat-1);
  double flo   = 1./tunit/ndat;
  double fhi   = (ndat/2)/tunit/ndat;
  size_t nfreq = size_t(floor(over*ndat/2+0.5));
  
  double freq[nfreq], chisq[nfreq];
  
  // Construct function object once

  Chisq chi(data);
  
  sinfit_chisq(data,flo,fhi,freq,chisq,nfreq);
  fbest = freq[0];
  cbest = chisq[0];

  for(long unsigned int i=1; i<nfreq; i++){
    if(chisq[i] < cbest){
      cbest = chisq[i];
      fbest = freq[i];
    }
  }
}    












