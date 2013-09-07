#include <cmath>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/rvanal.h"

/**
 * Given a semi-amplitude, phase and period, \c dpb
 * computes the probability of detection given a particular
 * set of times, and radial velocity uncertainties. It does
 * this using a random versus constant model chi**2 difference.
 * Since the random model gives chi**2 = 0, this reduces to
 * using the chi**2 after fitting a constant. 
 * \param data the data
 * \param period trial period
 * \param phase trial phase
 * \param k orbital velocity (unprojected); the program automatically averages
 * over randomly inclined orbits
 * \param ntrial if ntrial > 0 then the program will perform an averaging
 * over noise by adding gaussian uncertainties following the input errors
 * (see below). It does this ntrial times. The effect is usually small.
 * \param thresh threshold chi**2 to count as a detection
 * \param seed random number seed.
 */

float Rvanal::dpb(const Subs::Buffer1D<Subs::rv>& data, double period, double phase, float k, size_t ntrial, double thresh, Subs::INT4& seed){

  // Generate rvs for edge-on case. We will then use
  // these to generate a detection probability automatically
  // corrected for random orientation.
  //
  // We need to compute the coefficients C(v,v), C(v,n)
  // and C(n,n) where C stands for co-variance. The C(v,v)
  // can be computed outside the noise averaging loop to save
  // time.

  int l, ndat = 0;

  Subs::Buffer1D<float> v(data.size()), w(data.size());
  double mean = 0., sumw = 0.;

  for(l=0; l<data.size(); l++){
    if(data[l].z > 0.){
      ndat++;
      w[l] = 1./Subs::sqr(data[l].z);
      v[l]  = k*sin(Constants::TWOPI*(data[l].x/period-phase));
      sumw += w[l];
      mean += w[l]*v[l];
    }
  }
  mean /= sumw;
  double a = 0.;
  for(l=0; l<data.size(); l++) 
    if(data[l].z > 0.) 
      a += w[l]*Subs::sqr(v[l]-mean);
  
  double det;
  if(ntrial){

    // yes, we will average over noise
    
    Subs::Buffer1D<float> n(ndat);
    det = 0.;
    for(size_t nt=0; nt < ntrial; nt++){ 
      
      // Generate gaussian noise
      
      double sum = 0.;
      for(l=0;l<ndat;l++)
	if(data[l].z > 0.) 
	  sum += w[l]*(n[l] = data[l].z*Subs::gauss2(seed));

      sum /= sumw;
      double temp, b = 0., c = -thresh;
      for(l=0; l<ndat; l++){
	if(data[l].z > 0.){
	  b += w[l]*(v[l]-mean)*(temp = n[l]-sum);
	  c += w[l]*Subs::sqr(temp);
	}
      }

      // first treat a=0 case which although it should
      // be very rare is not impossible

      if(a == 0.0){

	// linear equation of the form 2*b*sini+c = 0
	// need to find over what part of range 0 to 1
	// this is > 0 ==> detection.
	
	if(b == 0.0){
	  if(c > 0.0) det++;
	}else{
	  double d = -c/(2.0*b);
	  if(d < 0.0){
	    if(b > 0.0) det++;
	  }else if(d > 1.0){
	    if(b > 0.0) det++;
	  }else{
	    if(b > 0.0){
	      det += sqrt(1.0-d*d);
	    }else{
	      det += 1.-sqrt(1.0-d*d);
	    }
	  }
	}
      }else{
	
	// quadratic equation of the form a*sini*sini+2*b*sini+c = 0
	// need to find over what part of range 0 to 1
	// this is > 0 ==> detection. There are 7 different cases
	// to sort out. Case 7 (s1 < 0, s2 > 1) does not appear
	// as it leads to no detections.
	
	double d = b*b-a*c;
	if(d < 0.){        
	  
	  // case 1: quadratic > 0 for all sini
	  
	  det++; 
	  
	}else{
	  
	  // evaluate quadratic roots with formula from
	  // NR to avoid roundoff
	  
	  double q = -(b+Subs::sign(sqrt(d),b));
	  double s1 = q/a;
	  double s2 = c/q;
	  if(s1 > s2) std::swap(s1, s2);
	  
	  // now go through 5 of remaining 6 cases, the sixth
	  // leading to no detections.
	  
	  if(s1 >= 1.0 || s2 <= 0.0){ 
	    
	    // cases 2 and 3: quadratic > 0 for sini = 0 to 1
	    
	    det++;
	    
	  }else if(s1 > 0.0 && s2 >= 1.0){
	    
	    // case 4: quadratic > 0 for sini = 0 to s1
	    
	    det += 1.0-sqrt(1.0-s1*s1); 
	    
	  }else if(s1 <= 0.0 && s2 < 1.0){
	    
	    // case 5: quadratic > 0 for sini = s2 to 1
	    
	    det += sqrt(1.0-s2*s2);
	    
	  }else if(s1 > 0.0 && s2 < 1.0){
	    
	    // case 6: quadratic > 0 for sini = 0 to s1 and s2 to 1
	    
	    det += 1.0-sqrt(1.0-s1*s1)+sqrt(1.0-s2*s2);
	    
	  }
	}
      }
    }
    det /= ntrial;
  }else{
    
    // no attempt to allow for noise, just compare average
    // expected value of chi**2 (= a*sini*sini+ndat-1)
    // with thresh to work out detection probability.
    
    double sinsqi;
    //    if(a > 0.0 && (sinsqi = (thresh-ndat+1)/a) < 1){
    if(a > 0.0 && (sinsqi = (thresh)/a) < 1){
      det = sqrt(1.-sinsqi);
    }else{
      det = 0.;
    }
  }

  return det;
}






