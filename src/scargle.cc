//
// Direct method for computing Lomb-Scargle periodogram. Advantage
// over other Press-Rybicki is that you do not have to compute
// whole range. This is worth it if there are very irregular data
// gaps. Uses trigonometric recurrence to compute cos and sin.
//

#include <stdlib.h>
#include <math.h>
#include "trm/subs.h"
#include "trm/rvanal.h"

void Rvanal::scargle(double x[], float y[], float e[], 
		     size_t n, double f1, double f2, 
		     double freq[], float pgram[], 
		     size_t nfreq, int &status){

  if(status) return;

  const double TWOPI = 6.2831853071795865;
  double xmin, xmax, xmid;

  if(status) return;

  for(i=0;(i<n && e[i]<=0.);i++);
  if(e[i]<=0.){
    printf("No valid values in input data\n");
    status = 1;
    return;
  }
  xmin = xmax = x[i];
  for(i=0;i<n;i++){
    if(x[i] < xmin) xmin = x[i];
    if(x[i] > xmax) xmax = x[i];
  }
  xmid = 0.5*(xmin+xmax);

  // Initialise recurrences

  double theta, *alpha, *beta, *czero, *szero, t, df;

  alpha = new double [n];
  beta  = new double [n];
  czero = new double [n];
  szero = new double [n];
  df  = (f2-f1)/(nfreq-1);
  for(size_t j=0;j<n;j++){
    if(e[j]>0.){
	theta    = df*(t = TWOPI*(x[j]-xmid));
	alpha[j] = -2.*Subs::sqr(sin(0.5*theta));
	beta[j]  = sin(theta);
	theta    = f1*t;
	czero[j] = cos(theta);
	szero[j] = sin(theta);
      }
  }

  // OK, off we go ...

  float cc, cs, ss, cy, sy, wgt, c, s;
  for(size_t i=0;i<nfreq;i++){
    freq[i] = f1+df*i;
    cc=cs=ss=cy=sy=0.;
    for(size_t j=0;j<n;j++){
      if(e[j]>0.){
	wgt = 1./Subs::sqr(e[j]);
	c = czero[j];
	s = szero[j];
        cc += wgt*c*c;
        cs += wgt*c*s;
	ss += wgt*s*s;
        cy += wgt*c*y[j];
        sy += wgt*s*y[j];

	// Update cos and sin. Brackets needed here to reduce roundoff
	// as far as possible.

	t        = czero[j];
	czero[j] = czero[j] + (alpha[j]*czero[j]-beta[j]*szero[j]);
	szero[j] = szero[j] + (alpha[j]*szero[j]+beta[j]*t);
      }
    }
    
    // Evaluate periodogram value

    pgram[i] = (ss*cy*cy-2.*cs*sy*cy+cc*sy*sy)/(cc*ss-cs*cs)/2.;
  }
  delete[] alpha;
  delete[] beta;
  delete[] czero;
  delete[] szero;
  return;
}



