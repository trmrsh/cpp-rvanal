#include <math.h>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_rvanal.h"

/**
 * m2min computes the minimum mass of the companion to 
 * a star of mass \c m1, observed to undergo rv changes of 
 * amplitude \c amp km/s at frequency \c freq cyc/day.
 * \param freq frequency in cycles/day
 * \param amp  semi-amplitude in km/s.
 * \param m1   primary mass (solar masses)
 */
float Rvanal::m2min(double freq, double amp, double m1){
  double f  = fm(freq,amp);
  double m2 = pow(Subs::sqr(m1)*f,0.333333);
  double mold = 0.;
  while(fabs(mold-m2) > 1.e-8*m2){
    mold = m2;
    m2   = pow(f*Subs::sqr(m1+m2),0.3333333333333333);
  }
  return m2;
}
