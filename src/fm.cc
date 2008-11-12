#include "trm_constants.h"
#include "trm_rvanal.h"

/**
 * fm computes the mass function equivalent to a
 * particular fequency and amplitude. The value is computed in
 * solar masses.
 * \param freq frequency in cycles/day
 * \param amp  semi-amplitude in km/s.
 */

float Rvanal::fm(double freq, double amp){
  return (Constants::DAY*pow(amp,3)/freq/(Constants::TWOPI*Constants::G*Constants::MSUN/1.e9));
}
