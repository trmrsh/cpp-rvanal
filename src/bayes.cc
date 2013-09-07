/*

!!begin 
!!title  Bayesian integration of period probnability
!!author T.R. Marsh
!!created 7 August 2002
!!revised 10 Aug 2005
!!index  bayes
!!root   bayes 
!!descr  bayes integrates the probability over a range of frequency
!!class  Program
!!css    style.css
!!head1  bayes integrates the probability over a range of frequency

!!emph{bayes} integrates the posterior probability of a period of a
sinusoid in a constant plus sinusoid model fit given the data.  It is
designed to be used to answer the question "what is the probability
that the period lies in the range a to b?"

To do this one first needs to integrate out the "nuisance" parameters
which are the constant and the cosine and sine amplitudes (c, a and b
say) at each period. One needs to define prior probabilities for
these. This is a difficult question which !!emph{bayes} ducks rather in just
assuming a uniform distribution over the volume of a cylinder in the
3D parameter space defined by c, a and b. This cylinder has its axis
along the 'c' axis (i.e. the constant or equivalently systemic
velocity for radial velocity applications) and extends from -gmax to
+gmax. It has radius kmax, i.e. the maximum likely semi-amplitude,
although in this case this is actually set by specification of the
maximum likely companion mass.  While this may well not be accurate,
it is hard to justify a different version, always a problem with
"priors" in Bayesian anlaysis. Thus it seems best to assume that the
answers given by 'bayes' are uncertain by some factor depending upon
the uncertainty in the prior. This probably means that one should be
completely confident in a period unless the chance of its !!emph{not} being
correct is pretty small, say 1 in 1000 or less. There is a similar
uncertainty in the prior over frequency which only confirms this
point. The uniform prior has the huge advantage that the nuisance
parameter integration can almost be done analytically.

The integration proceeds by trapezoidal integration over frequency,
each point of which may involve a Monte Carlo integral for the
nuisance parameter integration. The latter requirement is because the
integrand, which is a 3D gaussian, can extend outside the cylindrical
integration region, especially if there is degeneracy amongst the
parameters, and the Monte Carlo part corrects for this. With good data
however this does not usually make much difference because the 3D gaussian
is pretty much contained within the cylinder in most cases.  The integrand
over frequency is dominated by the factor exp(-Chi**2/2) and can therefore
be !!emph{extremely} peaky and so a fine sampling may be required, which
can take time. The result returned is the log to the base 10 of the
answer, but only makes sense in relative terms, i.e. when compared
against the integral over a different range.

!!head2 Arguments

!!table 
!!arg{ data  }  { Data file. Times must be in units of days, velocities in km/s for the
maximum mass function parameter to work properly.}
!!arg{ flo   }  { Lower frequency limt of integration, cycles/day}
!!arg{ fhi   }  { Upper frequency limt of integration, cycles/day}
!!arg{ gmax  }  { The half range of the constant factor. i.e. this is the maximum you believe
this should ever reach, perhaps 500 km/s for RV work (~ escape velocity from Galaxy)}
!!arg{ mmax  }  { Maximum mass function. This is how the maximum semi-amplitude is set,
being varied according to the orbital period. Basically this should be the maximum companion
mass you would believe in. It is measured in solar masses if the units of the data file are
days and km/s.}
!!arg{sigma}{Systematic uncertainty to add in quadrature to uncertainties}
!!arg{ncall}{Number of calls to trapezoidal integration (each one double the time taken, but you need
enough to sample the integrand properly, so the absolute minimum is approx log(4.*(fhi-flo)*T)/log(2)+1,
i.e. enough to make the frequency spacing fine enough that the frequency step represents a 1/4 cycle shift
over the baseline T of the data. For data taken over a couple of years, with a 20 cycle/day frequency
range, this is of order 16 (implying > 65,000 points), but this is a !!emph{minimum}, and more like 18 to 20
appear to be needed in practice. Every value is printed so that convergence can be judged. A value of 0 will 
go on until Ctrl-C is hit. ncall means the minimum number of calls if the convergence parameter 'change' > 0}
!!arg{change}{The program can be terminated when the value does not change by more than 'change'. If 'change' <=0
it is ignored. If it is > 0, then 'ncall' represents the minimum number of calls to make in order to prevent
spurious early convergence. Remember that 'change' represents the change in the log10 of the probability and
so useful values are << 1.}
!!arg{nmonte}{The number of samples to use per frequency point when carrying out the Monte Carlo correction
for degeneracy amongst the three nuisance parameters. If set to zero, no correction is made, but the
analytic result is taken directly. This is equivalent to assuming that probability is only significant
within the volume of integration and in general is an over-estimate, although for reasonable data it
might be very close to the mark, and it will be fastest.}
!!arg{seed}{Seed integer for random number generator, if nmonte > 0}
!!table

!!end 

*/

#include <climits>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <string>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/format.h"
#include "trm/input.h"
#include "trm/rvanal.h"

int main(int argc, char* argv[]){
  
  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("data",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("flow",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("gmax",      Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("mmax",      Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("ncall",     Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("change",    Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("nmonte",    Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("seed",      Subs::Input::LOCAL,   Subs::Input::PROMPT);

    // Get input

    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of times, velocities and uncertainties");

    // Read in data, skipping comments

    Subs::Buffer1D<Subs::rv> data;
    data.load_ascii(rvfile);
    
    if(data.size() < 5) throw std::string("Too few points for correct operation");

    double flo;
    input.get_value("flow", flo, 0.01, 1.e-5, 1.e5, "lower frequency limit (cycles/day)");
    double fhi;
    input.get_value("fhigh", fhi, std::min(1.e5,std::max(flo+10.,20.)), flo, 1.e5, "upper frequency limit (cycles/day)");
    double gmax;
    input.get_value("gmax", gmax, 500., DBL_MIN, DBL_MAX, "maximum conceivable systemic velocity");
    double mmax;
    input.get_value("mmax", mmax, 1.5, DBL_MIN, DBL_MAX, "maximum conceivable companion mass");
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");
    int ncall;
    input.get_value("ncall", ncall, 0, 0, INT_MAX, "number of calls to trapezoidal integrator (minimum if change > 0)");
    double change;
    input.get_value("change", change, 0.001, -DBL_MAX, 0.2, "maximum change in log10 to indicate convergence");
    int nmonte;
    input.get_value("nmonte", nmonte, 0, 0, INT_MAX, "number of Monte Carlo samples/frequency");
    long int seed;
    if(nmonte){
      input.get_value("seed", seed, 5678999L, LONG_MIN, LONG_MAX, "seed integer for random number generator");
      seed = -(long int)(fabs(double(seed)));
    }

    // Add sigma in quadrature, convert to MHJD

    int ngood = 0;
    for(int i=0; i<data.size(); i++){
      if(data[i].z > 0.){
	ngood++;
	data[i].z = sqrt(Subs::sqr(data[i].z)+Subs::sqr(sigma));
      }
    }

    // A rough guard against degeneracy
    if(ngood < 5) throw Rvanal::Rvanal_Error("Number of good points =" + Subs::str(ngood) +
					     " is too few for correct operation");

    // Loop through until enough calls have been made and possibly convergence has occurred
    int ncount = 0;
    double oldval = DBL_MAX, newval = -DBL_MAX;
    Subs::Format form;
    while(ncall == 0 || ncount < ncall || (change > 0. && fabs(newval-oldval) > change)){
      ncount++;
      oldval = newval;
      newval = Rvanal::bayes_trap(data, flo, fhi, gmax, mmax, nmonte, seed, ncount == 1);
      std::cout << "Call number " << ncount << ", log10(integral) = " << form(newval) << std::endl;
    }
  }
  catch(const std::string& msg){
    std::cerr << msg << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}









