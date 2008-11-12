/*

!!begin 
!!title  nthalias -- computes Nth best frequency and chi**2
!!author T.R. Marsh
!!created  April 2001
!!revised 10 Aug 2005
!!index  nthalias
!!root   nthalias
!!descr  nthalias - computes best frequency and chi**2 of sinfit
!!class  Program
!!css    style.css
!!head1  nthalias - performs search for best chi**2 frequency

!!emph{nthalias} searches over a user-defined frequency range
for the Nth frequency of minimum Chi**2. It does so by searching
over a grid of frequencies and minimising when triples with
minima are found. It is also possible to exclude very closely spaced
alias sets to force the program to look at competing alias clusters,
and to eliminate ones that give silly masses.

!!head2 Program call

nthalias data flow fhigh [over] width fmax sigma nalias

!!head2 Arguments

!!table 
!!arg{data} { Data file. Expects time, velocity, uncertainty. 
Lines must start with # to be regarded as comments.}
!!arg{flow  } { Lower frequency limt, cycles/unit of x}
!!arg{fhigh } { Upper frequency limt, cycles/unit of x}
!!arg{[over]} { Over-sampling factor to define the number of frequencies
to compute. If set to 1, the number of frequencies is such that from one
to the next there is a change of 0.1 cycles over the whole baseline of the 
input data. This should normally do.}
!!arg{width} {Fractional width in frequency within which higher chi**2
minima are ignored}
!!arg{fmax} {Maximum mass function to accept. i.e. this is the maximum
value of the minimum mass of the companion to the star of interest. Use this
to eliminate implausible 1000 solar mass companions etc. }
!!arg{sigma} {Systematic uncertainty to add to uncertainties, to soften the
effects of unrealistically small error bars that bright objects tend to
give. A typical value may be 0.1 pixels.}
!!arg{nalias}{The alias number to choose starting from 1 as the best}
!!table

!!end 

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <list>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_rvanal.h"

struct trough{
  double f, c;
};

int main(int argc, char* argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("data",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("flow",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("over",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("width",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("fmax",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("nalias",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);

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
    double over;
    input.get_value("over", over, 2., 0.1, 1000., "over sampling factor");
    double width;
    input.get_value("width", width, 0.01, 0., 1000., "avoidance width near a minimum (cycles/day)");
    double fmax;
    input.get_value("fmax", fmax, 50., 1.e-5, 10000., "maximum mass function to consider (solar masses)");
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");

    // Add sigma in quadrature, convert to MHJD

    int ngood = 0;
    for(int i=0; i<data.size(); i++){
      if(data[i].z > 0.){
	ngood++;
	data[i].z = sqrt(Subs::sqr(data[i].z)+Subs::sqr(sigma));
      }
    }

    if(ngood < 5)
      throw std::string("Number of good points =" + Subs::str(ngood) + " is too few for correct operation");

    int nalias;
    input.get_value("nalias", nalias, 1, 1, 1000, "alias number to compute");

    // Compute number of frequencies needed
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

    long unsigned int nfreq = (long unsigned int)(10.*over*(tmax-tmin)*(fhi-flo));
    nfreq = nfreq > 2 ? nfreq : 3;
    
    // do the search in chunks to avoid grabbing ludicrous amounts of memory

    const long unsigned int CHUNK = 100000;
    int nchunk = (nfreq-1)/CHUNK+1;
    long unsigned int ngrab = nfreq/nchunk+1;

    // Call main routine

    double freq[ngrab], chisq[ngrab];
    int j, k;
    double cbest[nalias], fbest[nalias], cmin, fmin, clow = 1.e30;

    // Construct function object once
    
    Rvanal::Chisq chi(data);

    std::list<trough> tlist;
    trough dat;
    bool ok;
    float con, amp;
    double f1, f2, t0;
    double tol = 1.e-5/(tmax-tmin);


    for(int n=0; n<nchunk; n++){
      f1 = flo + 0.99*(fhi-flo)*n/nchunk;
      f2 = flo + (fhi-flo)*(n+1)/nchunk;

      Rvanal::sinfit_chisq(data,f1,f2,freq,chisq,ngrab);

      for(long unsigned int i = 1; i<ngrab-1; i++){
	if(chisq[i] < chisq[i-1] && chisq[i] < chisq[i+1]){
	    
	  // Refine the peak value
	    
	  cmin = Subs::brent(freq[i], freq[i-1], freq[i+1], 
			     chi, tol, fmin);	

	  Rvanal::sinfit_chisq(data, fmin, con, amp, t0);

	  // Store if OK mass-wise

	  if(Rvanal::fm(fmin,amp) < fmax){
	    clow  = std::min(clow, cmin);
	    dat.f = fmin;
	    dat.c = cmin;
	    tlist.push_back(dat);
	  }
	}
      }
    }
      
    // Now apply condition on separation between peaks to select
    // the one wanted.
      
    std::list<trough>::iterator first = tlist.begin();
    std::list<trough>::iterator last  = tlist.end();
    std::list<trough>::iterator lit1, lit2;

    int nper=0;
    lit1 = first;
    while(lit1 != last){

      // Check that current minimum is not going to be bettered

      lit2 = first;
      ok   = true;
      while(lit2 != last && ok){
	if(lit2->c < lit1->c && fabs(lit1->f - lit2->f) < width*lit2->f){
	  ok = false;
	  break;
	}
	lit2++;
      }

      if(ok){
	for(j=0; j<nper && lit1->c > cbest[j]; j++);
	if(j < nalias){
	  for(k=std::min(nalias-1,nper); k>j; k--){
	    cbest[k] = cbest[k-1];
	    fbest[k] = fbest[k-1];
	  }
	  cbest[j] = lit1->c;
	  fbest[j] = lit1->f;
	  if(nper < nalias) nper++;
	}
      }
      lit1++;
    }

    if(nper == 1){
      if(nper < nalias)
	throw std::string("Only " + Subs::str(nper) + " alias found!");
    }else{
      if(nper < nalias)
	throw std::string("Only " + Subs::str(nper) + " aliases found!");
    }

    // Compute best fit parameters for selected alias
    
    Rvanal::sinfit_chisq(data,fbest[nalias-1],con,amp,t0);

    // Compute errors on fit

    Rvanal::sinerr err;
    err  = Rvanal::sinfit_errors(data,fbest[nalias-1]);

    std::cout << std::endl;
    std::cout << "Alias number " << nalias << "\n" << std::endl;
    std::cout << " Constant term = " 
	 << std::setprecision(6) << con << " +/- " << sqrt(err.vcc) << std::endl;
    std::cout << "Semi-amplitude = "
	 << std::setprecision(6) << amp << " +/- " << sqrt(err.vaa) << std::endl;
    std::cout << "    Phase zero = "
	      << std::setprecision(15) << t0 << " +/- " << std::setprecision(8) << sqrt(err.vtt) << std::endl;
    std::cout << "        Period = "
	 << std::setprecision(9) << 1./fbest[nalias-1] << " +/- " << sqrt(err.vpp) 
	 <<  std::endl;
    std::cout << std::endl;
    std::cout << "Correlation coefficient(c,a)  = " 
	 << err.vca/sqrt(err.vcc*err.vaa) << std::endl;
    std::cout << "Correlation coefficient(c,t0) = " 
	 << err.vct/sqrt(err.vcc*err.vtt) << std::endl;
    std::cout << "Correlation coefficient(c,p)  = " 
	 << err.vcp/sqrt(err.vcc*err.vpp) << std::endl;
    std::cout << "Correlation coefficient(a,t)  = " 
	 << err.vat/sqrt(err.vaa*err.vtt) << std::endl;
    std::cout << "Correlation coefficient(a,p)  = " 
	 << err.vap/sqrt(err.vaa*err.vpp) << std::endl;
    std::cout << "Correlation coefficient(t,p)  = " 
	 << err.vpt/sqrt(err.vtt*err.vpp) << std::endl;
    std::cout << std::endl;
    std::cout << "Chi**2 = " << cbest[nalias-1] << " for " << ngood 
	 << " valid points.\n\n" << std::endl;
  }
  catch(std::string mess){
    std::cerr << mess << std::endl;
  }
}








