/*

!!begin 
!!title  Period detection plotter
!!author T.R. Marsh
!!created March 2001
!!revised 10 Apr 2006
!!index  dprob
!!root   dprob 
!!descr  dprob plots detection probability versus period
!!class  Program
!!css    style.css
!!head1  dprob computes detection probability as function of period

!!emph{dprob} computes the probability of detecting binary motion based on
radial velocity measurements of a star at several time points as
a function of orbital period. It does this for user supplied primary
and secondary masses. It uses the routine pdb which itself is
based upon a simple constant versus random discrimination.

!!head2 Arguments

!!table 
!!arg{ data   } { file of time, radial velocity, uncertainty}
!!arg{ lplow  } { lower limit of log10(P (days))}
!!arg{ lphigh  } { upper limit of log10(P (days))}
!!arg{ nperiod    } { number of periods}
!!arg{ ndiv  } { sub-divisiion factor per period}
!!arg{ m1    } { primary mass (solar masses)}
!!arg{ m2    } { secondary mass (solar masses)}
!!arg{ lfprob} { log10 of false alarm probability}
!!arg{ nph   } { number of phases to use.}
!!arg{ ntrial} { number of trials/phase}
!!arg{ sigma } { uncertainty to add in quadrature}
!!arg{ seed  } { seed integer (only if ntrial > 0)}
!!arg{ device} { plot device, ignore to get print to stdout}
!!table

!!end 

*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <climits>
#include "cpgplot.h"
#include "trm/subs.h"
#include "trm/plot.h"
#include "trm/input.h"
#include "trm/constants.h"
#include "trm/rvanal.h"

int main(int argc, char* argv[]){
  
  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("data",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("lplow",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("lphigh",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nperiod",   Subs::Input::LOCAL , Subs::Input::PROMPT);
    input.sign_in("ndiv",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("m1",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("m2",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("lfprob",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nphase",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("ntrial",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("plot",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("device",    Subs::Input::GLOBAL, Subs::Input::PROMPT);

    // Get input

    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of HJDs, velocities and uncertainties");

    // Read in data, skipping comments and blank lines

    Subs::Buffer1D<Subs::rv> data;
    data.load_ascii(rvfile);
    
    if(data.size() < 2) 
      throw std::string("Number of points =" + Subs::str(data.size()) +
		   " is too few for correct operation");
    std::cerr << data.size() << " radial velocities loaded." << std::endl;

    float lplo, lphi;
    input.get_value("lplow",  lplo, -2.f, -5.f, 5.f, "log10(shortest period (d))");
    input.get_value("lphigh", lphi, std::max(lplo, 3.f), lplo, 5.f, "log10(longest period (d))");
    int   np, ndiv;
    input.get_value("nperiod", np,  100, 2, 10000, "number of period bins");
    input.get_value("ndiv",  ndiv,  10, 1, 10000, "number of sub-divisions/period bin");
    
    float m1, m2;
    input.get_value("m1",  m1, 0.5f, 0.00001f, 1000.f, "primary mass (solar masses)");
    input.get_value("m2",  m2, 0.5f, 0.00001f, 1000.f, "secondary mass (solar masses)");
    float lfprob;
    input.get_value("lfprob",  lfprob, -4.f, -100.f, 0.f, "log10(false alarm probability)");
    int nph, ntrial;
    input.get_value("nphase", nph,  20, 1, 10000, "number of phases to average over/period");
    input.get_value("ntrial", ntrial,  0, 0, 1000000, "number of noise trials/period");
    float sigma;
    input.get_value("sigma", sigma, 0.f, 0.f, 1000.f, "systematic uncertainty to add in quadrature to errors (km/s)");

    // Add sigma in quadrature

    int ngood = 0;
    double mt = 0.;
    for(int i=0; i<data.size(); i++){
      if(data[i].z > 0.){
	ngood++;
	data[i].z = sqrt(Subs::sqr(data[i].z)+Subs::sqr(sigma));
	mt += data[i].x;
      }
    }

    if(ngood < 2)
      throw std::string("Number of good points =" + Subs::str(ngood) +
		   " is too few for correct operation");

    mt /= ngood;
    for(int i=0; i<data.size(); i++) data[i].x -= mt;

    Subs::INT4 seed = 0;
    if(ntrial){
      input.get_value("seed", seed,  Subs::INT4(798792), Subs::INT4(INT_MIN), Subs::INT4(INT_MAX), "seed integer for random number generator");
      if(seed > 0) seed = -seed;
    }

    bool plot;
    input.get_value("plot", plot, true, "plot result (or send to stdout)?");
    std::string device;    
    if(plot){
      input.get_value("device", device, "/xs", "plot device");
    }else{
      std::cout << "#" << std::endl;
      std::cout << "# Output from 'dprob'" << std::endl;
      std::cout << "#" << std::endl;
      std::cout << "# data        = " << rvfile << std::endl;
      std::cout << "# lplow       = " << lplo   << std::endl;
      std::cout << "# lphigh      = " << lphi   << std::endl;
      std::cout << "# nperiod     = " << np     << std::endl;
      std::cout << "# ndiv        = " << ndiv   << std::endl;
      std::cout << "# m1          = " << m1     << std::endl;
      std::cout << "# m2          = " << m2     << std::endl;
      std::cout << "# lfprob      = " << lfprob << std::endl;
      std::cout << "# nphase      = " << nph    << std::endl;
      std::cout << "# ntrial      = " << ntrial << std::endl;
      std::cout << "# sigma       = " << sigma  << std::endl;
      if(ntrial)
	std::cout << "# seed        = " << seed  << std::endl;
      std::cout << "#" << std::endl;
      std::cout << " " << std::endl;
    }

    // Inputs done.

    // Compute chi**2 threshold equivalent to false alarm probability.

    float thresh = Subs::tchi(lfprob,ngood-1);

    std::cerr << "thresh = " << thresh << std::endl;

    // now start computation

    double lp, ph, period, k, detp;
    float xp[np], yp[np];

    for(int ip=0; ip < np; ip++){
      detp = 0.;
      for(int nd=0; nd < ndiv; nd++){
	lp = lplo + (lphi-lplo)*(ip+(nd-(ndiv-1)/2.)/std::max(1,ndiv-1))/(np-1);
	period = pow(10.,lp);
	
	// maximum orbital velocity (km/s) 
	
	k = pow(Constants::TWOPI*(Constants::G*Constants::MSUN)*(m1+m2)/(Constants::DAY*period),1./3.)*m2/(m1+m2)/1.e3;

	for(int iph=0; iph < nph; iph++){
	  if(ntrial){
	    ph = iph/(double)nph;
	  }else{
	    ph = iph/(double)nph/2.;
	  }
	  
	  detp += Rvanal::dpb(data, period, ph, k, ntrial, thresh, seed);
	  
	}
      }
      xp[ip] = lplo + (lphi-lplo)*ip/(np-1);
      yp[ip] = detp/nph/ndiv;
      
      if(!plot) std::cout << xp[ip] << " " << yp[ip] << std::endl;
    }

    if(plot){
      Subs::Plot panel(device);
      cpgsch(1.5);
      cpgscf(2);
      cpgslw(2);
      cpgsci(4);
      cpgenv(lplo,lphi,0.,1.0,0,0);
      cpgsci(2);
      cpglab("log\\d10\\u(Period (d))","Detection probability", " ");
      cpgsci(1);
      cpgline(np,xp,yp);
    }
  }
  catch(const std::string& mess){
    std::cerr << mess << std::endl;
    exit(EXIT_FAILURE);
  }
}








