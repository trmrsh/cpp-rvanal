/*

!!begin 
!!title  Period folder
!!author T.R. Marsh
!!created  March 2001
!!revised 10 Aug 2005
!!index  pfold
!!root   pfold
!!descr  pfold plots phase-folded data versus fit 
!!class  Program
!!css    style.css
!!head1  pfold plots phase-folded data versus fit

!!emph{pfold} searches for a specified alias and then plots the
data on that ephemeris, phase-folded.

!!head2 Arguments

!!table 
!!arg{data} { Data file. Expects time, velocity, uncertainty. 
Lines must start with # to be regarded as comments.}
!!arg{flo } { Lower frequency limt, cycles/unit of x}
!!arg{fhi } { Upper frequency limt, cycles/unit of x}
!!arg{over} { Over-sampling factor to define the number of frequencies
to compute. If set to 1, the number of frequencies is such that from one
to the next there is a change of 0.1 cycles over the whole baseline of the 
input data}
!!arg{width} {Fractional width in frequency within which higher chi**2
minima are ignored}
!!arg{fmax} {Maximum mass function to accept. i.e. this is the maximum
value of the minimum mass of the companion to the star of interest. Use this
to eliminate implausible 1000 solar mass companions etc. }
!!arg{sigma} {Systematic uncertainty to add to uncertainties, to soften the
effects of unrealistically small error bars that bright objects tend to
give. A typical value may be 0.1 pixels (specified in km/s)}
!!arg{nalias}{The alias number to choose starting from 1 as the best}
!!arg{ device} { plot device }
!!arg{ y1}      { Lower plot limit, km/s}
!!arg{ y2}      { Upper plot limit, km/s.}
!!table

!!end 

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include "cpgplot.h" 
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_plot.h"
#include "trm_input.h"
#include "trm_rvanal.h"

struct Trough{
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
    input.sign_in("device",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("ylow",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("yhigh",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);

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
    input.get_value("fmax", fmax, 50., 1.e-5, 10000., "maximum mass funcion to consider (solar masses)");
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");

    // Add sigma in quadrature

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
    std::string device;
    input.get_value("device", device, "/xs", "plot device");
    float yr1, yr2;
    input.get_value("ylow",  yr1, -200.f, -10000.f, 10000.f, "lower y limit (km/s)");
    input.get_value("yhigh", yr2, -200.f, -10000.f, 10000.f, "upper y limit (km/s)");

    // Compute number of frequencies
    
    // Compute number of frequencies
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

    std::list<Trough> tlist;
    Trough dat;
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
	  if(Rvanal::fm(fmin,amp) < fmax){
	    clow  = std::min(clow, cmin);
	    dat.f = fmin;
	    dat.c = cmin;
	    tlist.push_back(dat);
	  }
	}
      }
    }
      
    // Apply criteria to select out best minima
      
    std::list<Trough>::iterator first = tlist.begin();
    std::list<Trough>::iterator last  = tlist.end();
    std::list<Trough>::iterator lit1, lit2;

    int nper=0;
    lit1 = first;
    while(lit1 != last){

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

    std::cerr << "nper = " << nper << std::endl;
    // Compute best fit parameters for selected
    // alias
    
    Rvanal::sinfit_chisq(data,fbest[nalias-1],con,amp,t0);
    std::cout << "Selected alias: P = " << std::setprecision(9) << 1./fbest[nalias-1] 
	 << std::setprecision(5) << " d, chi = " << cbest[nalias-1] 
	 << ", c,a,t0 = " << con << ", " << amp << ", " 
	 << std::setprecision(10) << t0 << std::endl; 

    Subs::Plot plot(device);
    cpgsch(1.5);
    cpgscf(2);
    cpgslw(2);
    cpgsci(4);
    double p1 = 0., p2 = 2., pp, theta;
    cpgenv(p1,p2,yr1,yr2,0,0);
    cpgsci(2);
    std::string title = rvfile + ", P = " + Subs::str(1./fbest[nalias-1]) + " d, \\gx\\u2\\d = " +
      Subs::str(cbest[nalias-1]);

    cpglab("Orbital phase","Radial velocity (km/s)", title.c_str());
      
    const long unsigned int NPLOT = 1000;
    float xp[NPLOT], yp[NPLOT];
      
    cpgsci(3);
    long unsigned int i;
    for(i=0;i<NPLOT; i++){
      xp[i] = (pp = p1 + (p2-p1)*i/(NPLOT-1));
      theta = Constants::TWOPI*pp;
      yp[i] = con+amp*sin(theta);
    }
    cpgsls(2);
    cpgline(NPLOT,xp,yp);
    cpgmove(p1,con);
    cpgdraw(p2,con);
    cpgsls(1);
    cpgsch(1.);

    for(int i=0;i<data.size();i++){
      pp  = fbest[nalias-1]*(data[i].x-t0);
      pp  = pp - floor(pp);
      cpgsci(2);
      cpgmove(pp,data[i].y-data[i].z);
      cpgdraw(pp,data[i].y+data[i].z);
      cpgsci(1);
      cpgpt1(pp,data[i].y,17);
      pp++;
      cpgsci(2);
      cpgmove(pp,data[i].y-data[i].z);
      cpgdraw(pp,data[i].y+data[i].z);
      cpgsci(1);
      cpgpt1(pp,data[i].y,17);
    }
  }
  catch(std::string mess){
    std::cerr << mess << std::endl;
  }
}








