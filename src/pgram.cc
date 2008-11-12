/*

!!begin 
!!title  Periodogram
!!author T.R. Marsh
!!created March 2001
!!revised 10 Aug 2005
!!index  pgram
!!root   pgram 
!!descr  pgram plots a periodogram
!!class  Program
!!class  Display
!!css    style.css
!!head1  pgram computes a periodogram

!!emph{pgram} computes and plots the Chi**2 after fitting a sine wave +
constant to some x, y, error type data. It does this at regularly spaced
series of frequencies for speed using recurrence relations. It does not 
use the Press & Rybicki fast algorithm because this fails when data are 
very irregularly spaced. The lower the chi**2, the better the period. A
dashed line is plotted indicating the chi**2 of a constant only.

!!head2 Arguments

!!table 
!!arg{ data  }  { Data file}
!!arg{ flo   }  { Lower frequency limt, cycles/unit of x}
!!arg{ fhi   }  { Upper frequency limt, cycles/unit of x}
!!arg{ over }   { Oversampling factor to set the number 
of frequencies. 1 is such that there will be 0.1 cycles
change over the whole baseline from one frequency to the next,
which is about the minimum necessary to sample the periodogram properly.}
!!arg{sigma}{Systematic uncertainty to add in quadrature to uncertainties}
!!arg{ plot}  { yes to plot, no for output to stdout}
!!arg{ device}  { Plot device (if plot = true)}
!!arg{ y1}      { Lower plot limit (if plot = true)}
!!arg{ y2}      { Upper plot limit (if plot = true)}
!!arg{freq}     { Plot frequency, else period (if plot = true)}
!!arg{logx}     { Plot X scale logarithmically (if plot = true)}
!!table

Note that the program always spaces points uniformly in frequency, even if you
plot versus period because this allows fast recurrence relations to be applied
and peaks in periodograms tend to be spaced more uniformly in frequency than period.

!!end 

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <string>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_format.h"
#include "trm_plot.h"
#include "trm_input.h"
#include "trm_rvanal.h"

int main(int argc, char* argv[]){
  
  try{

    // Construct Input object
    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // Sign-in variables
    input.sign_in("data",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("flow",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("over",      Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("plot",      Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("device",    Subs::Input::GLOBAL,  Subs::Input::PROMPT);
    input.sign_in("ylow",      Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("yhigh",     Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("freq",      Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("logx",      Subs::Input::LOCAL,   Subs::Input::PROMPT);
    input.sign_in("xunits",    Subs::Input::LOCAL,   Subs::Input::NOPROMPT);

    // Get input

    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of times, velocities and uncertainties");

    // Read in data, skipping comments
    Subs::Buffer1D<Subs::rv> data;
    data.load_ascii(rvfile);
    
    if(data.size() < 4) throw std::string("Too few points for correct operation");

    double flo;
    input.get_value("flow", flo, 0.01, 1.e-5, 1.e5, "lower frequency limit (cycles/day)");
    double fhi;
    input.get_value("fhigh", fhi, std::min(1.e5,std::max(flo+10.,20.)), flo, 1.e5, "upper frequency limit (cycles/day)");
    double over;
    input.get_value("over", over, 2., 0.1, 1000., "over sampling factor");
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

    if(ngood < 4)
      throw std::string("Number of good points =" + Subs::str(ngood) + " is too few for correct operation");

    bool plot, logx, pfreq;
    input.get_value("plot", plot, false, "plot periodogram (else send to stdout)?");
    std::string device, xunits;
    float yr1, yr2;
    if(plot){
      input.get_value("device", device, "/xs", "plot device");
      input.get_value("ylow",   yr1, 0.f, -FLT_MAX, FLT_MAX, "lower y limit for plot");
      input.get_value("yhigh",  yr2, float(10.*ngood), -FLT_MAX, FLT_MAX, "upper y limit for plot");
      input.get_value("freq",   pfreq, true, "plot versus frequency (else period)?");
      input.get_value("logx",   logx, false, "make X axis logarithmic?");
      input.get_value("xunits", xunits, "d", "name for units of period");
    }

    // Compute number of frequencies
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
    long unsigned int nfreq = (long unsigned int)(10.*over*(tmax-tmin)*(fhi-flo));

    nfreq = nfreq > 2 ? nfreq : 3;

    // Do the computation in chunks to avoid grabbing ludicrous amounts of memory
    const long unsigned int CHUNK = 10000;
    int   nchunk = (nfreq-1)/CHUNK + 1;
    const long unsigned int NGRAB = nfreq/nchunk+1;

    // Call main routine
    double freq[NGRAB],  chisq[NGRAB];
    float  fplot[NGRAB], cplot[NGRAB];

    Subs::Plot panel;

    // Initialise plot
    if(plot){
      
      double s1, s2, s3, wgt;
      s1 = s2 = s3 = 0.;
      for(int i=0; i<data.size(); i++){
	if(data[i].z>0.){
	  wgt = 1./Subs::sqr(data[i].z);
	  s1 += wgt;
	  s2 += wgt*data[i].y;
	  s3 += wgt*Subs::sqr(data[i].y);
	}
      }
      float level = s3-s2*s2/s1;
      
      panel.open(device);
      cpgsch(1.5);
      cpgscf(2);
      cpgslw(2);
      cpgsci(4);
      if(logx && pfreq){
	cpgenv(log10(flo),log10(fhi),yr1,yr2,0,10);
      }else if(logx) {
	cpgenv(-log10(fhi),-log10(flo),yr1,yr2,0,10);
      }else if(pfreq) {
	cpgenv(flo,fhi,yr1,yr2,0,0);
      }else{
	cpgenv(1/fhi,1/flo,yr1,yr2,0,0);
      }
      cpgsci(2);
      std::string xlabel;
      if(pfreq)
	xlabel = "Frequency (cycles/" + xunits + ")";
      else
	xlabel = "Period (" + xunits + ")";

      cpglab(xlabel.c_str(), "\\gx\\u2\\d", " ");
      cpgsci(2);
      cpgsls(2);
      cpgmove(flo,level);
      cpgdraw(fhi,level);
      cpgsci(1);
      cpgsls(1);
    }

    double f1, f2, fl = 0, cl = 0;
    Subs::Format form;

    for(int n=0; n<nchunk; n++){

      // Choose ranges to ensure uniform steps
      f1 = flo + (fhi-flo)*(NGRAB*n)/(NGRAB*nchunk-1);
      f2 = flo + (fhi-flo)*(NGRAB*(n+1)-1)/(NGRAB*nchunk-1);

      Rvanal::sinfit_chisq(data,f1,f2,freq,chisq,NGRAB);

      // Connect between chunks
      if(plot){
	if(n){
	  cpgmove(fl, cl);
	  if(logx && pfreq){
	    cpgdraw(log10(freq[0]), chisq[0]);
	  }else if(logx){
	    cpgdraw(-log10(freq[0]), chisq[0]);
	  }else if(pfreq){
	    cpgdraw(freq[0], chisq[0]);
	  }else{
	    cpgdraw(1/freq[0], chisq[0]);
	  }
	}

	if(logx && pfreq){
	  fl = log10(freq[NGRAB-1]);
	}else if(logx){
	  fl = -log10(freq[NGRAB-1]);
	}else if(pfreq){
	  fl = freq[NGRAB-1];
	}else{
	  fl = 1/freq[NGRAB-1];
	}

	cl = chisq[NGRAB-1];
      }

      for(long unsigned int i=0; i<NGRAB; i++){
	if(plot){
	  if(logx && pfreq){
	    fplot[i] = float(log10(freq[i]));
	  }else if(logx){
	    fplot[i] = float(-log10(freq[i]));
	  }else if(pfreq){
	    fplot[i] = float(freq[i]);
	  }else{
	    fplot[i] = float(1/freq[i]);
	  }
	  cplot[i] = float(chisq[i]);

	}else{
	  std::cout << form(freq[i]) << " " << form(chisq[i]) << std::endl;
	}
      }
      
    
      if(plot)
	cpgline(NGRAB,fplot,cplot);
    }
  }
  catch(std::string msg){
    std::cerr << msg << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}









