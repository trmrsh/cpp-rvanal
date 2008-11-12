/*

!!begin 
!!title  Alias plotter
!!author T.R. Marsh
!!created March 2001
!!revised 02 Jan 2008
!!index  aliases
!!root   aliases
!!descr  aliases plots best periods
!!class  Program
!!css    style.css
!!head1  aliases finds and plots best periods

!!emph{aliases} computes the Chi**2 after fitting a sine wave +
constant to some x, y, error type data. It does so over a user-defined
frequency range and searches for the frequencies which give minima of
Chi**2. It refines each of these and then lists and (optionally) plots 
them over a defined time-interval. The mass function and minimum 
secondary mass (for primary=0.5 Msun) are also reported.

Best periods are those of minimum chi**2. !!emph{aliases} allows the
user to set limits on the number of best periods to report on and on
how far above the global best-fit to consider including a period. It
also allows one to select the best of a group of closely spaced
aliases rather than plotting them all. Typically one is only interested
in peaks within 20 or so of the best value of all.

Plots are drawn from sunset to sunrise. Times when the Sun is 15 degrees
below the horizon are indicated by vertical dashed lines. The date indicated
at the top of the plot corresponds to the UT date at sunset (i.e. the start
of the night). The sinusoids are only drawn over the period during which the
star is below the user-defined airmass limit. Airmasses are indicated along
the top axis. 
When a plot is specified, then the program will also compute uncertainties
on the predicted phases. If these are of order 1, then the phase is too
uncertain to be of any use and then any observation time is OK.

!!head2 Invocation

aliases data flow fhigh [over] delta nmax width fmax sigma plot (device date telescope
stardata star airmass y1 y2)
 
!!head2 Arguments

!!table 
!!arg{data} { Data file. Expects a series of lines starting with
the time (must be HMJD in order for sunrise/sunset times to be plotted correctly), 
velocity, uncertainty. Anything after these will be ignored. Lines must start with 
# or \n to be regarded as comments. Can have anything after the data.}
!!arg{flow } { Lower frequency limt, cycles/unit of x}
!!arg{fhigh } { Upper frequency limt, cycles/unit of x}
!!arg{over} { Over-sampling factor to define the number of frequencies
to compute. If set to 1, the number of frequencies is such that from one
to the next there is a change of 0.1 cycles over the whole baseline of the 
input data. This means that the program searches of order 10*(fhigh-flow)*T*over
frequencies altogether, where T is the difference in time from start to finish. Be warned
that this can be large and take a long time if you are not careful.}
!!arg{delta} {The maximum change in chi**2 compared with the best fit
period to bother with}
!!arg{nmax}  {The maximum number of minima to handle}
!!arg{width} {Fractional width in frequency within which higher chi**2
minima are ignored}
!!arg{fmax} {Maximum mass function to accept. i.e. this is the maximum
value of the minimum mass of the companion to the star of interest. Use this
to eliminate implausible 1000 solar mass companions etc. }
!!arg{sigma} {Systematic uncertainty to add to uncertainties, to soften the
effects of unrealistically small error bars that bright objects tend to
give. A typical value may be 0.1 pixels.}
!!arg{ scale  } { Poor errors leading to large chi**2 may make you think that you
have a single alias if selecting by difference in chi**2. This paramater if set = true
will force the best minimum to have a reduced chi**2 of 1 to account for this.}
!!arg{ plot  } { Make a plot or not. If yes then you need more parameters}
!!arg{ device} { plot device, ignore to get print to stdout}
!!arg{ date  } { Date for plot. This uses observatory parameters to plot 
a suitable visibility period. The date should be of the form 21/12/2001, i.e. dd/mm/yyyy}
!!arg{ telescope}{A string. e.g. WHT. If wrong you will get a print of
possible values.}
!!arg{ stardata}{A data file of star data, see below.}
!!arg{ star}{The name of the star of interest. Needs quotes if it contains
any blanks. This will be searched for (case sensitive) in the stars
of stardata.}
!!arg{ airmass}{Maximum airmass limit}
!!arg{ y1}      { Lower plot limit}
!!arg{ y2}      { Upper plot limit.}
!!table

!!head2 Notes

The program avoids problems with memory by doing the search for
aliases in discrete chunks which overlap slightly. Nevertheless it
does take a while to carry out this search and refinement. Normally
!!emph{over=1} should do but if you are worried that minima might
have been missed, increase the value of !!emph{over}.

The 'stardata' coordinate file must have entries such as:!!break

Star Name!!break
23 12 00 -00 34 22 2000!!break

i.e. the name on one line is followed by the RA, Dec and equinox
and then a blank line.

!!head2 Tips

A typical value of delta is about 20. Note that when you just have
a few data, a single new one can radically alter the alias structure.
This just means that the first few points are not enough to fix the
aliases with any certainty. The phase uncertainty alluded to above should
tell you that this is the case. However as new points are added, the aliases
should clean up to just a few and then it is worth paying attention to
exactly when you acquire the next point. If trying to distinguish between
2 possibilities, then go for the times when they are the most different.

!!end 

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include "cpgplot.h" 
#include "trm_subs.h"
#include "trm_plot.h"
#include "trm_input.h"
#include "trm_constants.h"
#include "trm_telescope.h"
#include "trm_date.h"
#include "trm_time.h"
#include "trm_star.h"
#include "trm_ephem.h"
#include "trm_observing.h"
#include "trm_rvanal.h"

struct Trough{
  double f, c;
};

int main(int argc, char* argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables

    input.sign_in("data",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("flow",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("over",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("delta",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("nmax",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("width",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fmax",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("scale",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("plot",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("device",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("date",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("telescope", Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("stardata",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("star",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("airmass",   Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("ylow",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("yhigh",     Subs::Input::LOCAL,  Subs::Input::PROMPT);

    // Get input

    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of HMJDs, velocities and uncertainties");

    // Read in data, skipping comments

    Subs::Buffer1D<Subs::rv> data;
    data.load_ascii(rvfile);
    
    if(data.size() < 5) throw std::string("Too few points for correct operation = " + Subs::str(data.size()));

    double flo;
    input.get_value("flow", flo, 0.01, 1.e-5, 1.e5, "lower frequency limit (cycles/day)");
    double fhi;
    input.get_value("fhigh", fhi, std::min(1.e5,std::max(flo+10.,20.)), flo, 1.e5, "upper frequency limit (cycles/day)");
    double over;
    input.get_value("over", over, 2., 0.1, 1000., "over sampling factor");
    double delta;
    input.get_value("delta", delta, 10., 0.0001, 10000., "maximum delta Chi**2 relative to best minimum");
    unsigned int nmax;
    input.get_value("nmax", nmax, 20u, 1u, 1000u, "maximum number of minima to compute");
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
	data[i].z = sqrt(Subs::sqr(data[i].z)+Subs::sqr(sigma));
	ngood++;
      }
    }

    if(ngood < 5)
      throw std::string("Number of good points =" + Subs::str(ngood) + " is too few for correct operation");

    bool scale;
    input.get_value("scale", scale, true, "scale to that best peak has reduced chi**2 = 1?");

    // Plotting stuff
    bool plot;
    input.get_value("plot", plot, false, "make plot of aliases for one night?");
    std::string device;
    Subs::Date date;
    Subs::Telescope telescope;
    Subs::Star star;
    double airmass;
    float y1, y2;
    if(plot){
      input.get_value("device", device, "/xs", "plot device");

      std::string sdate;
      input.get_value("date",  sdate, "17 Nov 1961", "date of night");
      date.set(sdate);

      std::string stelescope;
      input.get_value("telescope",  stelescope, "WHT", "telescope");
      telescope = Subs::Telescope(stelescope);
      
      std::string stardata;
      input.get_value("stardata",  stardata, "stars.lis", "file of star positions");

      std::ifstream file(stardata.c_str());
      if(!file) throw std::string("Could not open star data file = " + stardata);
      
      std::string sname;
      input.get_value("star",  sname, "IP Peg", "name of star");
      Subs::Ephem  eph;
      bool searching = true;
      while(searching && file >> star){
	if(star.name() == sname) searching = false;
	if(searching && !(file >> eph)){
	  if(file.bad()){
	    file.close();
	    throw std::string("Star data file has incorrect format");
	  }else{
	    file.clear();
	  }
	}
      }
      file.close();
      if(searching) throw std::string("Failed to find requested star = " + sname + " in file = " + stardata);

      input.get_value("airmass", airmass, 2., 1.001, 10., "maximum airmass to plot");

      double zd = fabs(star.dec()-telescope.latitude());
      if(zd > 360.*acos(1./airmass)/Constants::TWOPI)
	throw std::string("This star will never achieve the airmass limit");

      input.get_value("ylow",  y1, -200.f, -10000.f, 10000.f, "lower y limit for plot (km/s)");
      input.get_value("yhigh", y2, +200.f, -10000.f, 10000.f, "upper y limit for plot (km/s)");
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
    
    // Do the search in chunks to avoid grabbing ludicrous amounts of memory
    const long unsigned int CHUNK = 100000;
    int nchunk = (nfreq-1)/CHUNK+1;
    long unsigned int ngrab = nfreq/nchunk+1;

    // Call main routine
    double freq[ngrab], chisq[ngrab];
    long unsigned int j, k;
    double cbest[nmax], fbest[nmax], cmin, fmin, clow = 1.e30;

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
	  cmin = Subs::brent(freq[i], freq[i-1], freq[i+1], chi, tol, fmin);	
	  Rvanal::sinfit_chisq( data, fmin, con, amp, t0);
	  if(Rvanal::fm(fmin,amp) < fmax){
	    clow  = std::min(clow, cmin);
	    dat.f = fmin;
	    dat.c = cmin;
	    tlist.push_back(dat);
	  }
	}
      }
    }

    double sfac = 1.;
    if(scale) sfac = clow/(ngood-4);

    // Apply criteria to select out best minima
    std::list<Trough>::iterator first = tlist.begin();
    std::list<Trough>::iterator last  = tlist.end();
    std::list<Trough>::iterator lit1, lit2;

    unsigned int nper=0;
    lit1 = first;
    while(lit1 != last){

      if(lit1->c < clow + sfac*delta){
	lit2 = first;
	ok   = true;
	while(lit2 != last && ok){
	  if(lit2->c < lit1->c && fabs(lit1->f - lit2->f) < width*lit2->f) 
	    ok = false;
	  lit2++;
	}
	if(ok){
	  for(j=0; j<nper && lit1->c > cbest[j]; j++);
	  for(k=std::min(nmax-1,nper); k>j; k--){
	    cbest[k] = cbest[k-1];
	    fbest[k] = fbest[k-1];
	  }
	  if(j < nmax){
	    cbest[j] = lit1->c;
	    fbest[j] = lit1->f;
	    if(nper < nmax) nper++;
	  }
	}
      }
      lit1++;
    }

    std::cout << "Number of good data = " << ngood << std::endl;
    
    // Report best peaks
    long unsigned int nok;
    float c[nper], a[nper];
    double z[nper];
    for(nok=0; nok<nper && cbest[nok] < cbest[0]+sfac*delta; nok++){
      Rvanal::sinfit_chisq(data,fbest[nok],c[nok],a[nok],z[nok]);

      std::cout << "Peak " << nok+1 << ", P = " << std::setw(9) << std::setprecision(7) 
	   << 1./fbest[nok]
	   << ", chisq = " << std::setprecision(6) << cbest[nok] 
	   << ", amp = " << std::setprecision(5) << a[nok] << " km/s, fm = " 
	   << Rvanal::fm(fbest[nok],a[nok]) << " Msun, m2min = " 
	   << Rvanal::m2min(fbest[nok],a[nok],0.5) << " Msun."  << std::endl;
    }
    
    // Plot
    if(plot){
      
      long unsigned int i;

      // Sunrise/set times
      Subs::Time tim(date), sunset, twiend, twistart, sunrise;
      tim.add_hour(12.+telescope.longitude()/15.);
      
      if(!Observing::suntime(telescope, tim, -1., sunset))
	throw std::string("Could not find sunset!!");

      // Twilight start/end times
      if(!Observing::suntime(telescope, sunset, -15., twiend))
	throw std::string("Sun never gets to -15 in evening!!");

      tim   = twiend;
      tim.add_hour(0.1);   
      if(!Observing::suntime(telescope, tim, -15., twistart))
	throw std::string("Sun never gets to -15 in morning!!");

      if(!Observing::suntime(telescope, twistart, -1., sunrise))
	throw std::string("Could not find sunset!!");
      
      std::cout << "Sunset to sunrise: " << sunset << " to " << sunrise  << std::endl;
      std::cout << "       Sun < -15.: " << twiend << " to " << twistart << std::endl;

      // Precess to get airmasses better later on
      Subs::Position Sun;
      Sun.set_to_sun(sunset, telescope);

      // Now times when object is visible, limited by sunrise and set
      Subs::Time tfirst, tlast;
      if(!Observing::when_visible(star, telescope, sunset, sunrise, airmass, tfirst, tlast))
	throw std::string(star.name() + " is never visible between sunset and sunrise with airmass < " + Subs::str(airmass));

      double ts1  = sunset.hour();
      double ts2  = ts1 + 24.*(sunrise.mjd()-sunset.mjd());
      double mjd0 = floor(sunset.mjd());

      double tv1  = ts1 + 24.*(tfirst.mjd()-sunset.mjd());
      double tv2  = ts1 + 24.*(tlast.mjd()-sunset.mjd());

      double tw1  = ts1 + 24.*(twiend.mjd()-sunset.mjd());
      double tw2  = ts1 + 24.*(twistart.mjd()-sunset.mjd());

      // Compute heliocentric corrections with
      // single sun position.

      Subs::Time ttemp((tfirst.mjd()+tlast.mjd())/2.);
      Subs::Date *dfirst = &sunset;
      double off = star.tcorr_hel(ttemp, telescope)/Constants::DAY;

      std::cout << "           HMJD at start = " << std::setprecision(14) << tfirst.mjd()+off << std::endl;
      std::cout << "Heliocentric correction = " << std::setprecision(7)   << off << " days." << std::endl;

      double sphi, tm = mjd0 + (tv2+tv1)/48. + off;
      Rvanal::sinerr err;
      for(nok=0; nok<nper && cbest[nok] < cbest[0]+delta; nok++){
	err  = Rvanal::sinfit_errors(data,fbest[nok]);
	sphi = fbest[nok]*sqrt(err.vpp*Subs::sqr(fbest[nok]*(tm-z[nok]))+
			       2.*(tm-z[nok])*err.vpt*fbest[nok]+err.vtt);
	std::cout << "Peak " << nok+1 << ", P = " << std::setfill(' ') << std::setw(9) 
	     << std::setprecision(7) << 1./fbest[nok] << ", chisq = " 
	     << std::setprecision(6) << cbest[nok] 
	     << ", phase uncertainty (1 sigma) = " 
	     << std::setprecision(5) << sphi << " cycles " << std::endl;
      }

      // Now start the plot
      Subs::Plot panel(device);
      cpgsch(1.5);
      cpgscf(2);
      cpgslw(2);
      cpgsci(4);
      cpgvstd();
      Subs::ut_plot(ts1, ts2, y1, y2);
      cpgsci(2);
      std::string title = 
	star.name() + std::string(", ") +
	dfirst->str() + std::string(" (") + 
	telescope.telescope() + std::string(", ") + 
	telescope.site() + std::string(")");
      cpglab("UT (hours)","Radial velocity (km/s)", title.c_str());
      int ci = 1;
      
      const long unsigned int NPLOT = 1000;
      float xp[NPLOT], yp[NPLOT];
      double theta, tt;
      
      for(nok=0; nok<nper && cbest[nok] < cbest[0]+delta; nok++, ci++){
	for(i=0;i<NPLOT; i++){
	  xp[i] = (tt = tv1 + (tv2-tv1)*i/(NPLOT-1));
	  theta = Constants::TWOPI*fbest[nok]*((mjd0+tt/24.+off)-z[nok]);
	  yp[i] = c[nok]+a[nok]*sin(theta);
	}
	cpgsci(std::min(12,ci));
	cpgline(NPLOT,xp,yp);
	cpgsch(0.9);
	std::ostringstream str;
	str << nok+1;
	cpgptxt(tv1-(tv2-tv1)/60.,yp[0],0.,0.5,str.str().c_str());
	cpgptxt(tv2+(tv2-tv1)/60.,yp[NPLOT-1],0.,0.5,str.str().c_str());
	cpgsch(1.5);
      }
      for(int i=0;i<data.size();i++){
	cpgsci(2);
	tt  = 24.*(data[i].x-mjd0-off);
	cpgmove(tt,data[i].y-data[i].z);
	cpgdraw(tt,data[i].y+data[i].z);
	cpgsci(1);
	cpgpt1(tt,data[i].y,17);
      }
      cpgsls(2);
      cpgsci(2);
      cpgmove(tw1,y1);
      cpgdraw(tw1,y2);
      cpgmove(tw2,y1);
      cpgdraw(tw2,y2);
      float air;
      
      cpgsch(1);
      cpgsls(1);
      int nplot = int(floor(10.*(tv2-tv1)/(ts2-ts1)));
      nplot = nplot < 2 ? 2 : nplot;
      for(int i = 0; i<nplot; i++){
	tt = tv1 + (tv2-tv1)*i/(nplot-1);
	ttemp.set(mjd0+tt/24.);
	air = star.altaz(ttemp,telescope).airmass;
	air = floor(100.*air+0.5)/100.;
	std::ostringstream str;
	str << air;
	cpgptxt(tt,y2+(y2-y1)/30.,0.,0.5,str.str().c_str());
	cpgmove(tt,y2-(y2-y1)/40);
	cpgdraw(tt,y2-(y2-y1)/20);
      }
      cpgiden();
    }
  }

  catch(const std::string& mess){
    std::cerr << mess << std::endl;
  }

}








