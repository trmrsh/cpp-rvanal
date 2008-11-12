/*

!!begin 
!!title  Merger diagram plotter
!!author T.R. Marsh
!!created March 2007
!!index  mergers
!!root   mergers
!!descr  plots diagram of merger detection probability
!!class  Program
!!css    style.css
!!head1  mergers plots diagram of merger detection probability

!!emph{mergers} computes the probability of detecting binary motion 
as a function of mass and period assuming 

!!head2 Arguments

!!table 
!!arg{ data   } { file of time, radial velocity, uncertainty. Only the 
times and uncertainties actually matter}
!!arg{ spread } { Spread in times. Times will be sampled from a uniform
distribution of this width around the times read in from the file}
!!arg{ lplow  } { lower limit of log10(P (days))}
!!arg{ lphigh  } { upper limit of log10(P (days))}
!!arg{ nperiod    } { number of periods}
!!arg{ mlow    } { Lower limit of total mass (assumes equal mass ratio)}
!!arg{ mhigh   } { Upper limit of total mass (assumes equal mass ratio)}
!!arg{ nmass  } { number of masses }
!!arg{ lfprob} { log10 of false alarm probability}
!!arg{ nph   } { number of phases to use.}
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
#include "trm_subs.h"
#include "trm_plot.h"
#include "trm_input.h"
#include "trm_constants.h"
#include "trm_rvanal.h"

int main(int argc, char* argv[]){
  
  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("data",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("spread",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("lplow",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("lphigh",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nperiod",   Subs::Input::LOCAL , Subs::Input::PROMPT);
    input.sign_in("mlow",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("mhigh",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nmass",     Subs::Input::LOCAL , Subs::Input::PROMPT);
    input.sign_in("lfprob",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nphase",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
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
      throw std::string("Number of points =" + Subs::str(data.size()) + " is too few for correct operation");

    int ngood = 0;
    for(int i=0; i<data.size(); i++)
	if(data[i].z > 0) ngood++;
    if(ngood < 2) 
      throw std::string("Number of good points =" + Subs::str(ngood) + " is too few for correct operation");

    std::cerr << data.size() << " radial velocities loaded." << std::endl;

    double spread;
    input.get_value("spread", spread, 7., 0., 100000000., "spread in days to apply to the times");

    float lplo;
    input.get_value("lplow",  lplo, -2.f, -5.f, 5.f, "log10(shortest period (d))");
    float lphi;
    input.get_value("lphigh", lphi, std::max(lplo, 3.f), lplo, 5.f, "log10(longest period (d))");
    int   nlp;
    input.get_value("nperiod", nlp,  100, 2, 10000, "number of period bins");

    float mlo;
    input.get_value("mlow",  mlo, 0.2f, 0.01f, 10.f, "lowest total mass (solar))");
    float mhi;
    input.get_value("mhigh", mhi, std::max(mlo, 2.f), mlo, 10.f, "highest total mass (solar))");
    int   nm;
    input.get_value("nmass", nm,  100, 2, 10000, "number of mass bins");
    float lfprob;
    input.get_value("lfprob",  lfprob, -4.f, -100.f, 0.f, "log10(false alarm probability)");
    int nph;
    input.get_value("nphase", nph,  20, 1, 10000000, "number of phases to average over");

    if(ngood < 2)
      throw std::string("Number of good points =" + Subs::str(ngood) +
		   " is too few for correct operation");



    Subs::INT4 seed;
    input.get_value("seed", seed,  Subs::INT4(798792), Subs::INT4(INT_MIN), Subs::INT4(INT_MAX), "seed integer for random number generator");
    if(seed > 0) seed = -seed;

    std::string device;    
    input.get_value("device", device, "/xs", "plot device");

    // Inputs done.

    // Compute chi**2 threshold equivalent to false alarm probability.

    float thresh = Subs::tchi(lfprob,ngood-1);

    std::cerr << "Chi**2 threshold = " << thresh << std::endl;

    // now start computation

    double lp, m, period;
    float detp[nm*nlp];
    int nd = 0;
    double dm  = (mhi-mlo)/(nm-1);
    double dlp = (lphi-lplo)/(nlp-1);
    for(int im=0; im < nm; im++){
	m = mlo + (mhi-mlo)*im/(nm-1);
	for(int ip=0; ip < nlp; ip++, nd++){
	    lp = lplo + (lphi-lplo)*ip/(nlp-1);
	    double dp = 0.;
	    for(int nt=0; nt < nph; nt++){
		double mact  = m  + dm*(Subs::ran3(seed)-0.5);  
		double lpact = lp + dlp*(Subs::ran3(seed)-0.5);  
		period =  pow(10.,lpact);
		// Compute orbital speed in km/s
		float vorb  = pow(Constants::TWOPI*Constants::G*Constants::MSUN*mact/(Constants::DAY*period),1./3)/2000;

		// Random phase
		double phase = Subs::ran3(seed);

		// Compute velocities
		double sumw = 0., sumwy = 0., w;
		for(int nt=0; nt < data.size(); nt++){
		    if(data[nt].z > 0.){
			data[nt].y = vorb*sin(Constants::TWOPI*((data[nt].x+spread*(Subs::ran3(seed)-0.5)/period+phase)));
			w      = 1./Subs::sqr(data[nt].z);
			sumw  += w;
			sumwy += w*data[nt].y;
		    }
		}
		double mean = sumwy/sumw;

		// Compute chisq
		double chisq = 0.;
		for(int nt=0; nt < data.size(); nt++){
		    if(data[nt].z > 0.)
			chisq += Subs::sqr((data[nt].y-mean)/data[nt].z);
		}

		if(thresh < chisq)
		    dp += sqrt(1-thresh/chisq);
	    }
	    detp[nd] = 100.*dp/nph;
	}
    }

    Subs::Plot panel(device);
    cpgsch(1.5);
    cpgscf(2);
    cpgslw(4);
    cpgvstd();
    cpgswin(lplo-dlp/2,lphi+dlp/2,mlo-dm/2,mhi+dm/2);
    cpgsci(2);
    cpglab("log\\d10\\u(Period (d))","Total mass (M\\d\\(2281)\\u)", " ");
    cpgsci(1);
    float tr[] = {lplo-dlp,dlp,0,mlo-dm,0,dm};
    float c[] = {50,60,70,80,90};
//    cpggray(detp,nlp,nm,1,nlp,1,nm,100,0,tr);
    cpgcont(detp,nlp,nm,1,nlp,1,nm,c,5,tr);
    cpgconl(detp,nlp,nm,1,nlp,1,nm, 50., tr, "50", 30, 10);
    cpgconl(detp,nlp,nm,1,nlp,1,nm, 60., tr, "60", 30, 10);
    cpgconl(detp,nlp,nm,1,nlp,1,nm, 70., tr, "70", 30, 10);
    cpgconl(detp,nlp,nm,1,nlp,1,nm, 80., tr, "80", 30, 10);
    cpgsci(2);
    const int NPLOT = 200;
    float x[NPLOT], y[NPLOT];
    for(int i=0; i<NPLOT; i++){
	x[i] = lplo-dlp/2 + (lphi-lplo+dlp)*i/(NPLOT-1);
        y[i] = pow(4e-3,3./5.)*pow(24*pow(10,x[i]),8./5.);
    }
    cpgline(NPLOT,x,y);
    cpgsci(4);
    cpgbox("bcnst",0,0,"bcnst",0,0);
  }
  catch(const std::string& mess){
      std::cerr << mess << std::endl;
      exit(EXIT_FAILURE);
  }
}








