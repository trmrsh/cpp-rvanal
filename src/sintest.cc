/*

!!begin 
!!title   Sintest -- test of sine statistic
!!author  T.R. Marsh
!!created April 2001
!!revised 10 Aug 2005
!!index  sintest
!!root   sintest
!!descr  sintest - computes best frequency and Chi**2 of sinfit
!!class  Program
!!css    style.css
!!head1  sintest - performs search for best Chi**2 frequency

!!emph{sintest} searches over a user-defined frequency range
for the frequency of minimum Chi**2. It does so by searching
over a grid of frequencies and minimising when triples with
minima are found. It comes back with the change in chi**2 relative
to a constant model, the Monte Carlo model and the best fit frequency.

!!head2 Arguments

!!table 
!!arg{data} { Data file. Expects time (HJD), velocity, uncertainty. 
Lines must start with # to be regarded as comments.}
!!arg{pprob}{Period probability distribution, versus log10(P). It
must be monotonic in P and start and end with zero probability.}
!!arg{mass} {Effective mass of secondary = m2**3/(m1+m2)**2}
!!arg{flo } { Lower frequency limt, cycles/unit of x}
!!arg{fhi } { Upper frequency limt, cycles/unit of x. If fhi < flo,
then the search will be done at the single frequency flo}
!!arg{over} { Over-sampling factor to define the number of frequencies
to compute. If set to 1, the number of frequencies is such that from one
to the next there is a change of 0.1 cycles over the whole baseline of the 
input data. This should normally do.}
!!arg{sigma} {Systematic uncertainty to add to uncertainties, to soften the
effects of unrealistically small error bars that bright objects tend to
give. A typical value may be 0.1 pixels.}
!!arg{ntrial} {Number of trials.}
!!arg{seed} {Seed integer}
!!table

!!end 

*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <climits>
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_constants.h"
#include "trm_rvanal.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("data",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("pprob",     Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("mass",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("flow",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("over",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ntrial",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",      Subs::Input::LOCAL,  Subs::Input::PROMPT);

    // Get input

    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of times, velocities and uncertainties");

    // Read in data, skipping comments

    Subs::Buffer1D<Subs::rv> data;
    data.load_ascii(rvfile);
    
    if(data.size() < 5) throw Rvanal::Rvanal_Error("Too few points for correct operation");

    // Read in period distribution
    std::string ppfile;
    input.get_value("pprob", ppfile, "per.dist", "period distribution file");

    Subs::Buffer1D<Subs::xy<double,float> > pdist;
    pdist.load_ascii(ppfile);
    
    size_t ndist = pdist.size();
    if(ndist < 3) throw std::string("Too few points for correct operation");
    if(pdist[0].y != 0. && pdist[ndist-1].y != 0.)
      throw std::string("Distribution must start and end with PDF=0.");

    float  mass;
    input.get_value("mass", mass, 0.5f, 0.0000001f, 1000.f, "effective mass of secondary star"); 
    double flo, fhi;
    input.get_value("flow", flo, 0.01, 1.e-5, 1.e5, "lower frequency limit (cycles/day)");
    input.get_value("fhigh", fhi, std::min(1.e5,std::max(flo+10.,20.)), flo, 1.e5, "upper frequency limit (cycles/day)");
    double over;
    input.get_value("over", over, 2., 0.1, 1000., "over sampling factor");
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");
    int ntrial;
    input.get_value("ntrial", ntrial, 1000, 1, INT_MAX, "number of trials");
    Subs::INT4 seed;
    input.get_value("seed", seed, Subs::INT4(798792), Subs::INT4(INT_MIN), Subs::INT4(INT_MAX), "seed integer for random number generator");
    if(seed > 0) seed = -seed;
    int method;
    input.get_value("method", method, 0, 0, 2, "trial method");
    method = toupper(method);

    // Add sigma in quadrature, convert to MHJD

    double sumw = 0.;
    size_t nok=0;
    for(int i=0; i<data.size(); i++){
      if(data[i].z > 0.){
	data[i].z = sqrt(Subs::sqr(data[i].z)+Subs::sqr(sigma));
	sumw += 1./Subs::sqr(data[i].z);
	nok++;
      }
    }

    if(nok < 5)
      throw std::string("Number of good points =" + Subs::str(nok) + " is too few for correct operation");


    // Convert period distribution to cdf

    double sum = 0.;
    double cdf[ndist];
    cdf[0] = 0.;
    for(size_t j=1; j<ndist; j++){
      sum    += (pdist[j-1].y+pdist[j].y)*(pdist[j].x-pdist[j-1].x)/2.;
      cdf[j]  = sum; 
    }
    for(size_t j=1; j<ndist; j++) cdf[j] /= sum;

    unsigned long int k;
    double y, lp, period, vorb, cosi, phi, mean, chic;
    double cbest, fbest, t0;
    float con, amp;
    for(int i=0; i<ntrial; i++){

      // Generate period

      y = Subs::ran2(seed);
      k = Subs::locate(cdf,ndist,y);
      if(k == 0 || k == ndist) throw std::string("Period out of range");
      lp = (pdist[k-1].x*(cdf[k]-y)+pdist[k].x*(y-cdf[k-1]))/
	(cdf[k]-cdf[k-1]);
      period = pow(10.,lp);

      // Maximum orbital velocity (km/s) 
      
      vorb = pow(Constants::TWOPI*(Constants::G*Constants::MSUN)*mass/(Constants::DAY*period),0.333333333333)/1.e3;
      
      // Project

      cosi  = Subs::ran2(seed);
      vorb *= sqrt(1.-Subs::sqr(cosi));

      // Generate a phase

      phi = Subs::ran2(seed);

      // Generate fake data, evaluate chi**2 for best fit
      // constant model (=chic)

      mean= 0.;
      for(int j=0; j<data.size(); j++){
	if(data[j].z>0.){
	  data[j].y = vorb*sin(Constants::TWOPI*(data[j].x/period-phi))
	    + data[j].z*Subs::gauss2(seed);
	  mean += data[j].y/Subs::sqr(data[j].z);
	}
      }
      mean /= sumw;
      chic = 0.;
      for(int j=0; j<data.size(); j++)
	if(data[j].z>0.) chic += Subs::sqr((data[j].y-mean)/data[j].z);      

      // Search for best fit

      if(flo < fhi){
	Rvanal::bestsin(data, flo, fhi, over, fbest, cbest);
      }else{
	fbest = flo;
	cbest = Rvanal::sinfit_chisq(data, fbest, con,amp, t0); 
      }
      std::cout << chic-cbest << " " << 1/period << " " << fbest << std::endl;
    }
  }
  catch(const std::string& mess){
    std::cerr << mess << std::endl;
  }
  catch(const std::bad_alloc&){
    std::cerr << "memory problem in sintest" << std::endl;
  }
}


















