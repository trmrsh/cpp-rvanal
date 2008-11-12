/*

!!begin 
!!title  Sindet -- in a sinusoid detected or not?
!!author T.R. Marsh
!!created  April 2001
!!revised 10 Aug 2005
!!index  sindet
!!root   sindet
!!descr  works out sinusoid detection level
!!class  Program
!!css    style.css
!!head1  sindet - works out sinusoid detection level

!!emph{sindet} searches over a user-defined frequency range for the
frequency of minimum Chi**2. It does so by searching over a grid of
frequencies and minimising when triples with minima are found. It then
uses the Chi**2 found versus that for a constant alone as its test
statistic. i.e. we are looking for the improvement in Chi**2 that
occurs as a reult of adding the sinusoid terms. It then carries out
Monte Carlo trials to estimate the chance of such a level occuring at
random. It reports back the number of trials that exceeded the level
found for the data alone. If after 1000 trials, more than 10 exceeded 
the observed level, then we are below the 99% confidence level.

See below for more details of the operation of the routine.

!!head2 Arguments

!!table 
!!arg{data} { Data file. Expects time (HJD), velocity, uncertainty. 
Lines must start with # to be regarded as comments.}
!!arg{flo } { Lower frequency limt, cycles/unit of x}
!!arg{fhi } { Upper frequency limt, cycles/unit of x.}
!!arg{over} { Over-sampling factor to define the number of frequencies
to compute. If set to 1, the number of frequencies is such that from one
to the next there is a change of 0.1 cycles over the whole baseline of the 
input data. This should normally do.}
!!arg{sigma} {Systematic uncertainty to add to uncertainties, to soften the
effects of unrealistically small error bars that bright objects tend to
give. A typical value may be 0.1 pixels.}
!!arg{ntrial} {Number of trials.}
!!arg{nrep} {Progress reported every nrep trials. =0 to ignore}
!!arg{seed} {Seed integer}
!!arg{method}{Trial method: 'm' = Monte Carlo, using observed
errors, 'r' = randomisation. See below for more details.}
!!table

!!head2 Notes

!!head3 Methods

There are two methods for performing the trials. The Monte Carlo
method generates fake datasets with the same uncertainties as the
data. Randomisation takes the data values and interchanges the
times. It sets the errors = 1 and so cannot be directly compared with
the Monte Carlo results. Randomisation requires a reasonable number of
data points and is usually much more conservative than the Monte Carlo
method. It is in general much slower than the Monte
Carlo method because of the impossibility of time savings outlined
below. If you think your uncertainty estimates are good, then the
Monte Carlo method should be the method of choice. Randomisation
is so conservative in some cases that it leads to non-detection of
perfectly valid periods that are well-detected on the Monte Carlo
trials.

!!head3 Speed

The search over frequency is exhaustive but slow. A technique
used to speed things is only to search the noise trial data if
the Chi**2 compared with a constant exceeds the observed statistic.
This is because the statistic for the any set of data must always 
be less than or equal to the constant fit chi**2. This saves a huge
amount of time if the detection is very solid and it may not require
!!emph{any} searches over the fake data sets. For more marginal data
it is less useful and there can be a huge slow-down in such 
cases. This is when one should use !!emph{nrep} to keep track of
progress.

As a rough rule of thumb, if the observed Delta Chi**2 is less than 
about 100, you may encounter slow progress.

!!end 

*/

#include <climits>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
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
    input.sign_in("flow",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("over",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ntrial",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nreport",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("method",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

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
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");
    int ntrial, nrep;
    input.get_value("ntrial", ntrial, 1000, 1, INT_MAX, "number of trials");
    input.get_value("nreport", nrep, 0, 1, INT_MAX, "number of trials between progress reports");
    Subs::INT4 seed;
    input.get_value("seed", seed,  Subs::INT4(798792), Subs::INT4(INT_MIN), Subs::INT4(INT_MAX), "seed integer for random number generator");
    if(seed > 0) seed = -seed;
    char method;
    input.get_value("method", method, 'm', "mMrR", "trial method");
    method = toupper(method);

    // Add sigma in quadrature, convert to MHJD

    int ngood = 0;
    double sumw = 0.;
    for(int i=0; i<data.size(); i++){
      if(data[i].z > 0.){
	if(method == 'R'){
	  data[i].z = 1;
	}else{
	  data[i].z = sqrt(Subs::sqr(data[i].z)+Subs::sqr(sigma));
	}
	sumw += 1./Subs::sqr(data[i].z);
	ngood++;
      }
    }

    if(ngood < 5)
      throw std::string("Number of good points =" + Subs::str(ngood) + " is too few for correct operation");

    double mean, chic, cbest, fbest, frep;

    // Compute level for data

    double stat;
    mean= 0.;
    for(int j=0; j<data.size(); j++)
      if(data[j].z>0.)
	mean += data[j].y/Subs::sqr(data[j].z);
    mean /= sumw;
    chic = 0.;
    for(int j=0; j<data.size(); j++)
      if(data[j].z>0.) chic += Subs::sqr((data[j].y-mean)/data[j].z);      
    
    Rvanal::bestsin(data, flo, fhi, over, frep, cbest);
    stat = chic - cbest;

    std::cout << "Observed Delta Chi**2 = " << stat << std::endl;
    std::cout << "       Best frequency = " << frep << std::endl;

    // Now do same for noise

    Subs::Buffer1D<Subs::rv> fake(data.size());
    Subs::Buffer1D<bool> chosen(data.size());

    size_t nmore = 0, k = 0;
    for(int i=0; i<ntrial; i++){

      // Generate fake data

      if(method == 'M'){
	
	for(int j=0; j<data.size(); j++){
	  fake[j] = data[j];
	  if(fake[j].z>0.)
	    fake[j].y = fake[j].z*Subs::gauss2(seed);
	}

      }else if(method == 'R'){

	bool diff = false;
	while(!diff){
	  size_t nleft = data.size();
	  for(int j=0; j<data.size(); j++) chosen[j] = false;
	  for(int j=0; j<data.size(); j++, nleft--){
	    k = size_t(floor(nleft*Subs::ran2(seed)))+1;
	    for(int l=0, m=0; l<int(data.size()); l++){
	      if(!chosen[l]){
		m++;
		if(m == int(k)){
		  chosen[l] = true;
		  fake[j]   = data[l];
		  if(j != l) diff = true;
		}
	      }
	    }
	  }
	}
      }

      // Evaluate chi**2 for constant model (=chic)

      mean= 0.;
      for(int j=0; j<data.size(); j++)
	if(fake[j].z>0.)
	  mean += fake[j].y/Subs::sqr(fake[j].z);
      mean /= sumw;
      chic = 0.;
      for(int j=0; j<data.size(); j++)
	if(fake[j].z>0.) chic += Subs::sqr((fake[j].y-mean)/fake[j].z);  

      // If chic <= stat then there is no way that the noise
      // will exceed the level. Thus the (expensive) search
      // is only needed if chic > stat

      if(chic > stat){

	// Search for best fit

	Rvanal::bestsin(fake, flo, fhi, over, fbest, cbest);
	
	if(chic-cbest > stat) nmore++;
      }
      if(nrep)
	if((i+1) % nrep == 0) 
	  std::cout << nmore << " noise trials out of " << i+1
	       << " have exceeded the observed level." << std::endl;
    }
    if(nrep){
      std::cout << "Observed Delta Chi**2 = " << stat << std::endl;
      std::cout << "       Best frequency = " << frep << std::endl;
    }
    std::cout << nmore << " noise trials out of " << ntrial
	 << " have exceeded observed level." << std::endl;
  }
  catch(std::string mess){
    std::cerr << mess << std::endl;

  }
}


















