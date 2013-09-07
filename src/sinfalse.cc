/*

!!begin 
!!title  Sinfalse -- false alarm diagnostic
!!author T.R. Marsh
!!date   April 2001
!!index  sinfalse
!!root   sinfalse
!!descr  sinfalse - generate Chi**2 for noise
!!class  Program
!!css    style.css
!!head1  sinfalse - performs search for best Chi**2 frequency

!!emph{sinfalse} runs Monte Carlo trials in which noisy data is
searched for minimum Chi**2 over frequency. It then reports
a number that depends upon the difference between the minimum 
Chi**2 values based upon a constant versus a constant+sinusoid model 
along with the best fit frequency. The exact number reported depends 
upon the normalisation used, which is determined by the parameter method. 
Which Chi**2 to use depends upon the nature of the data.

!!emph{method=0} gives "no" normalisation,
i.e. the uncertainties on the data are taken to be correct.  In this
case the Chi**2 difference divided by 2 is reported. At a single frequency
this has a false alarm probability distribution = exp(-z) where z is the
value (this can be tested with fhi<flo, see below).

!!emph{method=1} gives normalisation by the sample variance
i.e. it is the Horne & Baliunas normalisation. See method=3 for
the distribution.

!!emph{method=2} gives normalisation by the variance determined relative
to the best sine fit.

All methods have roughly the same distribution for small values
but diverge significantly at large values. method=1 gives smaller,
while method=2 gives large false alarm probabilities than method=0.


!!head Invocation
 
sinfalse data [-n]/[-s over]/[flo fhi over] sigma ntrial seed method [file]";

!!head2 Arguments

!!table 
!!arg{data} { Data file. Expects time (HJD), velocity, uncertainty. 
Lines must start with # to be regarded as comments.}
!!arg{sampling}{'n' = natural fretquencies only, 'o' = natural frequencies
but over-sampled by a specific factor, 'f' = full minimisation}
!!arg{flow} {If sampling = 'f' is specified then the program will
do a search for the best frequency in a specific range.
flow is the lower frequency limt, cycles/unit of x}
!!arg{fhi } { Upper frequency limt, cycles/unit of x. If fhi < flo,
then the search will be done at the single frequency flo}
!!arg{over}{ Over-sampling factor to define the number of frequencies
to compute prior to maximisation. If set to 1, the number of 
frequencies is such that from one to the next there is a change of 0.1 
cycles over the whole baseline of the input data. This should normally do.
If set too small you may miss the best peak entirely.}
!!arg{sigma} {Systematic uncertainty to add to uncertainties, to soften the
effects of unrealistically small error bars that bright objects tend to
give. A typical value may be 0.1 pixels.}
!!arg{ntrial} {Number of trials.}
!!arg{seed} {Seed integer}
!!arg{method}{Normalisation method. 0:  use data variances directly;
1: normalise by sample variance (Horne & Baliunas); 2: normalise
by variance measured from best-fit sinusoid model.}
!!arg{binary}{If yes, output will be in binary format to save space}
!!arg{outfile}{If binary output, this is the name of the file to send it to}
!!table

!!end 

*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <climits>
#include "trm/subs.h"
#include "trm/input.h"
#include "trm/constants.h"
#include "trm/rvanal.h"

int main(int argc, char* argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("data",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("sampling",  Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("flow",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("fhigh",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("over",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("ntrial",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nreport",   Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("seed",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("method",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("binary",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("outfile",    Subs::Input::LOCAL,  Subs::Input::PROMPT);

    // Get input

    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of times, velocities and uncertainties");

    // Read in data, skipping comments

    Subs::Buffer1D<Subs::rv> data(rvfile);
    //    data.load_ascii(rvfile);
    
    if(data.size() < 5) throw std::string("Too few points for correct operation");

    char sampling;
    input.get_value("sampling", sampling, 'f', "nNoOfF", "sampling method ('n', 'o' or 'f')");
    sampling = toupper(sampling);

    double flo, fhi;
    if(sampling == 'F'){
      input.get_value("flow", flo, 0.01, 1.e-5, 1.e5, "lower frequency limit (cycles/day)");
      input.get_value("fhigh", fhi, std::min(1.e5,std::max(flo+10.,20.)), flo, 1.e5, "upper frequency limit (cycles/day)");
    }
    double over;
    if(sampling == 'O' || sampling == 'F')
      input.get_value("over", over, 2., 0.1, 1000., "over sampling factor");
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");
    int ntrial;
    input.get_value("ntrial", ntrial, 1000, 1, INT_MAX, "number of trials");
    Subs::INT4 seed;
    input.get_value("seed", seed,  Subs::INT4(798792), Subs::INT4(INT_MIN), Subs::INT4(INT_MAX), "seed integer for random number generator");
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

    bool fflag;
    input.get_value("binary", fflag, false, "send output as binary to a file?");
    std::string file;
    if(fflag) input.get_value("outfile", fflag, "output.bin", "name of binary format file for output");


    double mean, chic, cbest, fbest, t0;
    float con, amp, z, fb;
    std::ofstream fout;
    if(fflag){
      fout.open(file.c_str());
      fout.write((char*)&nok,sizeof(size_t));
    }
    for(int i=0; i<ntrial; i++){

      // Generate fake data, evaluate chi**2 for best fit
      // constant model (=chic)

      mean = 0.;
      for(int j=0; j<data.size(); j++){
	if(data[j].z>0.){
	  data[j].y = data[j].z*Subs::gauss2(seed);
	  mean += data[j].y/Subs::sqr(data[j].z);
	}
      }
      mean /= sumw;
      chic  = 0.;
      
      for(int j=0; j<data.size(); j++)
	if(data[j].z>0.) chic += Subs::sqr((data[j].y-mean)/data[j].z);      

      // Search for best fit

      if(sampling == 'N'){
	Rvanal::bestsin(data, fbest, cbest);
      }else if(sampling == 'O'){
	Rvanal::bestsin(data, over, fbest, cbest);
      }else{
	if(flo < fhi){
	  Rvanal::bestsin(data, flo, fhi, over, fbest, cbest);
	}else{
	  fbest = flo;
	  cbest = Rvanal::sinfit_chisq(data, fbest, con, amp, t0); 
	}
      }
      if(method == 0){
	z = (chic-cbest)/2.;
      }else if(method == 1){
	z = (nok-3)*(chic-cbest)/chic/2.;
      }else if(method == 2){
	z = (nok-3)*(chic-cbest)/cbest/2.;
      }
      fb = fbest;
      if(fflag){
	fout.write((char*)&z,sizeof(float));
      }else{
	std::cout << z << " " << fb << std::endl;
      }
    }
    if(fflag) fout.close();
  }
  catch(const std::string& mess){
    std::cerr << mess << std::endl;
  }
  catch(const std::bad_alloc&){
    std::cerr << "memory problem in sinfalse" << std::endl;
  }
}


















