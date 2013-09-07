/*

!!begin
!!author T.R.Marsh
!!title  vtest
!!descr  evaluates the simple variability test
!!root   vtest
!!index  vtest
!!css    style.css
!!head1  vtest - evaluates standard variability test

!!emph{vtest} computes chance of observing chi**2 as large
or larger than observed.

!!head2 Invocation

vtest data sigma

!!head2 Arguments

!!table
!!arg{ data } {Data file of times, RVs, uncertainties}
!!arg{ sigma} {km/s uncertainty to add in quadrature to allow 
for systematic uncertainty, especially for high S/N data}
!!table

!!end

*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "trm/subs.h"
#include "trm/format.h"
#include "trm/constants.h"
#include "trm/input.h"
#include "trm/rvanal.h"

int main(int argc, char *argv[]){

  try{

    // Construct Input object
    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // Sign-in variables
    input.sign_in("data",      Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("sigma",     Subs::Input::GLOBAL, Subs::Input::PROMPT);

    // Get input
    std::string rvfile;
    input.get_value("data", rvfile, "star.rv", "file of HJDs, velocities and uncertainties");

    // Read in data, skipping comments 
    Subs::Buffer1D<Subs::rv> data;
    data.load_ascii(rvfile);
    
    double sigma;
    input.get_value("sigma", sigma, 0., 0., 1000., "systematic uncertainty to add in quadrature to errors (km/s)");

    for(int i=0; i<data.size(); i++)
      if(data[i].z > 0.)
	data[i].z = sqrt(Subs::sqr(data[i].z) + Subs::sqr(sigma));

    double mean = 0., sumw=0.;
    int nok = 0;
    for(int i=0; i<data.size(); i++){
      if(data[i].z > 0.){
	nok++;
	sumw += 1./Subs::sqr(data[i].z);
	mean += data[i].y/Subs::sqr(data[i].z);
      }
    }
    mean /= sumw;

    double chisq = 0.;
    for(int i=0; i<data.size(); i++)
      if(data[i].z > 0.) 
	chisq += Subs::sqr((data[i].y-mean)/data[i].z);

    Subs::Format form;
    std::cout << "log10(pfalse) = " << form(log10(Subs::gammq((nok-1)/2.,chisq/2.))) << std::endl;
    
  }

  catch(std::string err){
    std::cerr << err << std::endl;
  }
}



