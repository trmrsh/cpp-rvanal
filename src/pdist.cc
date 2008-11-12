/*

!!begin
!!title pdist
!!author T.R. Marsh
!!date 31 May 2001
!!descr plots distributions of output from sinfalse
!!index pdist
!!root  pdist
!!class Programs
!!class Display
!!css   style.css
!!head1 pdist -- plots distributions of output 

!!emph{pdist} plots distributions of based upon binary files produced
by !!ref{sinfalse.html}{sinfalse}. These are designed to save time and memory
compared to use of perl. !!emph{pdist} plots the fraction of trials that
exceed a certain threshold as a funtion of that threshold. This has an
analytical form for a single frequency. !!emph{pdist} can plot this
function as well as its equivalent for independent frequencies. This then
allows one to determine the effective number of "independent" frequencies
although note that once independence is violated the functional form
deviates from the analytic version.

The analytic formula comes in three forms depending upon the normalisation
used. !!emph{pdist} can cope with three methods. 0 = no normalisation.
i.e. the data errors were assumed correct. 1= normalisation by the
variance a la Horne & Baliunas. 2= normalisation by the variance measured
after the best fitting sine has been removed. All are similar at small
thresholds but deviate very significantly beyond this.

!!emph{pdist} plots 1-sigma errors on the Monte Carlo points. These are
not independent from point to point, but give an idea of the uncertainty
on any one point.

!!head2 Invocation

pdist device x1 x2 y1 y2 (data method nind nbin)!!break

The arguments between brackets can be repeated any number of
times for different files.

!!head2 Arguments

!!table
!!arg{device}{Plot device}
!!arg{x1 x2}{X range}
!!arg{y1 y2}{Y range (in form of exponents, e.g. -2 0)}
!!arg{data}{data file produced by !!emph{sinfalse}}
!!arg{method}{0 = straight, 1=variance normalisation, 2=normalisation
by variance relative to best fit}
!!arg{nindep}{Number of "independent" frequencies}
!!arg{nbin}{Number of histogram bins}
!!arg{end}{End input of  previous four entries}
!!table

!!end

The arguments data, method nind and nbin can be repeated for multiple
files.

!!end

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>
#include <valarray>
#include "cpgplot.h"
#include "trm_subs.h"
#include "trm_input.h"
#include "trm_plot.h"
#include "trm_rvanal.h"
 
int main(int argc, char* argv[]){

  try{

    // Construct Input object

    Subs::Input input(argc, argv, Rvanal::RVANAL_ENV, Rvanal::RVANAL_DIR);

    // sign-in variables (equivalent to ADAM .ifl files)

    input.sign_in("device",    Subs::Input::GLOBAL, Subs::Input::PROMPT);
    input.sign_in("x1",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("x2",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y1",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("y2",        Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("data",      Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("method",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nindep",    Subs::Input::LOCAL,  Subs::Input::PROMPT);
    input.sign_in("nbin",      Subs::Input::LOCAL,  Subs::Input::NOPROMPT);
    input.sign_in("end",       Subs::Input::LOCAL,  Subs::Input::NOPROMPT);

    // Get input

    std::string device;
    input.get_value("device", device, "/xs", "plot device");
    float x1, x2, y1, y2;
    input.get_value("x1", x1, 0.f,   -FLT_MAX, FLT_MAX, "left X limit");
    input.get_value("x2", x2, 100.f, -FLT_MAX, FLT_MAX, "right X limit");
    input.get_value("y1", y1, 0.f,   -FLT_MAX, FLT_MAX, "lower X limit");
    input.get_value("y2", y2, 1.5f,  -FLT_MAX, FLT_MAX, "upper Y limit");

    std::vector<std::string> data;
    std::vector<int> method, nind, nbin;
    int meth, nnd, nbn;
    std::string datafile;

    bool stop = false;
    while(!stop){
      input.get_value("data", datafile, "output.bin", "data file from sinfalse");
      data.push_back(datafile);
      input.get_value("method", meth, 0, 0, 2, "variance normalisation method");
      method.push_back(meth);
      input.get_value("nindep", nnd, 100, 1,  INT_MAX, "number of 'independent frequencies'");
      nind.push_back(nnd);
      input.get_value("nindep", nbn, 100, 1,  10000, "number of histogram bins");
      nbin.push_back(nbn);
      input.get_value("end", stop, false, "input another file name?");
    }

    Subs::Plot plot(device);
    cpgsch(1.5);
    cpgslw(2);
    cpgsci(4);
    cpgscf(2);
    cpgenv(x1,x2,y1,y2,0,20);
    cpgsci(2);
    cpglab("\\gD\\gx\\u2\\d","False alarm probability",
	   "False alarm versus threshold");

    const size_t NPLOT = 200;
    std::valarray<float> xm(NPLOT), ym(NPLOT); 
    float e = 1., xt;
    for(size_t j=0; j<data.size(); j++){

      // Read and prepare the data

      std::vector<float> v;
      std::ifstream fin(data[j].c_str());
      if(!fin)
	throw std::string("Failed to open ") + data[j];
      
      float a;
      size_t ndata;
      if(!(fin.read((char*)&ndata,sizeof(size_t))))
	throw std::string("Failed to read number of data points from ") + data[j] ;
      while(fin.read((char*)&a,sizeof(float))){
	v.push_back(a);
      }
      fin.close();
      std::cout << v.size() << " points read from " << data[j] << std::endl;
      
      x2 = 1.1* (*max_element(v.begin(),v.end()));

      std::valarray<float> xh(nbin[j]), yh(nbin[j]);

      // initialise histogram plot arrays

      for(int i=0; i<nbin[j]; i++){
	xh[i] = x2*i/nbin[j];
	yh[i] = 0.;
      }
    
      // accumulate histogram
      
      for(size_t i=0; i<v.size(); i++)
	yh[size_t(floor(nbin[j]*v[i]/x2))]++;
      
      // convert to prob > xh[i]

      float sum = 0.;
      for(int i=nbin[j]-1; i>0; i--){
	sum  += yh[i];
	yh[i] = sum/v.size();
      }
      yh[0] = 1.;
      
      // Prepare and plot model (if possible)
      
      if(method[j] >=0 && method[j] <= 2){

	if(method[j] == 1){
	  xt = 0.9999*(ndata-1)/2.;
	}else{
	  xt = x2;
	}
	
	for(size_t i=0;  i<NPLOT; i++){
	  xm[i] = x1+(xt-x1)*i/(NPLOT-1);
	  if(method[j] == 0){
	    e     = exp(-xm[i]);
	  }else if(method[j] == 1){
	    e     = pow(1.-2.*xm[i]/(ndata-1),(ndata-3)/2.);
	  }else if(method[j] == 2){
	    e     = 1./pow(1.+2.*xm[i]/(ndata-3),(ndata-3)/2.);
	  }
	  if(e < 1.e-6){
	    ym[i] = log(nind[j]*e)/log(10.);
	  }else{
	    ym[i] = log(1.-pow((1.-e),nind[j]))/log(10.);
	  }
	}
	cpgsci(2);
	cpgline(NPLOT,&xm[0],&ym[0]);
      }

      // Plot data

      cpgsci(1);
      cpgsch(1.);
      float sigma;
      for(int i=0; i<nbin[j]; i++){
	if(yh[i] > 0.){
	  sigma = sqrt(yh[i]*(1.-yh[i])/v.size());
	  cpgmove(xh[i],log(yh[i]-sigma)/log(10.));
	  cpgdraw(xh[i],log(yh[i]+sigma)/log(10.));
	  cpgmove(xh[i]-(x2-x1)/(5.*nbin[j]),log(yh[i])/log(10.));
	  cpgdraw(xh[i]+(x2-x1)/(5.*nbin[j]),log(yh[i])/log(10.));
	  cpgpt1(xh[i],log(yh[i])/log(10.),1);
	}
      }
    }
  }
  catch(const std::string& msg){
    std::cerr << msg << std::endl;
    std::cerr << "Aborting pdist!" << std::endl;
  }
}
