// !!begin rverr rverr computes velocity uncertainty !!title rverr
// !!author T.R. Marsh
// !!head1 rverr computes velocity uncertainty
//
// rverr calculates the radial velocity uncertainty 
//
// !!head2 Arguments
//
// !!table 
// !!arg clo    !! lower limit on counts/pixel
// !!arg chi    !! upper limit on counts/pixel
// !!arg nc     !! number of points
// !!arg back   !! background noise (RMS counts)
// !!arg scale  !! electrons per count
// !!arg wlo    !! first wavelength
// !!arg whi    !! last wavelength
// !!arg nw     !! number of pixels
// !!arg device !! plot device
// !!table
//
// !!end 

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "../subs/subs.h"
#include "../subs/constants.h"
#include "../classes/graphics.h"

struct gauss{
  double w;
  float f, s;
};

int main(int argc, char* argv[]){
  
  const float EFAC = 4.*log(2.); 
  const int NG = 3;
  gauss g[NG];
  g[0].w = 6562.76;
  g[0].f = 10.;
  g[0].s = -0.1;
  g[1].w = 6562.76;
  g[1].f = 20.;
  g[1].s = -0.15;
  g[2].w = 6562.76;
  g[2].f = 40.;
  g[2].s = -0.3;
  
  // get inputs

  if(argc != 10){
    std::cerr << "usage: clo chi nc back scale wlo whi nw device"<< std::endl;
    exit(EXIT_FAILURE);
  }
  float clo   = atof(argv[1]);
  float chi   = atof(argv[2]);
  int   nc    = atoi(argv[3]);
  float back  = atof(argv[4]);
  float scale = atof(argv[5]);
  float wlo   = atof(argv[6]);
  float whi   = atof(argv[7]);
  int   nw    = atoi(argv[8]);
  char *device = argv[9];

  double s, w;
  float f, df, z, x, c = C/1.e3;
  float cts[nc], acc[nc];
  for(int ic=0; ic < nc; ic++){
    cts[ic] = clo + (chi-clo)*ic/(nc-1);
    s = 0.; 
    for(int iw=0; iw<nw; iw++){
      w  = wlo + (whi-wlo)*iw/(nw-1);
      df = 0.;
      f  = 1.;
      for(int ig=0; ig<NG; ig++){
	x   = (w-g[ig].w)/g[ig].f;
	z   = g[ig].s*exp(-EFAC*sqr(x));
	f  += z;
	df += 2.*EFAC*z*x*g[ig].w/c/g[ig].f;
      }
      s += sqr(cts[ic]*df)/(sqr(back) + cts[ic]*f/scale);
    }
    acc[ic] = 1./sqrt(s);
  } 

  graphics pt(device);
  pt.sch(1.5);
  pt.scf(2);
  pt.slw(2);
  pt.sci(4);
  pt.env(clo,chi,0.,100.,0,0);
  pt.sci(2);
  pt.lab("Counts/pixel","RV precision (km/s)", " ");
  pt.sci(1);
  pt.line(nc,cts,acc);
  pt.clos();
  exit(EXIT_SUCCESS);
}








