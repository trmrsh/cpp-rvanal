/*
!!START vtest_plot
!!title: vtest_plot
!!head: vtest_plot - evaluates standard variability test
!!usage: vtest_plot data pdist nruns ntrials mass1 mass2 seed x1 x2 y1 y2

<p>
vtest_plot computes chi**2 relative to a constant after
generating RVs from a user-supplied period distribution.
It then generates a plot of detection versus false alarm probability.
If only 1 run is specified, it also plots dashed lines indicating
1-sigma uncertainty. 



<p>
The arguments are as follows:

!!table_start
!!arg data     !! Data file of times, RVs, uncertainties
!!arg pdist    !! Period distribution versus log(period (days))
!!arg nruns    !! The number of tuns, each of ntrials 
!!arg ntrials  !! The number of trials per run
!!arg mass1    !! Mass of observable star (solar masses)
!!arg mass2    !! Mass of companion star (solar masses)
!!arg seed     !! seed for random number generator
!!arg x1       !! start ox x range for plot
!!arg x2       !! end of x range for plot
!!arg y1       !! start of y range for plot
!!arg y2       !! end of y range for plot
!!table_end

!!END */

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#define MODEL_PLT
#define FVECTOR_PLT
#include "../classes/dvector.h"
#include "../classes/fvector.h"
#include "../classes/data.h"
#include "../classes/model.h"
#include "../subs/subs.h"
#include "../subs/constants.h"
#include "../subs/qsort.C"

unsigned long int    data::ndef=20;
unsigned long int   model::ndef=20;
unsigned long int dvector::ndef=20;
unsigned long int fvector::ndef=20;

int main(int argc, char *argv[]){


  unsigned long int i, j;

  // Load arguments

  int narg = 1;

  narg++;
  char data_file[256];
  if(argc < narg){
    std::cout << "Name of data file: ";
    std::cin >> data_file;
  }else{
    strcpy(data_file,argv[narg-1]);
  }

  narg++;
  char theory_file[256];
  if(argc < narg){
    std::cout << "Name of theoretical period distribution: ";
    std::cin >> theory_file;
  }else{
    strcpy(theory_file,argv[narg-1]);
  }

  narg++;
  unsigned long int nruns;
  if(argc < narg){
    std::cout << "Number of runs: ";
    std::cin >> nruns;
  }else{
    nruns = atoi(argv[narg-1]);
  }

  narg++;
  unsigned long int ntrials;
  if(argc < narg){
    std::cout << "Number of trials: ";
    std::cin >> ntrials;
  }else{
    ntrials = atoi(argv[narg-1]);
  }

  narg++;
  float mass1;
  if(argc < narg){
    std::cout << "Primary mass: ";
    std::cin >> mass1;
  }else{
    mass1 = atof(argv[narg-1]);
  }

  narg++;
  float mass2;
  if(argc < narg){
    std::cout << "Secondary mass: ";
    std::cin >> mass2;
  }else{
    mass2 = atof(argv[narg-1]);
  }

  // Orbital speed in km/s for P = 1day

  float vscale = pow(2.*M_PI*G*MSUN*(mass1+mass2)/DAY,1./3.)
    *mass2/(mass1+mass2)/1.e3;

  narg++;
  long int seed;
  if(argc < narg){
    std::cout << "Seed integer: ";
    std::cin >> seed;
  }else{
    seed = atoi(argv[narg-1]);
  }
  seed = -(long int)fabs(seed);

  narg++;
  float x1;
  if(argc < narg){
    std::cout << "Start of x range for plot: ";
    std::cin >> x1;
  }else{
    x1 = atof(argv[narg-1]);
  }

  narg++;
  float x2;
  if(argc < narg){
    std::cout << "End of x range for plot: ";
    std::cin >> x2;
  }else{
    x2 = atof(argv[narg-1]);
  }

  narg++;
  float y1;
  if(argc < narg){
    std::cout << "Start of y range for plot: ";
    std::cin >> y1;
  }else{
    y1 = atof(argv[narg-1]);
  }

  narg++;
  float y2;
  if(argc < narg){
    std::cout << "End of y x range for plot: ";
    std::cin >> y2;
  }else{
    y2 = atof(argv[narg-1]);
  }

  try{

    // Load radial velocity data and period distribution
    // period distribution assumed in terms of log(p) where
    // p is measured in days. convert to cdf form.

    data  rvdata(data_file);
    model pdist(theory_file);    
    pdist = pdist.cdf();
    
    // Extract vectors from data

    dvector x = rvdata.get_vx();
    fvector v = rvdata.get_vy();
    fvector e = rvdata.get_ve();

    // Weights

    fvector w   = inv(sqr(e));
    double sumw = w.sum();

    unsigned long int npix = rvdata.get_npix();
    fvector rv(npix), vb(ntrials);
    float   logp, period, k, phi0, constant, thr;
    double  sum1, chisq;

    // Compute probability of noise showing as much variation as the data

    constant   = (w*v).sum()/sumw;
    std::cout << "logp = " << 
      log10(gammq((npix-1)/2.,(w*sqr(v-constant)).sum()/2.)) << std::endl;

    // Start plot

    graphics_open_device("/xs");
    graphics_set_character_height(1.5);
    graphics_set_font(roman);
    graphics_set_colour(blue);
    graphics_set_scale_draw_axes(x1,x2,y1,y2,0,loglin);
    graphics_set_colour(red);
    graphics_label_axes("False alarm probability", "Detection probability",
			" ");
    graphics_set_colour(white);

    fvector dprob(ntrials+1), thresh(ntrials+1), lfprob(ntrials+1);
    fvector err(ntrials+1);

    for(unsigned long int nrun=0; nrun<nruns; nrun++){
      for(unsigned long int ntrial=0; ntrial<ntrials; ntrial++){

	// Generate RV parameters, period in days
	// (so time must be in days). k in km/s (so
	// velocities must be in km/s), zero point phase.
	// Random inclinations are assumed.
	
	logp   = pdist.xeqv(ran2(seed));
	period = pow(10.,logp);
	k      = vscale*sqrt(1.-sqr(ran2(seed)))/pow(period,1./3.);
	phi0   = ran2(seed);
	
	// Set RV data array, binary motion plus noise
	
	rv   = tofvec(k*sin(2.*M_PI*(x/period-phi0))) + e*gauss2(npix,seed);
	
	// Evaluate best fit constant and Chi**2
	
	constant   = (w*rv).sum()/sumw;
	vb[ntrial] = (w*sqr(rv-constant)).sum();
      }
      
      // Sort into ascending order 
      
      vb.sort();

      // Convert into probability versus threhold level

      thresh[0]  = 0.;
      dprob[0]   = 1.;
      lfprob[0]  = 0.;
      err[0]     = 0.;

      for(i=0;i<ntrials;i++){
	thresh[i+1]  = vb[i];
	lfprob[i+1]  = log10(max(1.e-40F,gammq((npix-1)/2.,vb[i]/2.)));
	dprob[i+1]   = (ntrials-i-1)/(double)ntrials;
	err[i+1]     = sqrt(dprob[i+1]*(1.-dprob[i+1])/ntrials);
      }

      // Plot, with +/- 1sigma uncertainty if only 1

      plot_line(lfprob, dprob);
      if(nruns == 1){
	graphics_set_line_style(dashed);
	plot_line(lfprob,dprob-err);
	plot_line(lfprob,dprob+err);
      }
	
    }
    graphics_close_device();
  }

  catch(int i){
    std::cerr << "Aborted.\n";
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}



