/*
!!START stest
!!title: stest
!!head: stest - program for false alarm and detection probabilities
!!usage: stest data pdist ntrials flo fhi nf mass1 mass2 seed lfprob

<p>
<p>
stest computes difference in chi**2 between constant and constant+sine
models as a criterion for detecting sinusoidal signals. It first
computes the distribution of this value for pure noise, and then uses
this to compute detection probability.

<p>
Compare with vtest which works on constant model only. stest is much
slower than vtest and will give little or no advantage for small numbers 
of points; it will probably fail for 4 or fewer in fact.

<p>
The arguments are as follows:

!!table_start
!!arg data     !! Data file of times, RVs, uncertainties
!!arg pdist    !! Period distribution versus log(period (days))
!!arg ntrials  !! The number of trials per run
!!arg flo      !! Lowest frequency to search
!!arg fhi      !! Highest frequency to search
!!arg nf       !! Number of frequencies to search
!!arg mass1    !! Mass of observable star (solar masses)
!!arg mass2    !! Mass of companion star (solar masses)
!!arg seed     !! seed for random number generator
!!arg lfprob   !! log10 false alarm probability threshold
!!table_end

!!END */

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "../classes/dvector.h"
#include "../classes/fvector.h"
#include "../classes/data.h"
#include "../classes/model.h"
#include "../subs/subs.h"
#include "../subs/constants.h"
#include "../subs/qsort.C"

unsigned long int fvector::ndef=20;
unsigned long int dvector::ndef=20;
unsigned long int    data::ndef=20;
unsigned long int   model::ndef=20;

int main(int argc, char *argv[]){

  const double TWOPI = 2.*M_PI;
  unsigned long int i, j;
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
  unsigned long int ntrials;
  if(argc < narg){
    std::cout << "Number of trials: ";
    std::cin >> ntrials;
  }else{
    ntrials = atoi(argv[narg-1]);
  }

  narg++;
  double flo;
  if(argc < narg){
    std::cout << "Lowest frequency (cycles/unit time): ";
    std::cin >> flo;
  }else{
    flo = atof(argv[narg-1]);
  }
  flo *= 2.*M_PI;

  narg++;
  double fhi;
  if(argc < narg){
    std::cout << "Highest frequency (cycles/unit time): ";
    std::cin >> fhi;
  }else{
    fhi = atof(argv[narg-1]);
  }
  fhi *= 2.*M_PI;

  narg++;
  unsigned long int nfreqs;
  if(argc < narg){
    std::cout << "Number of frequencies: ";
    std::cin >> nfreqs;
  }else{
    nfreqs = atoi(argv[narg-1]);
  }
  if(nfreqs < 2){
    std::cerr << "Must specify at least two frequencies\n";
    exit(EXIT_FAILURE);
  }
  double df = (fhi-flo)/(nfreqs-1);

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
  float lfprob;
  if(argc < narg){
    std::cout << "log10(False alarm probability): ";
    std::cin >> lfprob;
  }else{
    lfprob = atof(argv[narg-1]);
  }

  try{

    // Load radial velocity data and period distribution
    // period distribution assumed in terms of log(p) where
    // p is measured in days. convert to cdf form.

    data  rvdata(data_file);
    model pdist(theory_file);    
    pdist = pdist.cdf();
    
    // Extract arrays

    unsigned long int npix = rvdata.get_npix();

    dvector  x = rvdata.get_vx();    
    fvector  y = rvdata.get_vy();    
    fvector  e = rvdata.get_ve();
    fvector  w = inv(sqr(e));
    float sumw = w.sum();

    // Trig function recurrences. NR p179. c=cos,s=sin

    double c[npix], ci[npix], s[npix], si[npix], ct;
    double alpha[npix], beta[npix];

    for(i=0;i<npix;i++){
      alpha[i] = 2.*sqr(sin(df*x[i]/2.));
      beta[i]  = sin(df*x[i]);
      ci[i]    = cos((flo-df)*x[i])/e[i];
      si[i]    = sin((flo-df)*x[i])/e[i];
    }

    // First of all compute the threshold equivalent to the
    // specified false alarm probability

    fvector ndchi(ntrials), ncchi(ntrials), rv(npix);
    float s1, s2, sc, ss, scc, ssc, sss, scy, ssy, sy;
    float a1, a2, a3, cf1, cf2, cf3, cf4, cf5, cf6;
    float t, t1, t2, chicon, det, dmax, dchi, k, period, phi0;
    float saver = 0., mu, eps = pow(10.,lfprob);
    unsigned long int ncheck = 20, noff;

    for(unsigned long int ntrial=0; ntrial<ntrials; ntrial++){

      // Time saver -- after a few trials we can guess
      // a level below which we need not bother with 
      // a frequency search -- "saver" 

      if(ntrial == ncheck){
	mu   = eps*ntrial;
	noff = (int)(mu + 6.*sqrt(mu) + 1);
	if(noff < ntrial){
	  fvector temp(ntrial);
	  for(i=0; i<ntrial; i++)
	    temp[i] = ndchi[i];
	  temp.sort();
	  saver = temp[ntrial-noff];
	}
	ncheck *= 2;
      }

      // Set RV data to pure noise, scale by uncertainties

      rv.gnoise(seed);
      rv *= e;

      // constant model chi-square value

      s1 = s2 = 0.;
      for(i=0; i<npix; i++){
	s1 += (t = w[i]*rv[i]);
	s2 += t*rv[i];
      }
      chicon = s2-s1*s1/sumw;

      if(chicon < saver){
	ndchi[ntrial] = chicon;
	ncchi[ntrial] = chicon;
      }else{

	// Now sine model which requires a search in frequency.
	// first initialize trig recurrences (which start a
	// little below start frequency so that when updated once 
	// they are correct)

	for(i=0; i<npix; i++){
	  c[i] = ci[i];
	  s[i] = si[i];
	}
	
	dmax = 0.;
	sy = s1;
	for(unsigned long int nfreq=0; nfreq<nfreqs; nfreq++){
	  
	  // Update trig recurrences and matrices (keep two copies
	  // since one (u) will be modified by svdcmp. Form
	  // normal equation sums.
	  
	  sc = ss = scc = ssc = sss = scy = ssy = 0.;
	  for(i=0; i<npix; i++){
	    c[i] -= (alpha[i]*(ct=c[i])+beta[i]*s[i]);
	    s[i] -= (alpha[i]*s[i]-beta[i]*ct);
	    sc   += (t1 = w[i]*c[i]);
	    ss   += (t2 = w[i]*s[i]);
	    scc  += t1*c[i];
	    ssc  += t2*c[i];
	    sss  += t2*s[i];
	    scy  += t1*rv[i];
	    ssy  += t2*rv[i];
	  }
	  cf1 = scc*sss-ssc*ssc;
	  cf2 = ssc*ss-sc*sss;
	  cf3 = sc*ssc-scc*ss;
	  cf4 = sumw*sss-ss*ss;
	  cf5 = sc*ss-sumw*ssc;
	  cf6 = sumw*scc-sc*sc;
	  det = sumw*cf1+sc*cf2+ss*cf3;
	  
	  a1  = (cf1*sy+cf2*scy+cf3*ssy)/det;
	  a2  = (cf2*sy+cf4*scy+cf5*ssy)/det;
	  a3  = (cf3*sy+cf5*scy+cf6*ssy)/det;
	  
	  // now compute chi-square values differenced with respect
	  // to constant case.
	  
	  for(dchi=chicon, i=0; i<npix; i++)
	    dchi -= w[i]*sqr(rv[i]-a1-a2*c[i]-a3*s[i]);
	  
	  // Keep track of maximum
	  
	  dmax = dmax > dchi ? dmax : dchi;
	  
	}
	ndchi[ntrial] = dmax;
	ncchi[ntrial] = chicon;
      }
      if((ntrial+1) % 1000 == 0) std::cout << "Reached trial " << ntrial+1 << 
				   ", time saver = " << saver << std::endl;
    }
    ndchi.sort();
    ncchi.sort();
    double dip;
    double frac = modf(ntrials*(1.-eps)-1,&dip);
    unsigned long int ip = (unsigned long int)dip;
    float thresh = ndchi[ip]*(1.-frac)+ndchi[ip+1]*frac;
    float tcomp  = ncchi[ip]*(1.-frac)+ncchi[ip+1]*frac;
    std::cout << "Threshold levels = " << thresh << ", " << tcomp << std::endl;

    // Now we do essentially the same but with binary motion
    // added. 

    unsigned long int ndetect = 0;
    for(unsigned long int ntrial=0; ntrial<ntrials; ntrial++){

      // Generate binary parameters

      period = pow(10.,pdist.xeqv(ran2(seed)));
      k = vscale*sqrt(1.-sqr(ran2(seed)))/
	  pow(period,0.3333333333333333);
      phi0   = ran2(seed);

      // Set radial velocities and compute constant
      // model chi**2

      s1 = s2 = 0.;
      for(i=0;i<npix;i++){
	rv[i] = k*sin(TWOPI*(x[i]/period-phi0)) + e[i]*gauss2(seed);
	s1 += (t = w[i]*rv[i]);
	s2 += t*rv[i];
      }
      chicon = s2-s1*s1/sumw;

      // If chicon is not high enough, we will never detect this
      // one so we can save ourselves some time by skipping the
      // frequency search

      if(chicon > thresh){

	// Now sine model which requires a search in frequency.
	// first initialize trig recurrences (which start a
	// little below start frequency so that when updated once 
	// they are correct)
	
	for(i=0; i<npix; i++){
	  c[i] = ci[i];
	  s[i] = si[i];
	}
	
	sy = s1;
	for(unsigned long int nfreq=0; nfreq<nfreqs; nfreq++){
	  
	  // Update trig recurrences and matrices (keep two copies
	  // since one (u) will be modified by svdcmp. Form
	  // normal equation sums.
	  
	  sc = ss = scc = ssc = sss = scy = ssy = 0.;
	  for(i=0; i<npix; i++){
	    c[i] -= (alpha[i]*(ct=c[i])+beta[i]*s[i]);
	    s[i] -= (alpha[i]*s[i]-beta[i]*ct);
	    sc   += (t1 = w[i]*c[i]);
	    ss   += (t2 = w[i]*s[i]);
	    scc  += t1*c[i];
	    ssc  += t2*c[i];
	    sss  += t2*s[i];
	    scy  += t1*rv[i];
	    ssy  += t2*rv[i];
	  }
	  cf1 = scc*sss-ssc*ssc;
	  cf2 = ssc*ss-sc*sss;
	  cf3 = sc*ssc-scc*ss;
	  cf4 = sumw*sss-ss*ss;
	  cf5 = sc*ss-sumw*ssc;
	  cf6 = sumw*scc-sc*sc;
	  det = sumw*cf1+sc*cf2+ss*cf3;
	  
	  a1  = (cf1*sy+cf2*scy+cf3*ssy)/det;
	  a2  = (cf2*sy+cf4*scy+cf5*ssy)/det;
	  a3  = (cf3*sy+cf5*scy+cf6*ssy)/det;
	  
	  // now compute chi-square values differenced with respect
	  // to constant case.
	  
	  for(dchi=chicon, i=0; i<npix; i++)
	    dchi -= w[i]*sqr(rv[i]-a1-a2*c[i]-a3*s[i]);
	  
	  // If ever threshold is exceeded, jump out of loop. This
	  // should save a fair bit of time.
	  
	  if(dchi > thresh){
	    ndetect++;
	    break;
	  }
	}
      }
    }

    std::cout << "Detection probability = " << 100.*ndetect/(double)ntrials << std::endl;
    
  }

  catch(int i){
    std::cerr << "Aborted.\n";
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}




