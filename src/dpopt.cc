/*

!!begin 
!!title  Detection optimisation
!!author T.R. Marsh
!!root   dpopt
!!descr  genetic algorithm optimisationof detection.
!!css   style.css
!!index  dpopt
!!head1 dpopt optimises detection probability 

!!emph{dpopt} carries out a genetic algorithm optimisation of mean detection 
probability for a particular number of observations picked out of a
discrete set of times. It does so by "evolving" a set of time sequences
to maximise their efficiency at picking up binaries. It calculates the
detection probability by integrating over a grid equally spaced in log(P)
and phase. The user specifies the number of points in each dimension. For
each grid point the program generates velocities assuming an edge-on orbit
and then calculates the inclination for which there would be no detection.
Detection is specified by looking for chi**2 to exceed a certain threshold
value defined by a use specified false alarm probability. The critical
inclination is then used to account for random orientations. Noise can be 
added too and as it does affect the detection probability, but it needs to
averaged over so there is also an optional number of trials per grid point. 
Note however that in this mode there is a component of the computation that 
grows as nperiod*nphase*ntrial which can get expensive. If ntrial=0 then 
noise is not added. This can be close to the right value and is much faster.

Times are chosen from a discrete set of possible times. Thus the sequence
3, 12, 34 means the time slots 3, 12 and 34. No sequence can have the same 
time slot repeated. Individuals of a new generation are made either by 
crossing two individuals of the old generation (selected depending upon 
their fitness) or by straight copying of old generation individuals. Crossing 
means taking some of the slots from one and some from another and combining 
them. Finally to allow flexibility, three different types of mutation are
supported: (1) a single slot is shifted within a small range, (2) a single
slot is shifted completely at random, (3) all slots are re-randomised just
as the population is set up initially. One would typically want to use more
than one of these.

!!head2 Arguments

!!table
!!arg{ndat   }{number of data points (at least 2)}
!!arg{sigma  }{uncertainty in km/s on each point}
!!arg{m1     }{primary (brightest) star's mass in solar masses}
!!arg{m2     }{secondary (faintest) star's mass in solar masses}
!!arg{lfprob }{log10(false alarm probabiliy)}
!!arg{nperiod}{number of grid points in log(P)}
!!arg{nphase }{number of grid points in phase}
!!arg{ntrial }{number of trials per grid point (0 for no noise correction)}
!!arg{npop   }{population size (kept constant)}
!!arg{ngen   }{number of generations}
!!arg{cross  }{fraction of new generation produced by crossing 
(typically > 50%)}
!!arg{mutate1}{rate of local mutations affecting a single slot }
!!arg{mutate2}{rate of random mutations affecting a single slot}
!!arg{mutate3}{rate of random mutations affecting all slots}
!!arg{mdiff  }{maximum slot shift for local mutations.}
!!arg{select }{selection pressure factor}
!!arg{seed   }{random number seed}
!!arg{pfrac  }{fraction to plot}
!!arg{device }{plot device}
!!table

Mutation rates should probably be no more than a few percent. Crossover rates
should be substantial. The selection pressure factor should be of order 0.5.

!!end 

*/

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "trm/subs.h"
#include "trm/constants.h"
#include "trm/sequence.h"
#include "dprob.h"

int main(int argc, char* argv[]){
  
  double TWOPI = 2.*M_PI;

  // get arguments

  int narg = 1;

  int ndat;
  if(narg < argc){
    ndat    = atoi(argv[narg]);
  }else{
    std::cerr << "Number of points (>1): ";
    std::cin >> ndat;
  }
  if(ndat < 2){
    std::cerr << "ndat should be at least 2" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  float sigma;
  if(narg < argc){
    sigma = atof(argv[narg]);
  }else{
    std::cerr << "Uncertainty per point (km/s, >0): ";
    std::cin >> sigma;
  }
  if(sigma <= 0.){
    std::cerr << "sigma must be > 0" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  float m1;
  if(narg < argc){
    m1 = atof(argv[narg]);
  }else{
    std::cerr << "Mass of primary (solar masses) (>0): ";
    std::cin >> m1;
  }
  if(m1 <= 0.){
    std::cerr << "m1 must be > 0" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  float m2;
  if(narg < argc){
    m2 = atof(argv[narg]);
  }else{
    std::cerr << "Mass of secondary (solar masses) (>0): ";
    std::cin >> m2;
  }
  if(m2 <= 0.){
    std::cerr << "m2 must be > 0" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  float lfprob;
  if(narg < argc){
    lfprob = atof(argv[narg]);
  }else{
    std::cerr << "Log10(false alarm probability) (<0): ";
    std::cin >> lfprob;
  }
  if(lfprob >= 0.){
    std::cerr << "lfprob must be < 0" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  int nperiod;
  if(narg < argc){
    nperiod    = atoi(argv[narg]);
  }else{
    std::cerr << "Number of periods to integrate over (>1): ";
    std::cin >> nperiod;
  }
  if(nperiod < 2){
    std::cerr << "nperiod should be at least 2" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  int nphase;
  if(narg < argc){
    nphase    = atoi(argv[narg]);
  }else{
    std::cerr << "Number of phases to integrate over (>0): ";
    std::cin >> nphase;
  }
  if(nphase < 1){
    std::cerr << "nphase should be at least 1" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  int ntrial;
  if(narg < argc){
    ntrial    = atoi(argv[narg]);
  }else{
    std::cerr << "Number of trials/period-phase point (>=0): ";
    std::cin >> ntrial;
  }
  if(ntrial < 0){
    std::cerr << "ntrial should be at least 1" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  int npop;
  if(narg < argc){
    npop    = atoi(argv[narg]);
  }else{
    std::cerr << "Population size (>1): ";
    std::cin >> npop;
  }
  if(npop < 2){
    std::cerr << "npop should be at least 2 and probably much higher" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  int ngen;
  if(narg < argc){
    ngen    = atoi(argv[narg]);
  }else{
    std::cerr << "Number of generations (>0): ";
    std::cin >> ngen;
  }
  if(ngen < 1){
    std::cerr << "ngen should be at least 1" << std::endl;
    exit(EXIT_FAILURE);
  }

  narg++;
  float cross;
  if(narg < argc){
    cross = atof(argv[narg]);
  }else{
    std::cerr << "Crossover rate (0 to 1, e.g. 0.8): ";
    std::cin >> cross;
  }
  if(cross < 0. || cross > 1){
    std::cerr << "cross must be between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  float mutate1;
  if(narg < argc){
    mutate1 = atof(argv[narg]);
  }else{
    std::cerr << "Local single slot mutation rate (0 to 1): ";
    std::cin >> mutate1;
  }
  if(mutate1 < 0. || mutate1 > 1){
    std::cerr << "mutate1 must be between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  float mutate2;
  if(narg < argc){
    mutate2 = atof(argv[narg]);
  }else{
    std::cerr << "Random single slot mutation rate (0 to 1): ";
    std::cin >> mutate2;
  }
  if(mutate2 < 0. || mutate2 > 1){
    std::cerr << "mutate2 must be between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  float mutate3;
  if(narg < argc){
    mutate3 = atof(argv[narg]);
  }else{
    std::cerr << "Random all slot mutation rate (0 to 1): ";
    std::cin >> mutate3;
  }
  if(mutate3 < 0. || mutate3 > 1){
    std::cerr << "mutate3 must be between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  int mdiff;
  if(narg < argc){
    mdiff = atoi(argv[narg]);
  }else{
    std::cerr << "Max shift for local mutations (>0): ";
    std::cin >> mdiff;
  }
  if(mdiff < 1){
    std::cerr << "mdiff must be > 0" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  float select;
  if(narg < argc){
    select = atof(argv[narg]);
  }else{
    std::cerr << "Selection pressure (0 to 1, e.g. 0.5): ";
    std::cin >> select;
  }
  if(select < 0. || select > 1){
    std::cerr << "select must be between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  long int seed;
  if(narg < argc){
    seed = atoi(argv[narg]);
  }else{
    std::cerr << "Random number seed: ";
    std::cin >> seed;
  }
  seed = -(int)fabs(seed);

  narg++;
  float pfrac;
  if(narg < argc){
    pfrac = atof(argv[narg]);
  }else{
    std::cerr << "Fraction of best to plot (0 to 1): ";
    std::cin >> pfrac;
  }
  if(pfrac < 0. || pfrac > 1){
    std::cerr << "pfrac be between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }  

  narg++;
  char device[256];
  if(narg < argc){
    strcpy(device,argv[narg]);
  }else{
    std::cerr << "Plot device: ";
    std::cin >> device;
  }

  // array of possible times in days.

  const int NTIMES = 10000;
  double time[NTIMES];
  bool plot[NTIMES];
  for(int i=0;i<NTIMES;i++){
    time[i] = 0.01*i;
  }

  // expected period distribution

  const int NPER = 100;
  double lper[NPER], pprob[NPER];
  double peak = 1.47712, fwhm = 0.03;
  double lp1=peak-2.*fwhm, lp2=peak+2.*fwhm;
  for(int i=0;i<NPER;i++){
    lper[i]  = lp1 + (lp2-lp1)*i/(NPER-1);
    pprob[i] = exp(-0.5*sqr((lper[i]-peak)/(fwhm/2.3548)));
  }

  // Compute chi**2 threshold equivalent to the false alarm probability.

  float thresh = tchi(lfprob,ndat-1);

  // Construct initial random population 

  bool ok;
  sequence pop[npop], desc[npop];
  
  for(int np=0; np < npop; np++){
    pop[np].set_nseq(ndat);
    pop[np].ran(NTIMES,seed);
    pop[np].order();
  }

  unsigned long int key[npop], j, j1, j2;
  double fit[npop], breed[npop], x, mfit, sum;

  // compute breeding probabilities 
  
  sum = 0.;
  for(int np=0; np < npop; np++){
    sum += 1.+select*(npop-1-2.*np)/(npop-1);
    breed[np] = sum; 
  }
  for(int np=0; np < npop; np++){
    breed[np] /= breed[npop-1];
  }

  graphics pt(device);
  pt.sch(1.5);
  pt.scf(2);
  pt.sci(4);
  float x1 = time[0];
  float x2 = time[NTIMES-1];
  float range = x2-x1;
  x1 -= range/20.;
  x2 += range/20.;
  pt.env(x1,x2,0.,ngen+1,0,0);
  pt.sci(2);
  pt.lab("Time (days)","Generation");
  pt.sci(1);
  pt.sch(1.);
  
  float kfac = pow(TWOPI*(G*MSUN)*(m1+m2)/DAY,0.333333333333)*m2/(m1+m2)/1.e3;
  float wgt, k, v[ndat], w[ndat], e[ndat], n[ndat];
  double period, ph, lp, detm, norm, detp, t[ndat], sumw = 0.;
  for(int nd=0; nd<ndat; nd++){
    e[nd] = sigma;
    w[nd] = 1./sqr(sigma);
    sumw += w[nd];
  }
  for(int ng=0; ng < ngen; ng++){ // loop over generations

    // Evaluate fitness parameter for each sequence
    // We integrate over a grid in period and phase
    // and attempt some averaging over noise.

    for(int np=0; np < npop; np++){
      for(int nd=0; nd<ndat; nd++){
	t[nd] = time[pop[np][nd]];
      }
      detm = norm = 0.;
      for(int nper=0; nper < nperiod; nper++){ // loop over period
	lp  = lper[0] + (lper[NPER-1]-lper[0])*nper/(nperiod-1);
	j   = locate(lper,NPER,lp);
	if(j > 0 && j < NPER){
	  wgt = ((lp-lper[j-1])*pprob[j]+(lper[j]-lp)*pprob[j-1])/
	    (lper[j]-lper[j-1]);
	  if(nper == 0 || nper == nperiod-1) wgt /= 2;
	  period = pow(10.,lp);
	  
	  // maximum orbital velocity (km/s) 
	  
	  k = kfac/pow(period,0.333333333333);
	  
	  detp = 0.;
	  for(int nph=0; nph < nphase; nph++){ // loop over phases
	    if(ntrial){
	      ph = nph/(double)nphase;
	    }else{
	      ph = nph/(double)nphase/2.;
	    }
	
	    detp += dpb(period, ph, k, ntrial, ndat, t, w, e, v, n, sumw, 
		    thresh, seed);
	  }
	  detp /= nphase;
	  detm += wgt*detp;
	  norm += wgt;
	}
      }
      fit[np] = 1.-detm/norm;
    }

    // rank by fitness 

    heaprank(fit,key,npop);
    mfit = fit[key[npop/2]];
    if(npop % 2 == 0) mfit = (mfit+fit[key[npop/2+1]])/2;

    // Progress report

    std::cout << "Generation " << ng+1 << ", median = " << 1-mfit << ", best = "
	 << 1.-fit[key[0]] << ", sequence [" << pop[key[0]] << "]" << std::endl;
    
    // plot the top pfrac, with an attempt to reduce the final plot size

    for(int nt=0; nt < NTIMES; nt++){
      plot[nt] = false;
    }
    for(int np=0; np < (int)(pfrac*npop); np++){
      for(int nd=0; nd < ndat; nd++){
	plot[pop[key[np]][nd]] = true;
      }
    }
    for(int nt=0; nt < NTIMES; nt++){
      if(plot[nt]) pt.pt1(time[nt],ng+1,1);
    }

    // breed a new population, transferring best unchanged

    desc[0] = pop[key[0]];
    for(int np=1; np < npop; np++){

      // generate new individuals

      if(ran2(seed) < cross){
	do{
	  j1 = 0;
	  while(j1 == 0 || j1 >= npop){
	    j1 = locate(breed,npop,(double)ran2(seed));
	  }
	  j2 = 0;
	  while(j2 == 0 || j2 >= npop){
	    j2 = locate(breed,npop,(double)ran2(seed));
	  }
	  desc[np] = splice(pop[key[j1]],pop[key[j2]],seed);
	}while(!desc[np].differ());
      }else{
	j = 0;
	while(j == 0 || j >= npop){
	  j = locate(breed,npop,(double)ran2(seed));
	}
	desc[np] = pop[key[j]];
      }

      // add mutations
      
      if(ran2(seed) < mutate1) desc[np].rmutate(NTIMES,seed);
      if(ran2(seed) < mutate2) desc[np].lmutate(NTIMES,mdiff,seed);
      if(ran2(seed) < mutate3) desc[np].ran(NTIMES,seed);

    }

    // copy over and order

    for(int np=0; np < npop; np++){
      pop[np] = desc[np];
      pop[np].order();
    }

  }
  pt.clos();
  exit(EXIT_SUCCESS);
}








