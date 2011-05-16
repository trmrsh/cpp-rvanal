#include <cstdlib>
#include <cmath>
#include <string>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_buffer2d.h"
#include "trm_rvanal.h"

/**
 * Given a set of data that can be fitted with a constant plus a sinusoid,
 * this function computes a number representing the probability that the frequency
 * of the sinusoid is a certain value. This is a Bayesian calculation which can
 * almost be done analytically except for degeneracy in the parameters. To cope
 * with this, a little Monte Carlo integration is needed.
 *
 * The calculation involves integrating over a 3D space spanned by the systemic velocity
 * and the coefficients of the cosine and sine components of the sinusoid. A prior over
 * this 3D space needs to be defined. For simplicity I take it to be uniform over a cylinder
 * with axis along the systemic velocity axis extending from -gmax to gmax and radius kmax.
 * The integrand is a 3D gaussian. If wholly enclosed within the cylinder there is an
 * analytic result for the value of the integral. If towards the edge, or if there is
 * degeneracy, the analytic result is an upper limit. A better value can then be obtained
 * with a bit of Monte Carlo integration. One has to be careful here however because it
 * is easy to come up with a very inefficient scheme since in many cases the integrand is only
 * significant over a tiny region of the available volume, so just randomly picking points
 * within the cylinder could be disastrous. Instead a weighted Monte Carlo is carried out 
 * in which points are generated in the most probable volume, but only contribute if they are
 * within the region of integration. This requires computation of the point of maximum probability
 * along with the eigen vectors and values of the three by three matrix describing the chi**2
 * variations near it (probability ~ exp(-chi**2/2). The user must choose how many points to use.
 * Very often every single one will lie in the volume of integration and the end result will be
 * the same as the analyticc estimate, but sometimes it won't.
 *
 * Note that Rvanal::sinfit_chisq(const Subs::Buffer1D<Subs::rv>&, double, float&, float&, double&) is essentially
 * a much simpler version of this routine that just returns the exponent at the best fit. It is however
 * not exactly the same.
 *
 * \param data the set of times, velocities and errors
 * \param f the frequency of interest in cycles per unit time
 * \param gmax the maximum gamma velocity
 * \param kmax the maximum semi-amplitude
 * \param nmonte the number of samples to take, 0 will just attempt the analytic result only.
 * \param seed the seed integer for the monte carlo integration for nmonte > 0
 * \return The function return the natural log of the probability. This only has meaning when compared
 * against another value for a different frequency.
 * \exception The routine may throw a Rvanal:Rvanal_Error
 */

double Rvanal::bayes_prob(const Subs::Buffer1D<Subs::rv>& data, double f, double gmax, double kmax, int nmonte, Subs::INT4& seed){

    // Compute mid-point
    double tmin, tmax, tmid;
    tmin = tmax = data[0].x;
    for(int i=1; i<data.size(); i++){
	if(data[i].z > 0.){
	    if(data[i].x < tmin) tmin = data[i].x;
	    if(data[i].x > tmax) tmax = data[i].x;
	}
    }
    tmid = 0.5*(tmin+tmax);

    // Compute sums
    double sumw=0., sumwc=0., sumws=0., sumwcc=0., sumwcs=0., sumwss=0.;
    double sumwy=0., sumwcy=0., sumwsy=0., sumwyy=0.;
    double theta, sine, cosine, wgt;
    for(int i=1; i<data.size(); i++){
	if(data[i].z > 0.){
	    theta   = Constants::TWOPI*f*(data[i].x-tmid);
	    cosine  = cos(theta);
	    sine    = sin(theta);
	    wgt     = 1./Subs::sqr(data[i].z);

	    sumw   += wgt;
	    sumwc  += wgt*cosine;
	    sumws  += wgt*sine;
	    sumwcc += wgt*cosine*cosine;
	    sumwcs += wgt*cosine*sine;
	    sumwss += wgt*sine*sine;

	    sumwyy += wgt*data[i].y*data[i].y;
	    sumwy  += wgt*data[i].y;
	    sumwcy += wgt*cosine*data[i].y;
	    sumwsy += wgt*sine*data[i].y;
	} 
    }

    // Determine eigenvectors and eigenvalues
    Subs::Buffer2D<double> a(3,3);
    a[0][0] = sumw;
    a[0][1] = a[1][0] = sumwc;
    a[0][2] = a[2][0] = sumws;
    a[1][1] = sumwcc;
    a[1][2] = a[2][1] = sumwcs;
    a[2][2] = sumwss;

    Subs::Buffer1D<double> d(3);
    Subs::Buffer2D<double> v(3,3);
    int nrot;
    Subs::jacob(a,d,v,nrot);

    if(d[0] <= 0. || d[1] <= 0. || d[2] <= 0.) 
	throw Rvanal_Error("Rvanal::bayes_prob(const Subs::Buffer1D<Subs::rv>&, double, double, double, int, long int): one or more eigenvalues not positive");


    // Determine best fit (solution of normal equations)
    Subs::Buffer2D<double> b(3,1);
    Subs::Buffer1D<double> amin(3);
    b[0][0] = sumwy;
    b[1][0] = sumwcy;
    b[2][0] = sumwsy;
    a[0][1] = a[1][0];
    a[0][2] = a[2][0];
    a[1][2] = a[2][1];
    Subs::gaussj(a, b);
    for(int i=0; i<3; i++) amin[i] = b[i][0];

    // Determine value of constant offset. Chi**2 = con + (quadratic form)
    double con = sumwyy - (sumwy*amin[0]+sumwcy*amin[1]+sumwsy*amin[2]);

    // If no Monte Carlo wanted, can save a little bit of effort ... the second term
    // is essentially a Jacobian, the third is a volume of integration factor
    if(nmonte == 0) return -con/2. - log(d[0]*d[1]*d[2])/2. - log(gmax*kmax*kmax);

    // Now Monte Carlo
    Subs::Buffer1D<double> scale(3), lambda(3), posvec(3);
    for(int i=0; i<3; i++) scale[i] = 1./sqrt(d[i]);

    int nok = 0;
    for(int n=0; n<nmonte; n++){

	// Generate multipliers
	for(int i=0; i<3; i++) lambda[i] = scale[i]*Subs::gauss2(seed);

	// Add in to form position
	for(int i=0; i<3; i++){
	    posvec[i] = amin[i];
	    for(int j=0; j<3; j++)
		posvec[i] += lambda[i]*v[i][j];
	}

	// Now test to see if it is within the cylinder of integration
	if(posvec[0] > -gmax && posvec[0] < gmax && 
	   sqrt(Subs::sqr(posvec[1])+Subs::sqr(posvec[2])) < kmax) nok++;
    }

    // Return the log of the answer. The first term comes from the exponent,
    // the second is the analytic version of the integral (ignoring factor of Pi
    // with a correction factor from the Monte Carlo part. The last term accounts for
    // the volume of integration. 
    return -con/2. + log(double(nok)/nmonte*scale[0]*scale[1]*scale[2]) - log(gmax*kmax*kmax);
}










