#include <cmath>
#include "trm_subs.h"
#include "trm_constants.h"
#include "trm_buffer2d.h"
#include "trm_rvanal.h"

/**
 * Computes covariances as an estimate of parameter uncertainties. The covariances
 * are returned via a structure defined in rvanal.h
 * The function carries out the linear minimisation but \b assumes
 * that the frequency has been optimised already.
 * \param data the data
 * \param freq the best-fit frequncy
 * \return The 10 covariances are returned in a single structure. In an obvious
 * notation they are vcc, vca, vcp, vct, vaa, vap, vat, vpp, vpt, vtt.
 */

Rvanal::sinerr Rvanal::sinfit_errors(const Subs::Buffer1D<Subs::rv>& data, double freq){

  size_t i, j, n = data.size();

  if(n < 4) throw std::string("Fewer than 4 points in sinfit_errors");

  float con, amp;
  double t0;
  sinfit_chisq(data, freq, con, amp, t0);

  // OK, off we go ... slightly pathetically cx is named
  // to avoid the c in the parameter list

  double cx, s;
  double wgt, theta;
  double ccc, cca, ccp, ccx, caa, cap, cax, cpp, cpx, cxx;
  double omega = Constants::TWOPI*freq;
  double dp, dt, dy;

  ccc=cca=ccp=ccx=caa=cap=cax=cpp=cpx=cxx=0;
  for(j=0; j<n; j++){
    if(data[j].z>0.){
      wgt   = 1./Subs::sqr(data[j].z);
      theta = omega*(data[j].x-t0);
      cx    = cos(theta);
      s     = sin(theta);
      dp    = -amp*theta*cx*freq; // derivative wrt period
      dt    = -amp*omega*cx;        // derivative wrt t0
      dy    = data[j].y-con-amp*s;    // data - fit

      // Second derivatives of chi**2 / 2.
      // each case = sum over weight*(fd1*fd2-dy*sd12)
      // where df1, fd2 are first derivatives of model wrt 
      // parameters 1 and 2 sd12 = second derivative of 
      // model wrt 1, 2

      ccc  +=  wgt;               
      cca  +=  wgt*s;
      ccp  +=  wgt*dp;
      ccx  +=  wgt*dt;
      caa  +=  wgt*s*s;
      cap  +=  wgt*(s*dp+dy*theta*cx*freq);
      cax  +=  wgt*(s*dt+dy*omega*cx);
      cpp  +=  wgt*(dp*dp+dy*(2.*dp*freq+amp*s*Subs::sqr(freq*theta)));
      cpx  +=  wgt*(dp*dt+dy*(dt+amp*theta*omega*s)*freq);
      cxx  +=  wgt*(dt*dt+dy*omega*omega*amp*s);
    }      
  }

  Subs::Buffer1D<size_t> indx(4);
  Subs::Buffer2D<double> m(4,4), ainv(4,4);
  m[0][0] = ccc;
  m[1][0] = m[0][1] = cca;
  m[2][0] = m[0][2] = ccp;
  m[3][0] = m[0][3] = ccx;
  m[1][1] = caa;
  m[2][1] = m[1][2] = cap;
  m[3][1] = m[1][3] = cax;
  m[2][2] = cpp;
  m[3][2] = m[2][3] = cpx;
  m[3][3] = cxx;

  double d;
  Subs::ludcmp(m,indx,d);

  Subs::Buffer1D<double> col(4);
  for(j=0;j<4;j++){
    for(i=0;i<4;i++) col[i] = 0.;
    col[j] = 1.;
    Subs::lubksb(m,indx,col);
    for(i=0;i<4;i++) ainv[i][j] = col[i];
  }
  sinerr err;
  err.vcc = ainv[0][0];
  err.vca = ainv[0][1];
  err.vcp = ainv[0][2];
  err.vct = ainv[0][3];
  err.vaa = ainv[1][1];
  err.vap = ainv[1][2];
  err.vat = ainv[1][3];
  err.vpp = ainv[2][2];
  err.vpt = ainv[2][3];
  err.vtt = ainv[3][3];
  return err;
}




