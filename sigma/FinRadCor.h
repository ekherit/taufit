#pragma once
#ifndef IBN_FIN_RAD_COR_H
#define IBN_FIN_RAD_COR_H
#include "Const.h"
#include <iostream>
#include <complex>
#include <vector>
#include <ibn/integral.h>
/* Final states corrections... 
 * \sigma_1 = SIGMA_CONST b(3-b)/2 * Fc*Fr
 * 
 * M.B.Voloshin The onset of e^+e^- \rightarrow \tau^+ \tau^- at threshold revisitted
 * hep-ph/012207
 */

#include <TMath.h> 
#include <TSpline.h>

inline double Fc(double v); // Кулоновское взаимодейсвтие в конечном состоянии 
inline double vFc(double v); // Fc * v
inline double Fr(double v);  //радиационные поправки (петлевые и всякая хрень)

// Некоторые вспомогательные функции 
static const double PRECISION=1e-10;
inline double Li2(double x);

inline double S(double b);
inline double h(double v, double mt = MTAU);

using namespace std;

// Colomb intaraction correction in final states (Voloshin) 
inline double Fc(double v)
{
    double z=PIALPHA/v;
    return z/(1-exp(-z));
}

// Fc*v function
inline double vFc(double v)
{
    if(v==0) return PIALPHA;
    return PIALPHA/(1. - exp(-PIALPHA/v));
}



inline double Li2(double x)	
{
    return TMath::DiLog(x); //use standart function
}

inline double S(double v)
{
    if(v == 0) return -4;
    return 1./v * ( - PI2/2. + (1.+v*v)*( PI2/6.0 + log((1.+v)/2.0)*log((1.+v)/(1.-v)) + 
		2. * ( Li2((1.-v)/(1.+v)) + Li2((1.+v)/2.) - Li2((1.-v)/2.) -
		    2.*Li2(v) ) + Li2(v*v) ) +
	    log((1.+v)/(1.-v)) * ( 11./8.*(1.+v*v) - 3.*v + sq(v*v)/(3-v*v)/2.) +
	    6.*v*log((1.+v)/2.) - 4.*v*log(v) + 3./4.*v*(5.-3.*v*v)/(3.-v*v) );
}
inline double h(double v, double mt)
{
  static const double b0 =  GAMMA_E + 1./6.;
  static const double b1  =  log(2.0) - 5.0/6.0; 
  if(v==0) return log(mt*ALPHA/ME) +  b0;
  double z = PIALPHA/v;
  double e = exp(-z);
  double h = (1.0 - (1.0+z)*e)/(1. - e);
  double u = ALPHA + (v-ALPHA)*pow(e,1./M_PI);
  double b = b0 + (b1-b0)*pow(e,1./M_PI);
  return h*(log(u*mt/ME) + b);
}

inline double Fr(double v)
{
    double fr = 1. + ALPHAPI* (S(v) + 2./3. * h(v) ) ;
    //double fr = 1. + ALPHAPI* (2./3. * h(v) ) ;
    return fr;
}




inline double hsub(double lambda, double t, double x, double v, double z)
{
  complex<double> il(0,lambda);
  complex<double> a(t, z*x/v);
  double im = ( pow(1.0/t+1.0, il) * pow( a, il-1.0 ) * pow ( a+1.0, - il  - 1.0) ).imag();
  double _x = pow(x,-1.0);
  double _x2 = _x*_x;
  return im*(1.0 + 0.5*_x2)*sqrt(1.0 - _x2)*_x;
};

inline double hint2(double v)
{
  double a[2] = {0,1.0};
  double b[2] = {1000,1000};
  double lambda = MTAU*ALPHA/v*0.5;
  double z = ME/MTAU;
  double relerr;
  return -2*lambda*ibn::dgaus( 
      [&lambda,&v,&z](double * x)->double { return hsub(lambda, x[0],x[1],v,z); },
      2, 
      a,
      b,
      1e-10, 
      relerr);
}


#endif
