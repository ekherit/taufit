/*
 * =====================================================================================
 *
 *       Filename:  vplepton.h
 *
 *    Description:  lepton vacuum polarization (changing tau)
 *
 *        Version:  1.0
 *        Created:  11.03.2012 01:41:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#ifndef IBN_VP_LEPTON_H
#define IBN_VP_LEPTON_H
#include <cmath>
#include "Const.h"
namespace vp
{
  inline double phi(double x)
  {
    if(x>1) return 0;
    if(x<0) return 0;
    return (2.+x)/6*sqrt(1-x);
  }

  inline double f(double x)
  {
    double aux=-5./9.-x/3.;
    double y=sqrt(fabs(1-x));
    if(x>1)
    {
      aux+=(2.+x)/3.*y*atan(1./y);
    }
    else
    {
      aux+=(2.+x)/6.*y*log((1+y)/fabs(1-y));    
    } 
    return aux; 
  }

  inline void lepton(double S,double *Re,double *Im, double mtau)
  {
    double Re1=0,Im1=0,Re2=0,Im2=0;
    double xmu,xtau,L;

    double xme = 4*ME*ME/S;
    Re1+=f(xme);
    Im1+=-M_PI*phi(xme);


    xmu  = 4.*MMU*MMU/S;
    Re1+=f(xmu);
    Im1+=-M_PI*phi(xmu);

    xtau = 4.*mtau*mtau/S;
    Re1+=f(xtau);
    Im1+=-M_PI*phi(xtau);


    *Re= ALPHAPI*Re1 + ALPHAPI*ALPHAPI*Re2;
    *Im= ALPHAPI*Im1 + ALPHAPI*ALPHAPI*Im2;
  }
}
#endif
