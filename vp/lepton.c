#include "inc.h"

double phi(double x)
 {
  if(x>1) return 0;
  if(x<0) return 0;
  return (2.+x)/6*sqrt(1-x);
 }

double f(double x)
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
 
void VP_Lepton(double S,double *Re,double *Im)
 {
   double Re1=0,Im1=0,Re2=0,Im2=0;
   double xmu,xtau,L;
	/*
	
   L=log(fabs(S)/me/me);
    Re1+=L/3.-5./9.;
    if(S>0) Im1+=-pi/3.;
     Re2+=L/4.+zeta3-5./24;
     if(S>0) Im2+=-pi/4.;
		*/ 
		 
	 
	 double xme = 4*me*me/S;
	 Re1+=f(xme);
	 Im1+=-pi*phi(xme);
		
		 
   xmu  = 4.*mmu*mmu/S;
    Re1+=f(xmu);
    Im1+=-pi*phi(xmu);

   xtau = 4.*mtau*mtau/S;
    Re1+=f(xtau);
    Im1+=-pi*phi(xtau);
		
		
   *Re= alpi*Re1 + alpi*alpi*Re2;
   *Im= alpi*Im1 + alpi*alpi*Im2;
 }
 
 
