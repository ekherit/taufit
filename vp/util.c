#include"inc.h"

void ERROR(char *mess)
 {
  fprintf(stderr,"ERROR -> %s\n",mess);
  exit(1);
 }

//------------------------------- Cross section of ee -> mu mu

double eemm(double s)
 {
  double a;
  double beta=1.-4*mmu*mmu/s;
  
  if(beta<=0) return 0;
  if(beta>=1) return 0;
  beta=sqrt(beta);
  a=beta*(6-2*beta*beta);
  return a*pi*alpha*alpha/3./s;  
 }  

//------------------------------- Breit-Wigner cross section

double bw(double s,double m,double Gamma,double Bee)
 {
  return 12*pi*Bee*Gamma*Gamma/((s-m*m)*(s-m*m)+m*m*Gamma*Gamma);
 }

//------------------------------- Interference between two resonances

double interference( double s,
                     double m1,double Gamma1,double Bee1,double alpha1,
                     double m2,double Gamma2,double Bee2,double alpha2
                   )
 {
   double a   = (alpha1-alpha2)/180.*pi;
   double x   = (s-m1*m1)*(s-m2*m2)+m1*m2*Gamma1*Gamma2;
   double y   = m1*Gamma1*(s-m2*m2)-m2*Gamma2*(s-m1*m1);
   double aux = (x*cos(a)+y*sin(a))/(x*x+y*y);
          aux = aux*12.*pi*2.*sqrt(Bee1*Bee2)*Gamma1*Gamma2;
   return aux;                                                                 
 }                                                   

//------------------------------- int_{treshold}^{+inf}(bw(s')/(s'-s))ds'
 
double bw_int(double s,double m,double Gamma,double Bee)
 {
  double aux;
  aux=((treshold2-m*m)*(treshold2-m*m)+m*m*Gamma*Gamma)/(s-treshold2)/(s-treshold2); 
  aux=0.5*log(aux);
  aux-=(s-m*m)/m/Gamma*(pi/2+atan((m*m-treshold2)/m/Gamma));
  aux*=bw(s,m,Gamma,Bee);
  return aux;
 }

