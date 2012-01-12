#include "inc.h"

//------------------------------- R without  resonances 

double VP_Rcont(double s)
 {
  double aux;
  aux=0;


  aux+=Interpolate_Tab(TabR,s);
  aux=aux/s/eemm(s);
  return aux;
 }

//------------------------------- Contribution of resonances to R

double VP_RRes(double s)
 {
   int i;
   double aux=0;

   for(i=0;i<nr;i++) if(i!=2) aux+=bw(s,Mr[i],Gr[i],Beer[i]);

   aux=aux/eemm(s);
   return aux;
 }


//------------------------------- R

double VP_R(double s)
 {
   if(s<treshold2) return 0;
   return VP_Rcont(s)  + VP_RRes(s);
 } 

