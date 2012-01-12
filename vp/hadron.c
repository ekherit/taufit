#include "inc.h"

void VP_Hadron(double S,double *Re,double *Im)
 {
  int i;


  *Re=0;
  *Im=0;
  if(S>Max2) return;

//-------------------------------------------------- narrow resonances
   for(i=0;i<nr;i++)  
    {
     if(i!=2) *Re+=bw_int(S,Mr[i],Gr[i],Beer[i]);
    } 

//-------------------------------------------------- Континуум

   FillCoef2_Tab(TabR,S);
   *Re+=GetIntegral_Tab(TabR);

//-------------------------------------------------- Континуум выше 10 ГэВ

//--------------------------------------------------

  *Re=-S/4./pi/pi/alpha*(*Re);
  *Im=-S/4./pi/alpha*VP_R(S)*eemm(S);
  
 }

//----------------------------- Счет ошибок и корреляций

void VP_Hadron_Err(double S, double *D_Re, double *D_Im, double *D_ReIm)
 {
  double C_Re, C_Re2, C_Im, C_Im2;
  *D_Re   =0;
  *D_Im   =0;
  *D_ReIm =0;

  if(S<(treshold2+deltaS)||(S>Max2))
   {
     return;
   }

   
   FillCoef2_Tab(TabR,S); 
   *D_Re+=GetDisp_Tab(TabR);

   Mark_Tab(TabR,S);

  *D_Im+=GetDispS_Tab(TabR);

   *D_ReIm+=GetDispRI_Tab(TabR);


  C_Re=-S/4./pi/pi/alpha;
  C_Im=-S/4./pi/alpha*eemm(S)   /S/eemm(S);

  C_Re2=C_Re*C_Re;
  C_Im2=C_Im*C_Im;

  *D_Re=C_Re2*(*D_Re);
  *D_Im=C_Im2*(*D_Im);
  *D_ReIm=S/(16*pi*pi*pi*alpha*alpha)*(*D_ReIm);


 } 
