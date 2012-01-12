#include "inc.h"

static double F(double x) 
 {
  double mx=fabs(x);
  if(mx<1.e-6) return 0;
  return x*log(mx);
 }


void FillCoef_Chain(Point v,double x0)
 {
  int i;
  double xa,xb,xc,FF,d_ba,d_cb,d_ca;
  Point pa,pb,pc;
  

  pb=v;       xb=pb->x;
  pc=v->next; xc=pc->x;    
   FF=F(x0-xb);
   d_cb=xc-xb;
    pb->c+= -0.5 + FF/d_cb - log(fabs(x0-xb));
    pc->c+=  0.5 - FF/d_cb;   

  pa=v;
  pb=pa->next;
  pc=pb->next;
  while(pc)  
   {
    xa=pa->x;  
    xb=pb->x;    
    xc=pc->x;  
     FF=F(x0-xb);
     d_ba=xb-xa; d_cb=xc-xb; d_ca=xc-xa;
      pa->c+= -0.5 - FF/d_ba;
      pb->c+=  FF*d_ca/d_cb/d_ba;
      pc->c+=  0.5 - FF/d_cb;   

    pa=pa->next;
    pb=pb->next;
    pc=pc->next;
   }

  xa=pa->x;  
  xb=pb->x;    
   FF=F(x0-xb);
   d_ba=xb-xa;
    pa->c+= -0.5 - FF/d_ba;
    pb->c+=  0.5 + FF/d_ba + log(fabs(x0-xb));   

 } 

void FillCoef2_Chain(Point v,double x0)
 {
  ClearCoef_Chain(v);
  FillCoef_Chain(v,0.);
  NegCoef_Chain(v);
  FillCoef_Chain(v,x0);
  MulCoef_Chain(v,1./x0);
 }

double GetIntegral_Chain(Point v)
 {
  Point aux;
  double I=0;
  
  aux=v;
  while(aux)
   {
    I+=aux->c*aux->y;
    aux=aux->next;
   } 
  return I;
 }
    
double GetDisp_Chain(Point v)
 {
  Point aux;
  double D=0;
  
  aux=v;
  while(aux)
   {
    D+=(aux->c*aux->ey)*(aux->c*aux->ey);
    aux=aux->next;
   } 
  return D;
 }


double GetDispS_Chain(Point v)
 {
  Point aux;
  double D=0;

  aux=v;
  while(aux)
   {
    D+=(aux->m*aux->ey)*(aux->m*aux->ey);
    aux=aux->next;
   }
  return D;
 }


double GetDispRI_Chain(Point v)
 {
  Point aux;
  double D=0;
  
  aux=v;
  while(aux)
   {
    D+=(aux->m*aux->ey)*(aux->c*aux->ey);
    aux=aux->next;
   }
  return D;
 }


void FillCoef_Tab(Tab t,double x0)
 {
  Tab aux;
  aux=t;
  while(aux)
   {
    FillCoef_Chain(aux->chain,x0);
    aux=aux->next;
   }
 }
 
void FillCoef2_Tab(Tab t,double x0)
 {
  Tab aux;
  aux=t;
  while(aux)
   {
    FillCoef2_Chain(aux->chain,x0);
    aux=aux->next;
   }
 }
 
double GetIntegral_Tab(Tab t)
 {
  Tab aux;
  double I=0;
  
  aux=t;
  while(aux)
   {
    I+=GetIntegral_Chain(aux->chain);
    aux=aux->next;
   } 
  return I; 
 }
 
double GetDisp_Tab(Tab t)
 {
  Tab aux;
  double D=0;
  
  aux=t;
  while(aux)
   {
    D+=GetDisp_Chain(aux->chain);
    aux=aux->next;
   } 
  return D;
 }

double GetDispS_Tab(Tab t)
 {
  Tab aux;
  double D=0;

  aux=t;
  while(aux)
   {
    D+=GetDispS_Chain(aux->chain);
    aux=aux->next;
   }
  return D;
 }
    

double GetDispRI_Tab(Tab t)
 {
  Tab aux;
  double D=0;

  aux=t;
  while(aux)
   {
    D+=GetDispRI_Chain(aux->chain);
    aux=aux->next;
   }
  return D;
 }

