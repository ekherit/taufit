#include"inc.h"

//---------------------------------------------------------------------- Point

Point CreatePoint(double x,double y,double ey)
 {
  Point aux;
  aux=malloc(sizeof(tpoint)); 
  if(!aux) ERROR("malloc");

  aux->x=x;
  aux->y=y;
  aux->ey=ey;
  aux->c=0;
  aux->m=0;
  aux->next=NULL;
  return aux;
 }

void ExchangePoints(Point p1, Point p2)
 {
  double aux;
  aux=p1->x;   p1->x  =p2->x;     p2->x  =aux;
  aux=p1->y;   p1->y  =p2->y;     p2->y  =aux;
  aux=p1->ey;  p1->ey =p2->ey;    p2->ey =aux;
  aux=p1->m;   p1->m  =p2->m;     p2->m  =aux;
 }
 
void LinkPoint(Point v,Point p)
 {
  Point aux,old;
  double a,b,ea,eb;

  if(p->x<v->x)
   {
    ExchangePoints(v,p);
    p->next=v->next;
    v->next=p;
    return;
   }  

  aux=v;
  old=v;
  while(aux)
   {
    if(p->x==aux->x)
     {
      a=aux->y; ea=aux->ey;
      b=p->y;   eb=p->ey;
      aux->y=( a/(ea*ea) +b/(eb*eb) )/( 1./(ea*ea) + 1./(eb*eb) );
      aux->ey=1./sqrt( 1./(ea*ea) + 1./(eb*eb) );
      free(p);
      return;
     }

    if( ( (p->x)>(old->x) ) && ( (p->x)<(aux->x) ) )
     {
      p->next=aux;
      old->next=p;
      return;
     } 
    old=aux;
    aux=aux->next; 
   }
   old->next=p;
 }

void AddPoint(Point v,double x,double y,double ey)
 {
  Point aux;
  aux=CreatePoint(x,y,ey);
  LinkPoint(v,aux);
 }
 
void PrintChain(Point v)
 {
  Point aux;
  aux=v;
  while(aux)
   {
    printf("%lf\t%lf\t%lf\t%lf\t%lf\n",aux->x,aux->y,aux->ey,aux->c,aux->m);
    aux=aux->next;
   }
 }

void ClearCoef_Chain(Point v)
 {
  Point aux;
  aux=v;
  while(aux)
   {
    aux->c=0;
    aux=aux->next;
   }
 }

void ClearMark_Chain(Point v)
 {
  Point aux;
  aux=v;
  while(aux)
   {
    aux->m=0;
    aux=aux->next;
   }
 }


void NegCoef_Chain(Point v)
 {
  Point aux;
  aux=v;
  while(aux)
   {
    aux->c=-(aux->c);
    aux=aux->next;
   }
 }

void MulCoef_Chain(Point v,double k)
 {
  Point aux;
  aux=v;
  while(aux)
   {
    aux->c=k*(aux->c);
    aux=aux->next;
   }
 }


double Interpolate_Chain(Point v,double x)
 {
  Point aux;
  double xx,yy;

  aux=v;
  if(x<aux->x) return 0;
  while(aux->next)
   {
    if( (x>aux->x) && (x<=aux->next->x) )
     {
      xx=x-aux->x;
      yy=aux->y + xx*(aux->next->y - aux->y)/(aux->next->x - aux->x); 
      return yy;
     }
    aux=aux->next;
   }
  return 0; 
 }


void Mark_Chain(Point v,double x)
 {
  Point aux;

  aux=v;
  if(x<aux->x) return ;
  while(aux->next)
   {
    if( (x>aux->x) && (x<=aux->next->x) )
     {
      aux->m       = (aux->next->x - x) / (aux->next->x - aux->x)   ;
      aux->next->m =    (x - aux->x)    / (aux->next->x - aux->x)   ;
     }
    aux=aux->next;
   }
 }
 
 
void FreeChain(Point v)
 {
  Point aux1,aux2;
  aux1=v;
  while(aux1)
   {
    aux2=aux1->next;
    free(aux1);
    aux1=aux2;
   }
 }




//-------------------------------------------------------------------- Tab

Tab CreateTab(Point p)
 {
  Tab aux;
  aux=malloc(sizeof(tTab)); 
  if(!aux) ERROR("malloc");

  aux->chain=p;
  aux->next=NULL;
  return aux;
 }

void AddToTab(Tab t,Point p)
 {
  Tab aux,new;
  new=CreateTab(p);
  aux=t;
  while(aux->next) aux=aux->next;
  aux->next=new;
 }

void FreeTab(Tab t)
 {
  Tab aux1,aux2;
  aux1=t;
  while(aux1)
   { 
    aux2=aux1->next;
    FreeChain(aux1->chain);
    free(aux1);
    aux1=aux2;
   }
 }

void PrintTab(Tab t)
 {
  Tab aux;
  aux=t;
  printf("Tab  -----------------------------------------------------------------\n");
  while(aux)
   {
    printf("---- Chain-----------------------------------------------------------\n");
    PrintChain(aux->chain);
    aux=aux->next;
   }
 } 
 
 double Interpolate_Tab(Tab t,double x)
 {
  Tab aux;
  double F=0;
  
  aux=t;
  while(aux)
   {
    F+=Interpolate_Chain(aux->chain,x);
    aux=aux->next;
   }
  return F; 
 }

void ClearMark_Tab(Tab t)
 {
  Tab aux;
  
  aux=t;
  while(aux)
   {
    ClearMark_Chain(aux->chain);
    aux=aux->next;
   }
 }

void Mark_Tab(Tab t,double x)
 {
  Tab aux;
  
  aux=t;
  ClearMark_Tab(aux);
  while(aux)
   {
    Mark_Chain(aux->chain,x);
    aux=aux->next;
   }
 }

