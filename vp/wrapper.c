#include"inc.h"

void vp_init()
 {
  VP_Init();
 }
 
void vp_done()
 {
  VP_Done();
 }
 
void vp_lepton(double *s, double *re, double *im)
 {
  double ss=*s;
  VP_Lepton(ss,re,im);
 }
void vp_hadron(double *s, double *re, double *im)
 {
  double ss=*s;
  VP_Hadron(ss,re,im);
 }

void vp_hadron_err(double *s, double *dre, double *dim, double *dreim)
 {
  double ss=*s;
  VP_Hadron_Err(ss,dre,dim,dreim);
 }


//------------------------------------- Для фортрана
void vp_init__()
 {
  vp_init();
 }
 
void vp_done__()
 {
  vp_done();
 }
 
void vp_lepton__(double *s, double *re, double *im)
 {
  vp_lepton(s,re,im);
 }

void vp_hadron__(double *s, double *re, double *im)
 {
  vp_hadron(s,re,im);
 }

void vp_hadron_err__(double *s, double *dre, double *dim, double *dreim)
 {
  vp_hadron_err(s,dre,dim,dreim);
 }
