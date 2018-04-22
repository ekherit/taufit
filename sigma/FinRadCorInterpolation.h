/*
 * =====================================================================================
 *
 *       Filename:  FinRadCorInterpolation.h
 *
 *    Description:  Final state radiative correction interpolation
 *
 *        Version:  1.0
 *        Created:  21.04.2018 14:34:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once
#include "FinRadCor.h"
#include "interpolate.h"

/* this part creates tables for final state radiative corrections */
class FSRC //with linear correction in the table
{
  double precision;
  /* 
   * 	maximum and minimal meaning of tau velocity
   * 	belong this value direct radiation corrections will be performed
   */
  double vmax, vmin; 
  double dv; //elementary step in velocity (it depends on precision )
  vector <double> rad_cor_array;
  double vel(unsigned i )	{	return (vmin  + i*dv); }
  bool isinit;
  public:
  FSRC(void){};
  FSRC( double p ,double  mx = 0.4, double mn = 0 )
  {
    vmin = mn;
    vmax = mx;
    precision = p;
    isinit=false;
  }	

  void Init(void)
  {
    cout << "Init FSRC table v=("<<vmin<<","<<vmax<<") ..." << flush ;
    //vmin = mn;
    //vmax = mx;
    if(vmin <  0 )
    { 
      cerr << "error: FinRadCor class: vmin below zero." << endl;
      exit(1);
    }
    if(vmax <= vmin)
    {
      cerr << "error: FinRadCor class: vmax <= vmin in Fin." << endl;
      exit(1);
    }
    if(precision==0)
    {
      cerr << "error: FinRadCor class: precision equal zero." << endl;
      exit(1);
    }
    unsigned N = unsigned (1./precision);
    double v;
    dv = (vmax-vmin)/N;
    rad_cor_array.resize(N+1);
    for ( unsigned i = 0 ; i< rad_cor_array.size() ; i++)
    {
      v= vel(i);
      rad_cor_array[i] = Fr(v);
      //cout << i << ": " << v << ": " << rad_cor_array[i] << endl;

    }
    isinit =true;
    std::cout << " OK\n";
  }

  double operator()(double v)
  {
    if(!isinit) Init();
    if(v >= vmin && v <=vmax)
    {
      double x=(v-vmin)/dv;
      size_t i=size_t(x);
      double r;
      if(i<rad_cor_array.size()-1)
      {
        r=rad_cor_array[i] + (x - double(i))*(rad_cor_array[i+1]-rad_cor_array[i]);
      }
      else
      {
        r=rad_cor_array[i];
      }
      return  r;
    } else
    {
      return Fr(v);
    }
  }
};


class FsrcSpline3
{
  TSpline3 * spline;
  public:
    FsrcSpline3(double precision=1e-5, double vmax=0.99)
    {
      auto f = [](double v) { return Fr(v)-1.0;};
      spline = interpolate(f,0, vmax, precision);
    }
    ~FsrcSpline3(void) { delete spline; }
    double operator()(double v) 
    {
      return spline->Eval(v) + 1.0;
    }
};

static FsrcSpline3 fsrc;
