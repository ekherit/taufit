/*
 * =====================================================================================
 *
 *       Filename:  SigmaInterpolation.h
 *
 *    Description:  This is wrapper for sigma interpolation
 *
 *        Version:  1.0
 *        Created:  21.04.2018 15:46:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#include "Sigma.h"

#include "TGraph2D.h"
#include <memory>

class SigmaInterpolated
{
  //map of the energy spread
  std::map<double, unique_ptr<TGraph2D>> graph;
  double MMIN, MMAX, WMIN, WMAX;
  public:
  SigmaInterpolated(double mmin=MTAU-1, double mmax=MTAU+1, double Wmin=2*MTAU-20, double Wmax= 2*MTAU+30)
  {
    MMIN = mmin;
    MMAX = mmax;
    WMIN = Wmin;
    WMAX = Wmax;
  }

  unique_ptr<TGraph2D>  Init(double delta, double prec=1e-10, long Nmax=20)
  {
    TGraph2D * g = new TGraph2D;
    long point = 0;
    std::cout << "Making TGraph2D for delta = " << delta << endl;
    for(long i = 0;i < Nmax; i++)
    {
      double M = MMIN + (MMAX-MMIN)/Nmax*i;
      for(long j=0;j<Nmax; j++)
      {
        double W = WMIN + (WMAX-WMIN)/Nmax*j;
        double sigma = sigma_total(W, delta,M,1e-10);
        cout << W/2-MTAU << "  " << M-MTAU << "   " << sigma << endl;
        g->SetPoint(point,M,W,sigma);
        point++;
      }
    }
    return unique_ptr<TGraph2D>(g);
  }

  double operator()(double W, double delta, double mt, double prec=1e-6)
  {
    auto it = graph.find(delta);
    if(it == graph.end())
    {
      it = graph.emplace(delta, Init(delta, prec)).first;
    }
    return it->second->Interpolate(W,mt);
  }
};
