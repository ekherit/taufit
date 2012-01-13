/*
 * =====================================================================================
 *
 *       Filename:  draw_sigma.cpp
 *
 *    Description:  Draw tau prodaction cross section.
 *
 *        Version:  1.0
 *        Created:  12.01.2012 16:13:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#include <TROOT.h>
#include <TApplication.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("draw_tau_cross_section","Draw tau cross section", initfuncs);

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TAxis.h>


#include "sigma/Sigma.h"

// For measurement calculation speed.
#include <ibn/timer.h>

int main(int argc, char ** argv)
{
	TApplication theApp("tau", &argc, argv);
  bool issigma=true;
  bool isdilog=true;
  if(issigma)
  {
    double sigma_W = 1.5;//beam energy spread in MeV
    double Emin = 1770; //min beam energy in MeV
    double Emax = 1800; //max beam energy in MeV
    double dE = 0.1;    //energy step in MeV
    unsigned N = unsigned(fabs(Emax-Emin)/dE); //number of points in graph
    TGraph * sigma_total_g = new TGraph(N);
    for(unsigned i=0;i<N;++i)
    {
      double E = Emin + dE*i; //calculate beam energy
      double sigma = sigma_total(2*E, sigma_W,1777,1e-5);
      sigma_total_g->SetPoint(i,E,sigma);
    }
    sigma_total_g->SetTitle("Tau cross section");
    sigma_total_g->Draw("ac");
    sigma_total_g->GetXaxis()->SetTitle("E, MeV");
    sigma_total_g->GetYaxis()->SetTitle("#sigma, pb");
    sigma_total_g->SetLineWidth(2);
  }
  if(isdilog)
  {
    TCanvas * c = new TCanvas("dilog","dilog");
    double NN=100000.;
    TH1D * h = new TH1D("dilog","dilog",200,-20,0);
    for(double i = 0;i<NN;++i)
    {
      double x = i/NN;
      double my_dilog = Li2_int_reg(x);
      double cern_dilog=TMath::DiLog(x);
      double diff = fabs((-my_dilog+cern_dilog)/cern_dilog);
//      cout << x << " " << diff << endl;
      h->Fill(log(diff)/log(10));
    }
    h->Draw();
  }
  theApp.Run();
  return 0;
}


