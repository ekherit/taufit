/*
 * =====================================================================================
 *
 *       Filename:  draw.h
 *
 *    Description:  Draw functions
 *
 *        Version:  1.0
 *        Created:  22.02.2013 16:27:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#ifndef IBN_TAUFIT_DRAW_H
#define IBN_TAUFIT_DRAW_H
#include <TGraphErrors.h>
#include <TGraph.h>
#include "TF1.h"

#include "ScanPoint.h"
#include <algorithm>

TGraphErrors * DataGraph(const list<ScanPoint_t> & SPL, string title="r_{i} \\varepsilon \\sigma(e^{+}e^{-}\\rightarrow \\tau^{+}\\tau^{-}) + \\sigma_{B}")
{
  TGraphErrors * gr = new TGraphErrors;
  int i=0;
  for(auto & sp : SPL)
  {
    gr->SetPoint(i, sp.energy.value, sp.Ntt/sp.luminosity.value/sp.effcor);
    gr->SetPointError(i, sp.energy.error, sqrt(sp.Ntt)/sp.luminosity.value/sp.effcor);
    i++;
  }
	gr->SetTitle(title.c_str());
  gr->SetLineWidth(2);
  gr->SetMarkerSize(1);
  gr->SetMarkerStyle(20);
	return gr;
}

TGraph * SigmaGraph(const std::list<ScanPoint_t> & SPL,  double MTAU, double  EFFECT, double BG, int NN )
{
  BG = fabs(BG);
  EFFECT = fabs(EFFECT);
  std::vector<double> S; //cross section
  std::vector<double> E; //energy
  S.reserve(NN);
  E.reserve(NN);
	double lum=1, efcor = 1,spread = 1;
  auto minmax = std::minmax_element(begin(SPL), end(SPL), [](const ScanPoint_t & p1, const ScanPoint_t & p2) {return p1.energy.value < p2.energy.value; });
  double EEmin =  minmax.first->energy.value;
  double EEmax =  minmax.second->energy.value;
  EEmin-=5;
  EEmax+=5;
  //calculate parameters for energy spread
  TF1 f("fun_spread","[0]*x*x");
  f.SetParameter(0,0);
  TGraph g;
  int i=0;
  for(auto & sp :SPL)
  {
    g.SetPoint(i, sp.energy.value, sp.energy_spread.value);
    i++;
  }
  g.Fit("fun_spread","goff");
  double Spar=f.GetParameter(0);
  for(int i=0; i< NN; i++)
  {
		E[i] = (EEmin + (EEmax - EEmin)/(NN-1)*i);
    spread = Spar*sq(E[i]);
		S[i]=(EFFECT*sigma_total(2*E[i],spread,MTAU,0.0001) + BG);
  }
	return new TGraph(NN,&E[0],&S[0]);
}

void draw_lum(const std::list<ScanPoint_t> & SPL)
{
  TGraphErrors * glum = new TGraphErrors;
  int i=0;
  for(auto & sp : SPL)
  {
    double r = double(sp.Nee)/double(sp.Ngg);
    glum->SetPoint(i,sp.energy.value, r);
    double dr = r*sqrt(1./sp.Ngg + 1./sp.Nee);
    glum->SetPointError(i, sp.energy.error, dr);
    i++;
  }
  TCanvas * clum = new TCanvas;
  glum->SetMarkerStyle(21);
  glum->Draw("ap");
}

extern string OUTPUT_FILE;
void draw_fitresult(TauMassFitter & fitter)
{
	TCanvas * sigma_c = new TCanvas("sigma", "sigma", 1.2*640,1.2*480); 
	sigma_c->SetGrid();
  TMultiGraph * mg = new TMultiGraph;
	TGraph * fit_gr = SigmaGraph(fitter.GetData(), fitter.M.value,fitter.EPS.value,fitter.BG.value,200);
  TGraph * data_gr = DataGraph(fitter.GetData(),"tit");
  fit_gr->SetLineWidth(2);
  fit_gr->SetLineColor(kRed);
  mg->Add(fit_gr,"c");
  mg->Add(data_gr,"p");
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("E, MeV");
  mg->GetYaxis()->SetTitle("#sigma_{obs}, pb");

  double x,y;
  x=fitter.M.value-5;
  y=15;
  std::vector <double> v(&data_gr->GetY()[0], &data_gr->GetY()[data_gr->GetN()]);
  std::sort(v.begin(),v.end());
  y=v.back();

  char texbuf[1024];
  sprintf(texbuf,"M_{#tau} = %8.3f^{%+4.2f}_{%+4.2f} MeV", fitter.M.value, fitter.errDM.first,fitter.errDM.second);
  TLatex * Mtex = new TLatex(x,y,texbuf);
  Mtex->Draw();
  sprintf(texbuf,"#varepsilon = %3.1f #pm %3.1f %%", fitter.EPS.value*100, fitter.EPS.error*100);
  TLatex * EPStex = new TLatex(x,y*0.9,texbuf);
  EPStex->Draw();
  sprintf(texbuf,"#sigma_{BG} = %3.1f^{%+4.2f}_{%+4.2f} pb", fitter.BG.value, fitter.errBG.first,fitter.errBG.second);
  TLatex *BGtex = new TLatex(x,y*0.8,texbuf);
  BGtex->Draw();

  sprintf(texbuf,"M_{#tau}-M_{PDG} = %5.3f #pm %4.2f MeV ",fitter.DM.value , sqrt(fitter.DM.error*fitter.DM.error + DMTAU_PDG*DMTAU_PDG));
  TLatex * DMtex = new TLatex(1783,y*0.3,texbuf);
  DMtex->Draw();
  
  sprintf(texbuf,"#chi^{2}/ndf = %2.1f/%d", fitter.CHI2, fitter.NDF);
  TLatex * chi2tex = new TLatex(1783,y*0.2,texbuf);
  chi2tex->Draw();

  sprintf(texbuf,"P(#chi^{2},ndf) = %3.1f", TMath::Prob(fitter.CHI2,fitter.NDF));
  TLatex * probtex = new TLatex(1783,y*0.1,texbuf);
  probtex->Draw();


  std::string pdf_filename = OUTPUT_FILE+".pdf";
  std::string root_filename = OUTPUT_FILE+".root";
  gPad->SaveAs(pdf_filename.c_str());
  gPad->SaveAs(root_filename.c_str());

  draw_lum(fitter.GetData());
}
#endif
