/*
 * =====================================================================================
 *
 *       Filename:  draw.h
 *
 *    Description:  Draw data graph, fit function and fit result with chi square
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
#include <algorithm>

#include <fmt/format.h>

#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TAxis.h>

#include "ScanPoint.h"

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
  TCanvas * clum = new TCanvas("luminosity","luminosity");
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
  
  sprintf(texbuf,"#chi^{2}/ndf = %1.2f/%d", fitter.CHI2, fitter.NDF);
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

void draw_fitresult(TauMassFitter2 & fitter, double MSHIFT=0)
{
	TCanvas * sigma_c = new TCanvas("sigma", "sigma", 1.2*640,1.2*480); 
	sigma_c->SetGrid();
  TMultiGraph * mg = new TMultiGraph;
	TGraph * fit_gr = SigmaGraph(fitter.GetData(), MTAU+fitter("M").value, fitter("EPS").value, fitter("BG").value,200);
  TGraph * data_gr = DataGraph(fitter.GetData(),"tit");
  TGraphErrors * diff_graph = new TGraphErrors;
  for(int i=0;i<data_gr->GetN();i++)
  {
    double E = data_gr->GetX()[i];
    double dE = data_gr->GetErrorX(i);
    double sigma_data = data_gr->GetY()[i];
    double sigma_fit = fit_gr->Eval(E);
    double diff = sigma_fit-sigma_data;
    diff_graph->SetPoint(i, E, diff);
    diff_graph->SetPointError(i,dE, data_gr->GetErrorY(i));
  }
  for(int i=0;i<fit_gr->GetN();i++) fit_gr->GetX()[i]-=MSHIFT;
  for(int i=0;i<data_gr->GetN();i++) data_gr->GetX()[i]-=MSHIFT;
  for(int i=0;i<diff_graph->GetN();i++) diff_graph->GetX()[i]-=MSHIFT;
  fit_gr->SetLineWidth(2);
  fit_gr->SetLineColor(kRed);
  mg->Add(fit_gr,"c");
  mg->Add(data_gr,"p");
  mg->Draw("a");

  auto smart_xasis_title = [&](TAxis * axis)
  {
    if(MSHIFT == 0) axis->SetTitle("E, MeV");
    else if( MSHIFT = MTAU_PDG ) axis->SetTitle("E-M_{#tau}^{PDG}, MeV");
    else axis->SetTitle(fmt::format("E-{}, MeV", MSHIFT).c_str());
  };
  mg->GetYaxis()->SetTitle("#sigma_{obs}, pb");
  smart_xasis_title(mg->GetXaxis());

  //if(MSHIFT == 0) mg->GetXaxis()->SetTitle("E, MeV");
  //else if( MSHIFT = MTAU_PDG ) mg->GetXaxis()->SetTitle("E-M_{#tau}^{PDG}, MeV");
  //else mg->GetXaxis()->SetTitle(fmt::format("E-{}, MeV", MSHIFT).c_str());

  double x,y;
  x=fitter("M").value-5+MTAU-MSHIFT;
  y=15;
  std::vector <double> v(&data_gr->GetY()[0], &data_gr->GetY()[data_gr->GetN()]);
  std::sort(v.begin(),v.end());
  y=v.back();

  char texbuf[1024];
  if(fitter.isminos) sprintf(texbuf,"M_{#tau} = %8.3f_{ %+4.2f}^{%+4.2f} MeV", fitter("M").value+MTAU, fitter("M").min,fitter("M").max);
  else sprintf(texbuf,"M_{#tau} = %8.3f #pm %4.2f MeV", fitter("M").value+MTAU, fitter("M").error);
  TLatex * Mtex = new TLatex(x,y,texbuf);
  //TLatex * Mtex = new TLatex(x,y,fmt::format("M_{#tau} = %8.3f_{ %+4.2f}^{%+4.2f} MeV", fitter.M.value, fitter.errDM.first,fitter.errDM.second).c_str());
  Mtex->Draw();
  sprintf(texbuf,"#varepsilon = %3.1f #pm %3.1f %%", fitter("EPS").value*100, fitter("EPS").error*100);
  TLatex * EPStex = new TLatex(x,y*0.9,texbuf);
  EPStex->Draw();
  sprintf(texbuf,"#sigma_{BG} = %3.1f_{%+4.2f}^{%+4.2f} pb", fitter("BG").value, 0.0 ,fitter("BG").error);
  TLatex *BGtex = new TLatex(x,y*0.8,texbuf);
  BGtex->Draw();

  sprintf(texbuf,"M_{#tau}-M_{PDG} = %5.3f #pm %4.2f MeV ",fitter("M").value , sqrt(pow(fitter("M").error,2.0) + DMTAU_PDG*DMTAU_PDG));
  TLatex * DMtex = new TLatex(1783-MSHIFT,y*0.3,texbuf);
  DMtex->Draw();
  
  sprintf(texbuf,"#chi^{2}/ndf = %1.2f/%d", fitter.CHI2, fitter.NDF);
  TLatex * chi2tex = new TLatex(1783-MSHIFT,y*0.2,texbuf);
  chi2tex->Draw();

  sprintf(texbuf,"prob(#chi^{2},ndf) = %3.2f", TMath::Prob(fitter.CHI2,fitter.NDF));
  TLatex * probtex = new TLatex(1783-MSHIFT,y*0.1,texbuf);
  probtex->Draw();


  std::string pdf_filename = OUTPUT_FILE+".pdf";
  std::string root_filename = OUTPUT_FILE+".root";
  gPad->SaveAs(pdf_filename.c_str());
  gPad->SaveAs(root_filename.c_str());

  draw_lum(fitter.GetData());
  auto diff_canvas = new TCanvas("diff_canavs", "Fit and data difference");
  diff_graph->SetMarkerStyle(21);
  diff_graph->Draw("ap");
  smart_xasis_title(diff_graph->GetXaxis());
  diff_graph->GetYaxis()->SetTitle("#sigma_{fit}  - #sigma_{data}, pb");
}
#endif
