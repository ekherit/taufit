/*
 * =====================================================================================
 *
 *       Filename:  draw_sigma.cpp
 *
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

#include <iostream>
#include <iomanip>
#include <map>
#include <chrono>

#include <boost/program_options.hpp>

#include <TROOT.h>
#include <TApplication.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("draw_tau_cross_section","Draw tau cross section", initfuncs);

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TF1.h>
#include <TSpline.h>

#include <fmt/printf.h>

#include "sigma/Sigma.h"
#include "interpolate.h"

// For measurement calculation speed.
#include <ibn/timer.h>
unsigned DEBUG=0;

TGraph * get_cross_section_graph(double sigma_W, double Emin=MTAU-6, double Emax=MTAU+16, double dE=0.2)
{
    unsigned N = unsigned(fabs(Emax-Emin)/dE); //number of points in graph
    TGraph * sigma_total_g = new TGraph(N);
    vp_init();
    for(unsigned i=0;i<N;++i)
    {
      double E = Emin + dE*i; //calculate beam energy
      double sigma = sigma_total(2*E, sigma_W, MTAU,1e-10);
      sigma_total_g->SetPoint(i,E-MTAU,sigma);
    }
    sigma_total_g->SetTitle("Tau cross section");
    sigma_total_g->Draw("ac");
    sigma_total_g->GetXaxis()->SetTitle("E, MeV");
    sigma_total_g->GetYaxis()->SetTitle("#sigma, pb");
    sigma_total_g->SetLineWidth(2);
    return sigma_total_g;
}

double sigma_approx(double * x, double *par)
{
  double E = x[0];
  double shift = par[0];
  double tau = par[1];
  double e =  (1.0+erf((E)/tau))*0.5;
  e = e*sqrt(1.0 + par[0]*e*E);
  //return e;
  //cout << E << "  " << e << endl;
  //e = 1;
  double cor = 0;
  for(int i=0;i<10; i++)
  {
    cor += par[2+i]*pow(E, i);
  }
  //std::cout << E << " " << e << "  " << cor << endl;
  return e*cor;
}

double fit_cross_section(void)
{
  TF1 * f = new TF1("f",sigma_approx,-6,+20,12);
  f->SetParameter(0,0);
  //f->FixParameter(0,0);
  f->SetParameter(1,1.47);
  f->SetParLimits(1,0,10);
  f->SetParameter(2,100);
  for(int i=1;i<10;i++) f->SetParameter(i,0);
  TGraph * g = get_cross_section_graph(1.47,MTAU-6, MTAU+16,0.01);
  TGraph * gs = get_cross_section_graph(1.47,MTAU-6, MTAU+16,0.1);
  TSpline3 * spline = new TSpline3("spline",gs->GetX(),gs->GetY(),gs->GetN());
  TCanvas * c =new TCanvas;
  c->Divide(1,2);
  c->cd(1);
  g->Draw("al");
  //f->Draw("same");
  //f->Draw("same");
  spline->SetLineColor(kRed);
  spline->Draw("same");
  //g->Fit("f");
  TGraph * gdiff = new TGraph(g->GetN());
  for(int i=0;i<g->GetN();i++)
  {
    //double y=f->Eval(g->GetX()[i]);
    double y=spline->Eval(g->GetX()[i]);
    double y0=g->GetY()[i];
    gdiff->SetPoint(i, g->GetX()[i], (y-y0)/y0);
  }
  c->cd(2);
  gdiff->Draw("al");
}



void draw_cross_section(void)
{
    double sigma_W = 1.5;//beam energy spread in MeV
    double Emin = 1750; //min beam energy in MeV
    double Emax = 1810; //max beam energy in MeV
    double dE = 0.1;    //energy step in MeV
    unsigned N = unsigned(fabs(Emax-Emin)/dE); //number of points in graph
    TGraph * sigma_total_g = new TGraph(N);
    vp_init();
    for(unsigned i=0;i<N;++i)
    {
      double E = Emin + dE*i; //calculate beam energy
      double sigma = sigma_total(2*E, sigma_W,MTAU,1e-10);
      sigma_total_g->SetPoint(i,E-MTAU,sigma);
    }
    sigma_total_g->SetTitle("Tau cross section");
    sigma_total_g->Draw("ac");
    sigma_total_g->GetXaxis()->SetTitle("E, MeV");
    sigma_total_g->GetYaxis()->SetTitle("#sigma, pb");
    sigma_total_g->SetLineWidth(2);
}

void draw_cross_section_diff(double dm, double sigma_W=1.47)
{
    double Emin = MTAU-6; //min beam energy in MeV
    double Emax = MTAU+16; //max beam energy in MeV
    double dE = 0.2;    //energy step in MeV
    unsigned N = unsigned(fabs(Emax-Emin)/dE); //number of points in graph
    TGraph * sigma_total_g = new TGraph(N);
    vp_init();
    for(unsigned i=0;i<N;++i)
    {
      double E = Emin + dE*i; //calculate beam energy
      double sigma0 = sigma_total(2*E, sigma_W,MTAU,1e-10);
      double sigma1 = sigma_total(2*(E+dm), sigma_W,MTAU+dm,1e-10);
      sigma_total_g->SetPoint(i,E-MTAU,(sigma1-sigma0)/sigma0);
    }
    sigma_total_g->SetTitle("Tau cross section");
    sigma_total_g->Draw("ac");
    sigma_total_g->GetXaxis()->SetTitle("E, MeV");
    sigma_total_g->GetYaxis()->SetTitle("#sigma, pb");
    sigma_total_g->SetLineWidth(2);
}


void draw_vp(void)
{
  vp_init();
  TCanvas * c = new TCanvas("vp","vp");
  TGraph * g =new TGraph();
  int i=0;
  double m=1777;
  for(double E = MTAU_PDG - 20; E<MTAU_PDG+20; E+=0.1)
  {
    double S = 4*E*E;
    double vpold = Vp(S,m);
    double vp = Vp2(S,m);
    cout << E << " " << vp << " " << (vp-Vp_old(S))/vp <<  " "  << (vp-VP1777(S,m))/vp<< endl;
    g->SetPoint(i,E,vp);
    i++;
  }
  g->Draw("al");
  g->SetLineWidth(2);

}

void draw_vp_hadron(void)
{
  vp_init();
  TCanvas * c = new TCanvas("vp","vp");
  c->Divide(2,0);
  TGraph * gre =new TGraph();
  TGraph * gim =new TGraph();
  int i=0;
  cout << setw(15) << "E, MeV" << setw(15) << "ReH" << setw(15) << "ImH" << endl;
  for(double E = MTAU_PDG - 50; E<MTAU_PDG+50; E+=0.1)
  {
    double S = 4*E*E;
    double reh = ReHadron(S);
    double imh = ImHadron(S);
    //cout << E << " " << vp << " " << (Vp(S)-Vp_old(S))/vp << endl;
    cout << setw(15) << E << setw(10) << reh << setw(15) << imh << endl;
    gre->SetPoint(i,E,reh);
    gim->SetPoint(i,E,imh);
    i++;
  }
  c->cd(1);
  gre->Draw("al");
  gre->SetLineWidth(2);
  c->cd(2);
  gim->Draw("al");
  gim->SetLineWidth(2);
}

int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  opt_desc.add_options()
    ("help,h","Print this help")
    ("spline", "Interpolate cross section")
    ("cross-section", "Draw cross section")
    ("cross-section-diff", "Draw cross section difference")
    ("cross-section-fit", "Fit the cross section")
    //("vp", po::value<std::string>()->default_value("reh"),"Vacuum polarization")
    ("vp", "Vacuum polarization")
    ("vph","Hadronic vacuum polarization")
    ("fsrc-cmp","Time comparizon of variose algorithm fsrc calculation")
    ("hint","h correction thru integral")
    ("h","h correction (my first interpolation)")
    ;
  po::positional_options_description pos;
  //pos.add("data",-1);
  po::variables_map opt; //options container
  try
  {
    po::store(po::command_line_parser(argc, argv).options(opt_desc).positional(pos).run(), opt);
    //std::ifstream config_file("fit.cfg");
    //po::store(po::parse_config_file(config_file,opt_desc,true), opt);
    po::notify(opt);
  } 
  catch (boost::program_options::error & po_error)
  {
    cerr << "WARGNING: configuration: "<< po_error.what() << endl;
  }

  if(opt.count("help"))
  {
    std::clog << opt_desc;
    return 0;
  }
	TApplication theApp("tau", &argc, argv);
  bool issigma=true;
  bool isdilog=false;
  if(opt.count("cross-section")) draw_cross_section();
  if(opt.count("cross-section-diff")) draw_cross_section_diff(1,1.47);
  if(opt.count("cross-section-fit")) fit_cross_section();
  if(opt.count("vp")) draw_vp();
  if(opt.count("vph")) draw_vp_hadron();
  if(opt.count("spline"))
  {
    auto f = [](double E) { return sigma_total(2*E, 1.47,MTAU,1e-10); };
    auto spline = interpolate(f,MTAU-5, MTAU+16, 1e-2);
    spline->SaveAs("spline.C");
    TGraph * g1 = new TGraph;
    TGraph * g2 = new TGraph;
    for(double E=MTAU-5; E<MTAU+16;E+=0.1)
    {
      double s0 = f(E);
      double s = spline->Eval(E);
      g1->SetPoint(g1->GetN(), E-MTAU, s-s0);
      g2->SetPoint(g2->GetN(), E-MTAU, (s-s0)/s0);
    }
    auto c = new TCanvas;
    c->Divide(1,2);
    c->cd(1);
    g1->Draw("a*l");
    c->cd(2);
    g2->Draw("a*l");
  }
  if(opt.count("fsrc"))
  {
    auto f = [](double v) { return Fr(v);};
    auto spline = interpolate(f,0, 0.99, 1e-6);
    spline->SaveAs("fsrc_minus_one.C");
    TGraph * g0 = new TGraph;
    TGraph * g01 = new TGraph;
    TGraph * g1 = new TGraph;
    TGraph * g2 = new TGraph;
    for(double x=0; x<0.9;x+=0.01)
    {
      double s0 = f(x);
      double s = spline->Eval(x);
      g0->SetPoint(g0->GetN(), x, s0);
      g01->SetPoint(g01->GetN(), x, s);
      g1->SetPoint(g1->GetN(), x, s-s0);
      g2->SetPoint(g2->GetN(), x, (s-s0)/s0);
    }
    auto c = new TCanvas;
    c->Divide(1,3);
    c->cd(1);
    g0->Draw("al");
    g01->SetLineColor(kRed);
    g01->Draw("l");
    c->cd(2);
    g1->Draw("al");
    g1->GetYaxis()->SetLabelSize(0.07);
    g1->GetXaxis()->SetLabelSize(0.07);
    c->cd(3);
    g2->Draw("al");
    g2->GetYaxis()->SetLabelSize(0.07);
    g2->GetXaxis()->SetLabelSize(0.07);
  }
  if(opt.count("fsrc-cmp"))
  {
    FSRC fsrc1(1e-6,0.4,0.0);
    FsrcSpline3 fsrc2(1e-6);
    FsrcSpline3 fsrc3(1e-5);
    using clock = std::chrono::system_clock;
    double step=1e-5;
    auto start1 = clock::now();
    for(double v=0; v<0.4; v+=step) fsrc1(v);
    auto stop1 = clock::now();
    std::cout << "FSRC: " <<  (stop1-start1).count()*1e-9 <<endl;

    auto start2 = clock::now();
    for(double v=0; v<0.4; v+=step) fsrc2(v);
    auto stop2 = clock::now();
    std::cout << "FsrcSpline3: " <<  (stop2-start2).count()*1e-9 <<endl;

    auto start3 = clock::now();
    for(double v=0; v<0.4; v+=step) fsrc3(v);
    auto stop3 = clock::now();
    std::cout << "FsrcSpline3(1e-5): " <<  (stop3-start3).count()*1e-9 <<endl;

    auto start4 = clock::now();
    for(double v=0; v<0.4; v+=step) Fr(v);
    auto stop4 = clock::now();
    std::cout << "Fr(v) " <<  (stop4-start4).count()*1e-9 <<endl;

  }
  if(opt.count("hint"))
  {
    TGraph * g1 = new TGraph;
    for(double v=0.015; v<=0.5;v+=0.01)
    {
      double z = PIALPHA /v;
      //double h0 = (1. - (1.+z)*exp(-z))/(1. - exp(-z))*log(2.0*MTAU*v/ME);
      //double h0 = 0.5*z*(log(2.0*MTAU*v/ME)-5.0/6.0);
      double h = hint2(v);
      g1->SetPoint(g1->GetN(), v, h);
      std::cout << v << "   " << h << endl;
    }
    g1->Draw("al");
  }
  if(opt.count("h"))
  {
    TGraph * g1 = new TGraph;
    for(double v=0.0; v<=0.5;v+=0.0001)
    {
      double hv = h(v,MTAU);
      g1->SetPoint(g1->GetN(), v, hv);
    }
    g1->Draw("al");
  }
  theApp.Run();
  return 0;
}


