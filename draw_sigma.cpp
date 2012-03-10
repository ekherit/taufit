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

#include <iostream>
#include <iomanip>

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


#include "sigma/Sigma.h"

// For measurement calculation speed.
#include <ibn/timer.h>
unsigned DEBUG=0;


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


void draw_dilog(void)
{
  TCanvas * c = new TCanvas("dilog","dilog");
  double NN=1e6;
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
  cout << "Measure time for my dilog..." << endl;
  ibn::timer timer;
  for(double i = 0;i<NN;++i)
  {
    double x = i/NN;
    Li2_int_reg(x);
  }
  double myt=timer.elapsed();
  cout << "Measure time for my dilog..." << endl;
  timer.restart();
  for(double i = 0;i<NN;++i)
  {
    double x = i/NN;
    TMath::DiLog(x);
  }
  double cernt=timer.elapsed();
  cout << "my dilog time = "<<myt << ", cern dilog time = " << cernt << endl;

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
    ("cross-section", "Draw cross section")
    //("vp", po::value<std::string>()->default_value("reh"),"Vacuum polarization")
    ("vp", "Vacuum polarization")
    ("vph","Hadronic vacuum polarization")
    ("dilog", "draw Dilogarithm function")
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
  if(opt.count("dilog")) draw_dilog();
  if(opt.count("vp")) draw_vp();
  if(opt.count("vph")) draw_vp_hadron();
  theApp.Run();
  return 0;
}


