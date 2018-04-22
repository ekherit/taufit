/*
 * =====================================================================================
 *
 *       Filename:  sim.cpp
 *
 *    Description:  Find optimized point's layout for tau scan
 *
 *        Version:  1.0
 *        Created:  15.04.2018 08:27:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#include <vector>
#include <iostream>
#include <chrono>

#include <boost/program_options.hpp>

#include <fmt/printf.h>

using namespace std;

#include <sigma/Sigma.h>

#include <TF1.h>


#include <TROOT.h>
#include <TH1.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TAxis.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TPad.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TFormula.h>
#include <TString.h>
#include <TVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TApplication.h>
#include <TRandom.h>
#include <TFile.h>
R__EXTERN TSystem *gSystem;
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("ROOTPROGRAM","ROOTPROGRAM", initfuncs);

//const double MTAU=1776.86; //PDG2018
//const double DMTAU=0.12;

unsigned DEBUG=0;
double cm_energy(double E1, double E2)
{
  return 2*sqrt(E1*E2)*cos(0.011);
}

//observed cross section
double sigma_obs(double E, double m, double eps, double bg, double Sw)
{
  return sigma_total(2*E,Sw, m,1e-10)*eps + bg;
}

struct point_t
{
  double E;
  double dE;
  double S; //spread
  double dS; 
  double N;
  double L; //lum
};

vector<point_t> Data;
void fcn1(Int_t& n, Double_t*, Double_t&f, Double_t*par, Int_t)
{

  double m = par[0];
  double eps = par[1];
  double bg = par[2];
  double chi2 = 0;
  int i = 0;
  for(auto & point: Data)
  {
    //expected number of events for estimated parameters
    double nu = point.L * sigma_obs(point.E, MTAU+m, eps, bg, point.S);
    chi2 += 2*(nu - point.N + point.N*log(std::max(point.N,1.0)/nu));
    i++;
  }
  if(!std::isnormal(chi2)) chi2=1e100;
  f=chi2;
}

void fcn2(Int_t& n, Double_t*, Double_t&f, Double_t*par, Int_t)
{
  double Ltot = par[4];
  double eps = par[7];
  double bg = par[8];
  double Sw = par[9];
  double DM23 = par[10];

  TMinuit minuit(3);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(fcn1);
  minuit.DefineParameter(0, "dm",     0,  0.037,  -1,  +1);
  minuit.DefineParameter(1, "eff", eps,  0.0018, 0.0, 1.0);      
  minuit.DefineParameter(2, "bg.", bg,  0.14, 0,  10);      

  Data.resize(5);
  //fill energy points
  Data[0].E = par[0]+MTAU;
  Data[1].E = par[1]+MTAU;
  Data[2].E = Data[1].E+DM23;
  Data[3].E = par[2]+MTAU;
  Data[4].E = par[3]+MTAU;

  //fill lums
  Data[0].L = Ltot*par[5];
  Data[1].L = Ltot*par[6];
  Data[2].L = Data[1].L*2/3;
  Data[4].L = (Ltot - Data[0].L - Data[1].L - Data[2].L)/1.5;
  Data[3].L = Data[4].L * 0.5;

  double mass_shift=0, mass_shift_error=0;
  auto print = [&](void)
  {
    cout << "fcn2="<< f << " , mass shift = " << mass_shift*1e3 << " ± "  << mass_shift_error*1e3 <<  " keV" << endl;
    fmt::printf("%5s  %10s  %10s %10s %10s\n", "#","dM, MeV", "Sw,MeV", "L,1/pb","N");
    for(int i=0;i<5;i++)
    {
      fmt::printf("%5d  %10.3f  %10.3f  %10.3f %10.3f\n", i, Data[i].E-MTAU, Data[i].S, Data[i].L, Data[i].N);
    }
  };

  if(Data[4].L < 0) 
  {
    f=1e100*Data[4].L*Data[4].L;
    print();
    return;
  }
  //calculate number of events
  for(int i=0;i<5;i++)
  {
    Data[i].S = Sw*pow(Data[i].E/MTAU, 2.0);
    Data[i].N = Data[i].L*sigma_obs(Data[i].E, MTAU, eps, bg, Data[i].S);
  }

  double pp[3] = {0,eps,bg};
  double * grad = 0;
  double chi2;
  int iflag;
  //minuit.Eval(3,grad, chi2, pp,iflag);
  //std::cout << "chi2 in minimum " << chi2 << endl;
  //minuit.Migrad();
  minuit.mnhess();
  minuit.mnmatu(1);
  //minuit.mnmnos();
  //minuit.mnprin(4,0);
  minuit.GetParameter(0,mass_shift, mass_shift_error);
  static double dm0 = mass_shift_error;
  //cout << m << " " << dm << endl;
  f=(mass_shift_error-dm0)*100;
  print();
  //fmt::printf("%5s  %10s  %10s  %10s\n", "#","dM, MeV", "L,1/pb","N");
  //for(int i=0;i<5;i++)
  //{
  //  fmt::printf("%5d  %10.3f  %10.3f  %10.2f\n", i, Data[i].E-MTAU, Data[i].L, Data[i].N);
  //}
}

void fit(double Ltot=100, double eps=0.06, double bg=0.3, double Sw=0.922, double DM23=0)
{
  TMinuit minuit(11);
  minuit.SetFCN(fcn2);
  minuit.DefineParameter(0 , "dE1" , -5   , 0   , 0  , 0);
  minuit.DefineParameter(1 , "dE2" , -0.3, 0.1, -5 , 3.5);
  minuit.DefineParameter(2 , "dE4" , 3.5  , 0   , 0  , 0);
  minuit.DefineParameter(3 , "dE5" , 15   , 0   , 0  , 0);
  minuit.DefineParameter(4 , "L"   , Ltot , 0   , 0  , 0);
  minuit.DefineParameter(5 , "L1"  , 0.15,  0.08,  0  , 1);
  minuit.DefineParameter(6 , "L2"  , 0.35,  0.06, 0  , 1);
  minuit.DefineParameter(7 , "eps" , eps  , 0   , 0  , 0);
  minuit.DefineParameter(8 , "bg"  , bg   , 0   , 0  , 0);
  minuit.DefineParameter(9 , "Sw"  , Sw   , 0   , 0  , 0);
  minuit.DefineParameter(10 , "DM23", DM23, 0   , 0  , 0);
  minuit.SetMaxIterations(10000);

  minuit.Migrad();
  //minuit.mnsimp();
  minuit.mnhess();
  minuit.mnmatu(1);
  Double_t aminRes,edmRes,errdefRes;
  Int_t nvparRes,nparxRes,icstatRes;
  minuit.mnstat(aminRes,edmRes,errdefRes,nvparRes,nparxRes,icstatRes);
  minuit.mnprin(3,aminRes);
}

int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  double jpsi_spread;
  double tau_spread;
  double psip_spread;
  double energy_spread=1.1*pow(2*MTAU/MJPSI,2.0);
  double bg;
  double eps;
  double Ltotal;
  double DM23;
  opt_desc.add_options()
    ("help,h","Print this help")
    ("jpsi-spread", po::value<double>(&jpsi_spread), "Energy spread on J/psi, MeV")
    ("tau-spread", po::value<double>(&tau_spread), "Energy spread on tau threshold, MeV")
    ("psip-spread", po::value<double>(&psip_spread), "Energy spread on psi(2S), MeV")
    ("bg", po::value<double>(&bg)->default_value(0.3), "Background, pb")
    ("eps", po::value<double>(&eps)->default_value(0.06), "Registration efficiency")
    ("dm23", po::value<double>(&DM23)->default_value(0), "threshold point energy split value, MeV")
    ("luminocity,L", po::value<double>(&Ltotal)->default_value(100), "Total integrated luminocity")
    ;
  po::positional_options_description pos;
  //pos.add("input",-1);
  po::variables_map opt; //options container
  try
  {
    po::store(po::command_line_parser(argc, argv).options(opt_desc).positional(pos).run(), opt);
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


  if(opt.count("jpsi-spread"))
  {
    cout << "Using jpsi-spread " << jpsi_spread << " MeV" << endl;
    energy_spread = jpsi_spread*pow(2*MTAU/MJPSI, 2.0);
  }
  else if (opt.count("psip-spread"))
  {
    cout << "Using psi(2S) energy spread " << psip_spread << " MeV" << endl;
    energy_spread = psip_spread*pow(2*MTAU/MPSIP, 2.0);
  }
  else if (opt.count("tau-spread"))
  {
    cout << "Using tau energy spread " << tau_spread << " MeV" << endl;
    energy_spread = tau_spread;
  }
  std::cout << "Energy spread on tau threshold " <<  energy_spread << " MeV" << endl;
  std::cout << "Background cross section " << bg << " pb" << endl;
  std::cout << "Registration efficiency " << eps << endl;
  std::cout << "Threshold point split energy " << DM23 << endl;
  fit(Ltotal, eps, bg, energy_spread, DM23);
}
