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
#include <fstream>
#include <chrono>

#include <boost/program_options.hpp>

#include <fmt/printf.h>

using namespace std;

#include <sigma/Sigma.h>
#include <ScanPoint.h>
#include "fit/TauMassFitter.h"
#include "draw.h"

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

std::string OUTPUT_FILE="fitscan.txt";
TauMassFitter2 * TAUMASSFITTER = nullptr;

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
  double cbs_energy_error=0.01;
  std::string scenario_file_name;
  opt_desc.add_options()
    ("help,h","Print this help")
    ("jpsi-spread", po::value<double>(&jpsi_spread), "Energy spread on J/psi, MeV")
    ("tau-spread", po::value<double>(&tau_spread), "Energy spread on tau threshold, MeV")
    ("psip-spread", po::value<double>(&psip_spread), "Energy spread on psi(2S), MeV")
    ("bg", po::value<double>(&bg)->default_value(0.3), "Background, pb")
    ("eps", po::value<double>(&eps)->default_value(0.06), "Registration efficiency")
    ("dm23", po::value<double>(&DM23)->default_value(0), "threshold point energy split value, MeV")
    ("luminocity,L", po::value<double>(&Ltotal)->default_value(100), "Total integrated luminocity")
    ("input", po::value<std::string>(&scenario_file_name), "Scenario file name")
    ("cbs-energy-error",po::value<double>(&cbs_energy_error), "Energy measurement error")
    ("minos", "Calculate MINOS errors")
    ;
  po::positional_options_description pos;
  pos.add("input",-1);
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

  std::ifstream scenario_file(scenario_file_name);
  if(!scenario_file) { cerr << "Unable to open scenario file: " << scenario_file_name << endl; return -1;}
  double E, L, S, Ntt;
  int n;
  int i = 0;
  std::vector<ScanPoint_t> SP;
  ScanPoint_t sp;
  TRandom R;
  while(scenario_file >>sp.n >>  sp.energy.value >> sp.energy_spread.value >>  sp.luminosity.value >> sp.Ntt)
  {
    sp.energy.value += MTAU;
    double sigma = sigma_obs(sp.energy.value, MTAU,eps, bg, sp.energy_spread.value);
    double mu = sigma*sp.luminosity.value;
    sp.Ntt = R.Poisson(mu);
    std::cout << "E = " << sp.energy.value << " L = " << sp.luminosity.value <<  "  sigma_obs " << sigma << "  Ntt = " << sp.Ntt << endl;
    SP.push_back(sp);
  }
  TauMassFitter2 fitter;
  fitter.Fit(SP);
  if(opt.count("minos")) fitter.Minos();
  draw_fitresult(fitter);
}
