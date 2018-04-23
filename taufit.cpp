/*
 * =====================================================================================
 *
 *       Filename:  taufit.cpp
 *
 *    Description:  Tau mass fit
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


#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TRandom.h>

#include <regex>

#include "fit/TauMassFitter.h"

#include <TROOT.h>
#include <TApplication.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("Fit tau mass","Fit tau mass", initfuncs);

#include "sigma/Sigma.h"

// For measurement calculation speed.
#include <ibn/timer.h>
#include <ibn/averager.h>


//#include "fit/FitFunctions.h"
#include "draw.h"


unsigned DEBUG=0;
std::string INPUT_FILE="scan.txt";
std::string OUTPUT_FILE="fitscan.txt";
std::string EFFCOR_FILENAME="effcor.dat";
std::string FILTER="";
double ENERGY_VARIATION=0.01;
double EFFICIENCY_VARIATION=0.005;
std::string RESULT_FILE="fitresult.txt";
int SEED=0;

void print_result(std::ostream & , TauMassFitter2 & , string comment="");

double correct_energy(std::list<ScanPoint_t> & SPL, double dmjpsi, double dmpsi2s)
{
  double dm1 = dmjpsi/2;
  double dm2 = dmpsi2s/2;
  double m2=3686.109/2;
  double m1=3096.916/2;

  for(auto & sp : SPL)
  {
    double E = sp.energy.value;
    double dE = dm1 + (E-m1) *  (dm2-dm1) / (m2-m1);
    sp.energy.value = E - dE;
    cout << E << " - " << dE << " = " << sp.energy.value << endl;
  }
}


void correct_luminocity2(std::list<ScanPoint_t> & SPL, unsigned lumopt=LUM_GG)
{
  using namespace std;
  SPL.sort([](ScanPoint_t & sp1, ScanPoint_t & sp2) { return sp1.energy.value < sp2.energy.value;});
  double IL=0;
  double INgg=0;
  double INee=0;
  double E0 = SPL.front().energy.value;
  for(auto & sp: SPL)
  {
    IL+=sp.luminosity.value;
    INgg+=sp.Ngg*sq(sp.energy.value/E0);
    INee+=sp.Nee*sq(sp.energy.value/E0);
  }
  double ggS0 = INgg/IL;
  double bbS0 = INee/IL;
  cout << "Estimation gamma-gamma cross section: " << ggS0 << " nb" << endl;
  cout << "Estimation Bhabha cross section: " << bbS0 << " nb" << endl;
  cout << "Correct luminosity using ";
  string lumstr;
  switch(lumopt)
  {
    case LUM_GG:
      lumstr="gamma gamma envents";
      break;
    case LUM_BB:
      lumstr="Bhabha envents";
      break;
  }
  cout << lumstr << endl;
  cout << setw(15) << "energy" << setw(15) << "old lum" << setw(15) << "change" << setw(15) << " new lum " << endl;
  for(auto & sp : SPL)
  {
    double old_lum = sp.luminosity.value;
    switch(lumopt)
    {
      case LUM_GG:
        sp.luminosity.value = sp.Ngg*sq(sp.energy.value/E0)/ggS0;
        sp.luminosity.error = sqrt(sp.Ngg)*sq(sp.energy.value/E0)/ggS0;
        break;
      case LUM_BB:
        sp.luminosity.value = sp.Nee*sq(sp.energy.value/E0)/bbS0;
        sp.luminosity.error = sqrt(sp.Nee)*sq(sp.energy.value/E0)/bbS0;
        break;
    }
    cout << setw(15) << sp.energy.value << setw(15) << old_lum << setw(15) << sp.luminosity.value - old_lum << setw(15) << sp.luminosity.value << endl;
  }
}


using namespace std;
int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  int opt_minos;
  double ESHIFT=0;
  double bg,eps;
  double cbs_energy_error=0.1;
  double jpsi_spread;
  double tau_spread;
  double psip_spread;
  double energy_spread=0;
  opt_desc.add_options()
    ("help,h","Print this help")
    ("input", po::value<std::string>(&INPUT_FILE)->default_value("scan.txt"), "File with data")
    ("output",po::value<std::string>(&OUTPUT_FILE)->default_value("fitscan.txt"), "Output file *.root and *.pdf")
    ("lum",po::value<string>()->default_value("gg"), "luminosity  gg - gamma-gamma, bb-Bhabha, online")
    ("correct-energy","Apply energy correction to CBS energies")
    ("correct-efficiency","Apply efficiency corrections from file effcor.dat")
    ("noeffcor","Not correct efficiency")
    ("effcor",po::value<string>(&EFFCOR_FILENAME), "Efficiency correction file. Default value is effcor.dat ")
    ("spread",po::value<double>()->default_value(1.374),"energy spread" )
    ("mjpsi",po::value<double>()->default_value(0.112), "Difference between measurement and PDG mass of the JPSI" )
    ("mpsi2s",po::value<double>()->default_value(0.238),"Difference between measurement and PDG mass of the PSI2S" )
    ("dmjpsi",po::value<double>()->default_value(0.064), "Error of the jpsi difference")
    ("dmpsi2s",po::value<double>()->default_value(0.113),"Error of the psi2s difference")
    ("ivanos","Draw lfcn")
    ("novp", "No vacuum polarization")
    ("exit", "exit after fitting")
    ("variate-energy",po::value<double>(&ENERGY_VARIATION), "Variate point energy")
    ("variate-efficiency",po::value<double>(&EFFICIENCY_VARIATION), "Variate efficiency correction")
    ("seed",po::value<int>(&SEED), "Random seed")
    ("fitresult",po::value<std::string>(&RESULT_FILE), "Result file which accumulating with result")
    ("filter",po::value<std::string>(&FILTER),"Regex to filter the data")
    ("minos", "Calculate minos errors" )
    ("pdgshift", "Shift energy drawing on PDG mass value")
    ("sim", "Simulation using scenario")
    ("bg", po::value<double>(&bg)->default_value(0.3), "Background, pb")
    ("eps", po::value<double>(&eps)->default_value(0.06), "Registration efficiency")
    ("cbs-energy-error",po::value<double>(&cbs_energy_error), "Energy measurement error")
    ("jpsi-spread", po::value<double>(&jpsi_spread), "Energy spread on J/psi, MeV")
    ("tau-spread", po::value<double>(&tau_spread), "Energy spread on tau threshold, MeV")
    ("psip-spread", po::value<double>(&psip_spread), "Energy spread on psi(2S), MeV")
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


  double Sw = opt["spread"].as<double>();
  if(opt.count("novp")) IS_VP_COR=false;


  std::list<ScanPoint_t> SPL;
  if(opt.count("sim"))
  {
    std::ifstream scenario_file(INPUT_FILE);
    if(!scenario_file) { cerr << "Unable to open scenario file: " << INPUT_FILE << endl; return -1;}
    ScanPoint_t sp;
    TRandom R;
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
    while(scenario_file >> sp.n >>  sp.energy.value >> sp.energy_spread.value >>  sp.luminosity.value >> sp.Ntt)
    {
      sp.energy.value += MTAU;
      sp.energy.error = cbs_energy_error;
      if(energy_spread !=0 )
      {
        sp.energy_spread.value = energy_spread*pow(sp.energy.value/MTAU,2.0);
      }
      double sigma = sigma_total(2*sp.energy.value, sp.energy_spread.value,  MTAU, 1e-10)*eps + bg;
      sp.energy.value = R.Gaus(sp.energy.value, cbs_energy_error*0.5);
      sp.energy_spread.value = Sw*pow(sp.energy.value/MTAU,2.0);
      double mu = sigma*sp.luminosity.value;
      sp.Ntt = R.Poisson(mu);
      //std::cout << "E = " << sp.energy.value << " L = " << sp.luminosity.value <<  "  sigma_obs " << sigma << "  Ntt = " << sp.Ntt << endl;
      SPL.push_back(sp);
    }
    print(SPL);
  }
  else 
  {
    SPL = read_data(INPUT_FILE, Sw);
  }



  if(opt.count("correct-energy"))
  {
    cout << "Correct CBS energies: " << endl;
    correct_energy(SPL, opt["mjpsi"].as<double>(), opt["mpsi2s"].as<double>());
  }

  if(opt.count("effcor"))
  {
    cout << "Correct efficiency from file: " << EFFCOR_FILENAME << endl;
    ifstream file(EFFCOR_FILENAME.c_str());
    if(!file)
    {
      cerr << " No efficiency correction file: " << "effcor.dat" << endl;
      exit(1);
    }
    for(auto & sp :SPL)
    {
      file >> sp.effcor;
      cout <<  sp.energy.value << " " << sp.effcor << endl;
    }
  }

  if(opt.count("noeffcor"))
  {
    cout << "Suppress efficency correction." << endl;
    for(auto & sp :SPL)
    {
      sp.effcor = 1;
    }
  }

  if(opt.count("lum"))
  {
    string s = opt["lum"].as<string>();
		cout << "LUMINOSITY=" << s << endl;
    if(s=="gg") LUMINOCITY = LUM_GG;
    if(s=="bb"|| s=="ee") LUMINOCITY = LUM_BB;
    if(s=="online") LUMINOCITY = LUM_DEFAULT;
		if(LUMINOCITY == LUM_BB || LUMINOCITY == LUM_GG) correct_luminocity2(SPL, LUMINOCITY);
  }

	
	

  if(opt.count("variate-energy"))
  {
    cout << "Variate energy in each point on: " << ENERGY_VARIATION << " MeV" << endl;
    cout << "Random seed is " << SEED << endl;
    TRandom r(SEED);
    cout << setw(5) << "point" << setw(15) << "variation, MeV"  << setw(15) << "new energy, MeV" << endl;
    for(auto & sp : SPL)
    {
      double dE = ENERGY_VARIATION*r.Gaus();
      sp.energy.value = sp.energy.value + dE;
      cout << setw(5) << sp.n << setw(15) << dE << setw(15) << sp.energy.value << endl;
    }
  }

  if(opt.count("variate-efficiency"))
  {
    cout << "Variate efficiency in each point on: " << 100*EFFICIENCY_VARIATION << "%" << endl;
    cout << "Random seed is " << SEED << endl;
    TRandom r(SEED);
    cout << setw(5) << "point" << setw(15) << "variation"  << setw(15) << "new effcor" << endl;
    for(auto & sp : SPL)
    {
      double dEPS = EFFICIENCY_VARIATION*r.Gaus();
      sp.effcor = sp.effcor+sp.effcor*dEPS;
      cout << setw(5) << sp.n << setw(15) << dEPS << setw(15) << sp.effcor << endl;
    }
  }

  TauMassFitter2 fitter;
  fitter.Fit(SPL);
  if(opt.count("minos")) fitter.Minos();
  print_result(std::cout, fitter,"# ");
	
  std::ifstream input_file(INPUT_FILE, std::ios::binary);
  std::ofstream output_file(OUTPUT_FILE, std::ios::binary);
  output_file << input_file.rdbuf();
  print_result(output_file, fitter, "#");

  ofstream result_file(RESULT_FILE.c_str(), fstream::app);
  if(!result_file)
  {
    cerr << "Unable to open file " << RESULT_FILE << endl;
  }
  char buf[65535];
  sprintf(buf,"%8.3f  %+5.3f  %+5.3f  %4.2f  %+4.2f  %+4.2f  %4.2f  %+4.2f  %+4.2f",
      fitter("M").value,         fitter("M").min,       fitter("M").max,     
      fitter("EPS").value,     fitter("EPS").min*100, fitter("EPS").max*100,
      fitter("BG").value,       fitter("BG").min,      fitter("BG").max   
      );
  result_file << buf << endl;
  result_file.close();

	TApplication theApp("tau", &argc, argv);
  if(opt.count("pdgshift"))
  {
    ESHIFT = MTAU;
    std::cout << "Shift the graphic with " << ESHIFT << "  MeV " << endl;
  }
  draw_fitresult(fitter, ESHIFT);

  if(!opt.count("exit")) theApp.Run();
  return 0;
}

void print_result(std::ostream & os, TauMassFitter2 & fitter, std::string comment)
{
  char buf[1024];
  sprintf(buf, "MTAU   = %8.3f +- %5.3f MeV  = %8.3f %+5.3f %+5.3f MeV", fitter("M").value+MTAU , fitter("M").error, fitter("M").value+MTAU, fitter("M").min,fitter("M").max);
  os << comment << buf << endl;

  sprintf(buf, "M-MPDG = %8.3f +- %5.3f MeV", fitter("M").value, fitter("M").error);
  os << comment << buf << endl;
  
  sprintf(buf,"EPS = %3.1f %+3.1f %+3.1f %%", fitter("EPS").value*100, fitter("EPS").min*100,fitter("EPS").max*100);
  os << comment << buf << endl;

  sprintf(buf,"BG = %4.2f %+4.2f %+4.2f pb", fitter("BG").value, fitter("BG").min, fitter("BG").max);
  os << comment << buf << endl;

  sprintf(buf,"chi2/ndf = %6.5f/%d", fitter.CHI2, fitter.NDF);
  os << comment << buf << endl;

  sprintf(buf,"P(chi2,ndf) = %3.1f", TMath::Prob(fitter.CHI2,fitter.NDF));
  os << comment << buf << endl;
}

