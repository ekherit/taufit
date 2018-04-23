/*
 * =====================================================================================
 *
 *       Filename:  Fitter.h
 *
 *    Description:  Fit 
 *
 *        Version:  1.0
 *        Created:  15.02.2013 17:52:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
#ifndef IBN_TAUMASS_FITTER_H
#define IBN_TAUMASS_FITTER_H
#include <vector>
#include <list>
#include <chrono>

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/FunctionMinimum.h>


#include "../ScanPoint.h"

#include <sigma/Sigma.h>
#include <sigma/SigmaInterpolation.h>

#include <TMinuit.h>

class TauMassFitter :  public  ROOT::Minuit2::FCNBase 
{
  std::list<ScanPoint_t> * SPL; //list of scan points
  ROOT::Minuit2::MnUserParameters inipar; //initial values of parameters
  ROOT::Minuit2::MnUserParameters minpar; //optimal value of parameters
  const double MTAUSHIFT = MTAU_PDG2018;
  public:
  TauMassFitter(void)
  {
    inipar.Add("M",      0.09, 0.01); //tau mass - M_PDG, MeV
    inipar.Add("EPS",    0.06, 0.2); //efficienc
    inipar.Add("BG",     1,   0.5,  0,  20); //background, pb
    //inipar.Fix(2);
  }

  mutable SigmaInterpolated fSigmaInterpolated;

  //double sigma_total(double W, double spread, double M, double prec) const
  //{
  //  return fSigmaInterpolated(W,spread,M,prec);
  //}


  double operator() (const std::vector<double> & par) const
  {
    double  dM = par[0]; //mass shift respect to PDG value
    double  eps = par[1]; //efficiency
    double  bg = par[2]; //background
    double chi2=0;
    for(const ScanPoint_t & sp : *SPL)
    {
      double cross_section = sigma_total(2*sp.energy.value, sp.energy_spread.value, dM+MTAUSHIFT,PRECISION);
      double visible_cross_section = eps*sp.effcor*cross_section + bg;
      double nu = sp.luminosity.value * visible_cross_section; //expected number of events
      //cout << sp.energy.value << " " << sp.energy_spread.value << " " << MTAUSHIFT << " " << endl;
      chi2 += 2*(nu - sp.Ntt + sp.Ntt*log(std::max(sp.Ntt,1.0)/nu));
			//cout << chi2 << " M=" << dM+MTAUSHIFT <<  " eps=" << eps << " bg=" << bg << " cross=" << cross_section << "  nu=" << nu << " chi2="<< chi2 << endl;
    }
    if(std::isnan(chi2)) 
		{
			//cout << "Bad chi2 " << chi2 << endl;
			return 1e100;
		}
    return chi2;
  } 
  double Up() const { return 1.; }
  
  
  void Fit(std::list<ScanPoint_t> & spl)
  {
    SPL=&spl;
    using namespace ROOT::Minuit2;
    MnMigrad migrad(*this, inipar);
    FunctionMinimum minimum = migrad();
    std::cout << minimum << std::endl;
    minpar=minimum.UserParameters();
    DM =ibn::valer<double>(minpar.Value(0), minpar.Error(0));
    EPS=ibn::valer<double>(minpar.Value(1), minpar.Error(1));
    BG =ibn::valer<double>(minpar.Value(2), minpar.Error(2));
    M = DM + MTAUSHIFT;
    CHI2 = (*this)(minpar.Params());
    NDF = SPL->size() - migrad.VariableParameters();
    MnMinos minos(*this, minimum);
    errDM = minos(0,1000000);
    errEPS = minos(1,1000000);
    errBG = minos(2,1000000);
  }
  ibn::valer<double> DM;
  ibn::valer<double> M;
  ibn::valer<double> EPS;
  ibn::valer<double> BG;
  //negative and positive errors for parameters
  std::pair<double,double> errDM;
  std::pair<double,double> errEPS;
  std::pair<double,double> errBG;
  double CHI2;
  int NDF;

  const list<ScanPoint_t> & GetData(void) { return *SPL; }
};

class TauMassFitter2;

TauMassFitter2 * TAUMASSFITTER;

class TauMassFitter2
{
  std::unique_ptr<TMinuit> minuit;
  std::vector<ScanPoint_t> SP;
  struct parinfo_t
  {
    std::string name;
    double value; //initial value
    double error; //error or first step;
    double min;
    double max;
    bool fixed;
    bool limited;
  };


  std::vector<parinfo_t> inipar; //initial parameter value
  std::vector<parinfo_t> minpar; //optimal parameter value


  public:

  bool isminos=false;
  ibn::valer<double> DM;
  ibn::valer<double> M;
  ibn::valer<double> EPS;
  ibn::valer<double> BG;
  //negative and positive errors for parameters
  std::pair<double,double> errDM;
  std::pair<double,double> errEPS;
  std::pair<double,double> errBG;
  double CHI2;
  int NDF;

  TauMassFitter2(void)
  {
    inipar.push_back(  {"M"   , 0    , 0.1 , -1 , +1  , false , true} );
    inipar.push_back(  {"EPS" , 0.06 , 0.1 , 0  , 1   , false , true} );
    inipar.push_back(  {"BG"  , 0  , 0.1  , 0  , 5 , false , true} );
  }


  //observed cross section
  double sigma_obs(double E, double Sw, double m, double eps, double bg, double effcor=1.0)
  {
    double sigma = sigma_total(2*E, Sw, m, 1e-10);
    return sigma*eps*effcor + bg;
  }

  double log_likelyhood(double nu, double N)
  {
    //nu - expected number of events from fit model
    //N - measured number of events
    return -(nu - N + N*log(std::max(N,1.0)/nu));
  }

  double GetChi2(double m, double eps, double bg)
  {
    double chi2 = 0;
    for(auto & p: SP)
    {
      //expected number of events for estimated parameters
      double E = p.energy.value;
      double S = p.energy_spread.value;
      double L = p.luminosity.value;
      double nu = L*sigma_obs(E, S,  MTAU+m, eps, bg, p.effcor);
      chi2 += -2*log_likelyhood(nu, p.Ntt);
    }
    //if(!std::isnormal(chi2)) chi2=1e100;
    return chi2;
  }

  static void fcn(Int_t& n, Double_t*, Double_t&f, Double_t*par, Int_t)
  {
    f = TAUMASSFITTER->GetChi2(par[0],par[1],par[2]);
  }

  void Fit(void)
  {
    minuit.reset(new TMinuit(inipar.size()));
    int i=0;
    for( auto & p : inipar)
    {
      minuit->DefineParameter(i, p.name.c_str(), p.value,  p.error, p.min , p.max);
      i++;
    }
    TAUMASSFITTER = this;
    minuit->SetFCN(TauMassFitter2::fcn);
    isminos = false;
    minuit->Migrad();
    minuit->mnhess();
    minuit->mnmatu(1);
    minuit->mnprin(4,0);
    NDF = minuit->GetNumFreePars();

    minuit->GetParameter(0, DM.value,  DM.error);
    minuit->GetParameter(1, EPS.value, EPS.error);
    minuit->GetParameter(2, BG.value,  BG.error);
    M.value = DM.value + MTAU;
    M.error = DM.error;
    CHI2 = GetChi2(DM.value, EPS.value, BG.value);
  }

  void Fit(const std::list<ScanPoint_t> & sp)
  {
    SP.reserve(sp.size());
    for(auto & p : sp) SP.push_back(p);
    Fit();
  }
  void Fit(const std::vector<ScanPoint_t> & sp)
  {
    SP = sp;
    Fit();
  }
  
  list<ScanPoint_t>  GetData(void)
  {
    list<ScanPoint_t> spl;
    for(auto & p: SP ) spl.push_back(p);
    return spl;
  }
  //Calculate minos errors
  void Minos(void) 
  {
    isminos = true;
    double eparab, gcc;
    minuit->mnmnos();
    minuit->mnerrs(0, errDM.second, errDM.first,eparab, gcc);
    //std::cout << "M : " << errDM.first << " " << errDM.second << "  " << eparab << "  " << gcc << endl;
    minuit->mnerrs(1, errEPS.second, errEPS.first,eparab, gcc);
    //std::cout << "EPS : " << errEPS.first << " " << errEPS.second << "  " << eparab << "  " << gcc << endl;
    minuit->mnerrs(2, errBG.second, errBG.first,eparab, gcc);
    //std::cout << "BG : " << errBG.first << " " << errBG.second << "  " << eparab << "  " << gcc << endl;
  }
};
#endif //TAUMASS_FITTER
