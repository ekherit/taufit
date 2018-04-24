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
#include <algorithm>

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
    //std::cout << minimum << std::endl;
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

extern TauMassFitter2 * TAUMASSFITTER;

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
    bool minos=false;
  };


  std::vector<parinfo_t> inipar; //initial parameter value
  std::vector<parinfo_t> minpar; //optimal parameter value


  public:

  bool isminos=false;
  bool is_free_energy=false;
  double CHI2;
  int NDF;

  TauMassFitter2(void)
  {
    inipar.push_back(  {"M"   , 0    , 0.1 , -1 , +1  , false , true} );
    inipar.push_back(  {"EPS" , 0.06 , 0.02 , 0  , 1   , false , true} );
    inipar.push_back(  {"BG"  , 0.1  , 0.05  , 0  , 5 , false , true} );
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

  double GetChi2(double m, double eps, double bg, double * par=nullptr)
  {
    double chi2 = 0;
    int i=0;
    for(auto & p: SP)
    {
      //expected number of events for estimated parameters
      double E = p.energy.value;
      double S = p.energy_spread.value;
      double L = p.luminosity.value;
      if(par != nullptr) 
      {
        //std::cout  << i << "  " << par[i] << std::endl;
        E+=par[i];
      }
      double nu = L*sigma_obs(E, S,  MTAU+m, eps, bg, p.effcor);
      chi2 += -2*log_likelyhood(nu, p.Ntt);
      i++;
    }
    //if(!std::isnormal(chi2)) chi2=1e100;
    return chi2;
  }

  static void fcn(Int_t& n, Double_t*, Double_t&f, Double_t*par, Int_t)
  {
    if(TAUMASSFITTER->is_free_energy)
    {
      f = TAUMASSFITTER->GetChi2(par[0],par[1],par[2], &par[3]);
      for(int i=0;i<TAUMASSFITTER->SP.size();++i)
      {
        f+=pow( par[3+i]/TAUMASSFITTER->SP[i].energy.error, 2.0);
      }
    }
    else
    {
      f = TAUMASSFITTER->GetChi2(par[0],par[1],par[2]);
    }
  }

  void Fit(void)
  {
    if(is_free_energy) minuit.reset(new TMinuit(inipar.size()+SP.size()));
    else minuit.reset(new TMinuit(inipar.size()));
    int i=0;
    for( auto & p : inipar)
    {
      minuit->DefineParameter(i, p.name.c_str(), p.value,  p.error, p.min , p.max);
      i++;
    }
    if(is_free_energy)
    {
      for(i=inipar.size();i<inipar.size()+SP.size();i++)
      {
        int point = i-inipar.size();
        minuit->DefineParameter(i, ("DE"+std::to_string(point+1)).c_str(), 0, SP[point].energy.error, -5 , 5);
      }
    }
    TAUMASSFITTER = this;
    minuit->SetFCN(TauMassFitter2::fcn);
    isminos = false;
    minuit->Migrad();
    minuit->mnhess();
    minuit->mnmatu(1);
    minuit->mnprin(4,0);
    if(is_free_energy) NDF = 2*SP.size() - minuit->GetNumFreePars();
    else NDF = SP.size() - minuit->GetNumFreePars();

    //std::cout << "minuit->GetNumPars()  = " << minuit->GetNumPars() << std::endl;
    minpar.resize(minuit->GetNumPars());
    for(int i=0;i<minpar.size();i++) 
    {
      if(i<inipar.size()) minpar[i].name = inipar[i].name;
      minuit->GetParameter(i, minpar[i].value, minpar[i].error);
      minpar[i].min = -minpar[i].error;
      minpar[i].max = minpar[i].error;
    }
    //std::cout << "minpar0 = " << minpar[0].value << " 1=" << minpar[1].value << "  2=" << minpar[2].value << endl;
    CHI2 = GetChi2(minpar[0].value, minpar[1].value, minpar[2].value);
  }



  const parinfo_t & operator()(std::string nm) 
  {
    auto & p =  *find_if(minpar.begin(), minpar.end(), [&nm](const parinfo_t & p) { return  p.name == nm; } );
    //std::cout << "in operator: " << p.name << "  " << nm << std::endl;
    return p;
  };



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
  void Minos(int par=-1) 
  {
    isminos = true;
    double eparab, gcc;
    minuit->mnmnos();
    if(par<0 || par > minpar.size()-1) for(int i=0;i<minpar.size();i++) minuit->mnerrs(i, minpar[i].max, minpar[i].min,eparab, gcc);
    else minuit->mnerrs(par, minpar[par].max, minpar[par].min,eparab, gcc);
  }

  void Hesse(void)
  {
  }

  void SetFreeEnergy(bool b=true)
  {
    is_free_energy = b;
  }
};
#endif //TAUMASS_FITTER
