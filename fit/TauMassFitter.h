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
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/FunctionMinimum.h>

#include <list>

#include "../ScanPoint.h"

#include <sigma/Sigma.h>

class TauMassFitter :  public  ROOT::Minuit2::FCNBase 
{
  std::list<ScanPoint_t> * SPL; //list of scan points
  ROOT::Minuit2::MnUserParameters inipar; //initial values of parameters
  ROOT::Minuit2::MnUserParameters minpar; //optimal value of parameters
  const double MTAUSHIFT = MTAU_PDG2011;
  public:
  TauMassFitter(void)
  {
    inipar.Add("M",      0.09, 0.01); //tau mass - M_PDG, MeV
    inipar.Add("EPS",    0.06, 0.2); //efficienc
    inipar.Add("BG",     1,   0.5,  0,  20); //background, pb
    //inipar.Fix(2);
  }

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
      chi2 += 2*(nu - sp.Ntt + sp.Ntt*log(std::max(sp.Ntt,1ul)/nu));
			//cout << chi2 << " M=" << dM+MTAUSHIFT <<  " eps=" << eps << " bg=" << bg << " cross=" << cross_section << "  nu=" << nu << " chi2="<< chi2 << endl;
    }
    if(std::isnan(chi2)) 
		{
			cout << "Bad chi2 " << chi2 << endl;
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
#endif //TAUMASS_FITTER
