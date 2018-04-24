/*
 * =====================================================================================
 *
 *       Filename:  sim.h
 *
 *    Description:  simulation helper function
 *
 *        Version:  1.0
 *        Created:  24.04.2018 15:11:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once
#include <list>
#include <vector>

#include "ScanPoint.h"
#include "sigma/Sigma.h"

class ScanSimulator
{
  const std::list<ScanPoint_t> & scenario;
  bool variate_number_of_events=true;
  bool variate_cbs_energy = true;
  ibn::valer<double> tau_mass = {MTAU,0};
  ibn::valer<double> energy_spread = {0,0};
  ibn::valer<double> eps = {0.06,0};
  ibn::valer<double> bg = {0.3,0};
  double cbs_energy_error;
  TRandom R;
  public:
    ScanSimulator(const std::list<ScanPoint_t> & sc, long seed=0) : scenario(sc)
    {
      R.SetSeed(seed);
    }

    std::list<ScanPoint_t> simulate(void)
    {
      std::list<ScanPoint_t> SPL;
      for(auto & sc : scenario)
      {
        ScanPoint_t sp;
        sp = sc;
        double mt = R.Gaus(tau_mass.value, tau_mass.error); //variated tau mass
        sp.energy.value += mt;
        double spread = energy_spread.value == 0 ?  sp.energy_spread.value : energy_spread.value*pow(sp.energy.value/mt,2.0);
        spread = R.Gaus(spread, energy_spread.error);
        double sigma = sigma_total(2.0*sp.energy.value, spread,  mt, 1e-10)*R.Gaus(eps.value, eps.error) + R.Gaus(bg.value, bg.error);
        double mu = sigma*sp.luminosity.value;

        sp.Ntt = variate_number_of_events ? R.Poisson(mu) : mu;

        //variate energy (simulation of CBS measurement) This error for one bunch
        double cbs_error = cbs_energy_error == 0 ? sp.energy_spread.error : cbs_energy_error; 
        //this energy value Wcm/2
        sp.energy.value = variate_cbs_energy ? R.Gaus(sp.energy.value, cbs_error) : sp.energy.value;
        SPL.push_back(sp);
      }
      return SPL;
    }
    void SetTauMass(double value) { tau_mass.value = value; }
    void SetTauMassError(double error) { tau_mass.error = error; }

    void SetEfficiency(double value) { eps.value = value; }
    void SetEfficiencyError(double error) { eps.error = error; }

    void SetBackground(double value) { bg.value = value; }
    void SetBackgroundError(double error) { bg.error = error; }

    void SetEnergySpread(double value) { energy_spread.value = value; }
    void SetEnergySpreadError(double error) { energy_spread.error = error; }
};

