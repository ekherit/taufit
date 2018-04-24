/*
 * =====================================================================================
 *
 *       Filename:  fit.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  24.04.2018 16:42:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#pragma once
#include <iostream>
#include <fstream>

#include "draw.h"
#include "fit/TauMassFitter.h"

extern std::string INPUT_FILE;
extern std::string OUTPUT_FILE;
extern std::string RESULT_FILE;
extern std::string EFFCOR_FILENAME;
extern std::string FILTER;
extern double ENERGY_VARIATION;
extern double EFFICIENCY_VARIATION;

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

class Fitter
{
  public:

    bool minos=false;
    bool free_energy = false;
    bool draw_fit_result = true;

    double energy_shift = 0;
    Fitter(void)
    {
    }

    void Fit(const std::list<ScanPoint_t> & SPL)
    {
      TauMassFitter2 fitter;
      if(free_energy) fitter.SetFreeEnergy(true);
      else fitter.SetFreeEnergy(false);
      fitter.Fit(SPL);

      if(minos) { fitter.Minos(); }
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
      if(draw_fit_result) draw_fitresult(fitter, energy_shift);
    }

};
