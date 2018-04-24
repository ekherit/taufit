/*
 * =====================================================================================
 *
 *       Filename:  ScanPoint.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  22.02.2013 15:54:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#ifndef IBN_TAUFIT_SCAN_POINT_H
#define IBN_TAUFIT_SCAN_POINT_H

#include <ibn/valer.h>
#include "sigma/Const.h"

#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <boost/format.hpp>
#include <fmt/printf.h>
#include <fmt/format.h>

using namespace std;
enum 
{
	LUM_DEFAULT, 
  LUM_GG,
  LUM_BB
};

extern unsigned LUMINOCITY;

struct ScanPoint_t
{
  int n; //point id
  ibn::valer<double> energy; //beam energy
  ibn::valer<double> energy_spread; //center of mass energy spread
  ibn::valer<double> luminosity; //luminosity
  double Ntt=0.0; //number of selected tau pair
  double Ngg=1.0; //number of gamma-gamma events
  double Nee=1.0; //number of Bhabha events
  double effcor=1.0; //correction to registration efficiency
};

inline std::list<ScanPoint_t> read_data(std::string fname, double sigmaW_mtauPDG) 
{
  std::list<ScanPoint_t> SPL;
  ifstream file(fname);
  if(!file) 
  {
    std::cerr << "Unable to open file : " << std::endl;
    return SPL;
  }
  std::string FILTER = "";
  std::regex filter(FILTER);
  std::regex comment(" *#.*");
  std::regex empty(" *");

  int colw=15;
  std::cout << boost::format("%6s%24s%15s%15s%15s%15s")%"POINT"%"E[MeV]"%"DELTA[MeV]"%"LUM[1/pb]"%"EVENT"%"EFCOR" << std::endl; 
  int i=0;
  typedef boost::format fmt;
  double total_luminosity=0;
  while(!file.eof())
  {
    std::string line;
    getline(file,line);
    std::smatch match;
    if(std::regex_match(line,match,filter)) continue;
    if(std::regex_match(line,match,comment)) continue;
    if(std::regex_match(line,match,empty)) continue;
    ScanPoint_t sp;
    std::istringstream is(line);
    double tmp;
    is >> tmp;
    is >> sp.luminosity.value; 
		sp.luminosity.error=0;
		is >> sp.luminosity.error;
    is >> sp.energy.value;
    is >> sp.energy.error;
    is >> sp.energy_spread.value;
    is >> sp.energy_spread.error;
    is >> sp.Ntt;
    is >> sp.Nee;
    is >> sp.Ngg;
    is >> sp.effcor;
    sp.energy.value=sp.energy.value/2; //convert to beam energy
    sp.energy.error=sp.energy.error/2; //convert to beam energy
    sp.energy_spread.value = sigmaW_mtauPDG*pow(sp.energy.value/MTAU_PDG2011,2);
    sp.energy_spread.error = 0;
    sp.luminosity.value=sp.luminosity.value/1000; //convert lum into picobarn
    sp.luminosity.error=sp.luminosity.error/1000; //convert lum into picobarn
    sp.n = i;
    SPL.push_back(sp);
    std::cout << boost::format("%6d%15.3f +- %5.3f%15.3f%15.3f%15d%15.4f") 
      % (i+1) 
      % sp.energy.value % sp.energy.error 
      % sp.energy_spread.value 
      % sp.luminosity.value 
      % sp.Ntt
      % sp.effcor << std::endl;
    i++;
    total_luminosity += sp.luminosity.value;
  }
  std::cout << "TOTAL LUMINOCITY=" << total_luminosity <<  " pb" << std::endl;
  //sort by energy order
  SPL.sort([](ScanPoint_t & p1, ScanPoint_t &p2) { return p1.energy.value < p2.energy.value; });
  return SPL;
}

inline void print(const std::list<ScanPoint_t> & SPL, std::ostream & os=std::cout) 
{
  //os << boost::format("%6s%24s%15s%15s%15s%15s")%"point"%"e[mev]"%"delta[mev]"%"lum[1/pb]"%"event"%"efcor" << std::endl; 
  os << fmt::format("{:5}{:>24}{:>15}{:>15}{:>15}{:>15}", 
      "#PNT", "E[MeV]", "DELTA[MeV]", "∫L[1/pb]", "EVENT", "EFFCOR") << std::endl; 
//σw
  int i = 0;
  for (auto &sp : SPL)
  {
    os << boost::format("%5d%15.3f +- %5.3f%15.3f%15.3f%15d%15.4f") 
      % (i+1) 
      % sp.energy.value % sp.energy.error 
      % sp.energy_spread.value 
      % sp.luminosity.value 
      % sp.Ntt
      % sp.effcor << std::endl;
    i++;
  };
}

#endif
