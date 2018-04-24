/*
 * =====================================================================================
 *
 *       Filename:  draw.h
 *
 *    Description:  Draw data graph, fit function and fit result with chi square
 *
 *        Version:  1.0
 *        Created:  22.02.2013 16:27:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
#pragma once

#ifndef IBN_TAUFIT_DRAW_H
#define IBN_TAUFIT_DRAW_H
#include <string>
#include <list>

#include "ScanPoint.h"
#include "fit/TauMassFitter.h"

class TGraphErrors; 
class TGraph;

extern std::string OUTPUT_FILE;

extern TGraphErrors * DataGraph(const list<ScanPoint_t> & SPL, string title="r_{i} \\varepsilon \\sigma(e^{+}e^{-}\\rightarrow \\tau^{+}\\tau^{-}) + \\sigma_{B}");
extern TGraph * SigmaGraph(const std::list<ScanPoint_t> & SPL,  double MTAU, double  EFFECT, double BG, int NN );
extern void draw_lum(const std::list<ScanPoint_t> & SPL);
extern void draw_fitresult(TauMassFitter & fitter);
extern void draw_fitresult(TauMassFitter2 & fitter, double MSHIFT=0);
#endif
