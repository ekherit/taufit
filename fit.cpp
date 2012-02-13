/*
 * =====================================================================================
 *
 *       Filename:  draw_sigma.cpp
 *
 *    Description:  Draw tau prodaction cross section.
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

#include <TROOT.h>
#include <TApplication.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("draw_tau_cross_section","Draw tau cross section", initfuncs);

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TAxis.h>

#include <fstream>
#include <iostream>

#include "sigma/Sigma.h"

// For measurement calculation speed.
#include <ibn/timer.h>


#include "fit/FitFunctions.h"

using namespace std;
int main(int argc, char ** argv)
{
	TApplication theApp("tau", &argc, argv);
	ifstream file("BESTAU-6.dat");
	if(!file) { cerr << "cant open file BESTAU.dat" << endl; exit(0);}
	cout << "Reading BES data...\n";
	FillData(file);

	/*
	TCanvas * sigma_c = new TCanvas("sigma_c", "Sigma", 800,900);	
	sigma_c->SetGrid();
	sigma_c->SetFillColor(0);
	sigma_c->Divide(1,2);
	sigma_c->cd(1);
	TGraphErrors * data_gr = DataGraphN("BES data");
	sigma_c->cd(2);
	TGraphErrors * data_gr2 = DataGraph("BES data");
	sigma_c->Update();
	*/
	
	TMinuit *myMin = new TMinuit(3); 
	myMin->SetFCN(Lfcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

 	myMin->SetErrorDef(0.5);
	myMin->mnparm( 0, "MTAU", 0,   0.2,   0,    0, ierflg);
	myMin->mnparm( 1, "EFF",  0.0426,  0.02,  0,  0.1, ierflg);
	myMin->mnparm( 2, "BG",    0,  0.4,   0,   10 , ierflg);
	//myMin->FixParameter(2);
	//myMin->FixParameter(1);
	// Now ready for minimization step
	arglist[0]=2;
	myMin->mnexcm("SET STR", arglist ,1,ierflg);
	myMin->Migrad();
	myMin->mnhess();
	myMin->mnhess();
	myMin->mnmnos();
  //	arglist[0]=10000;
  //myMin->mnexcm("MINOS", arglist ,1,ierflg);
	// Print results
	//Double_t amin,edm,errdef;
	//Int_t nvpar,nparx,icstat;
	//myMin->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//myMin->mnprin(4,amin);
	//gr=ivanos(myMin, -0.7, +0.9, 0.005);
	//ivanos(myMin,0,0.7);
	//exit(0);
	double par[3],erpar[3];
	get_parameters(myMin,par, erpar);
	TCanvas * sigma_c = new TCanvas("sigma", "sigma", 800,900); 
	sigma_c->SetGrid();
	TGraph * sigma_gr = SigmaGraph(par[0],par[1],par[2],40);
	TGraphErrors * data_gr = DataGraph("r_{i} \\varepsilon \\sigma(e^{+}e^{-}\\rightarrow \\tau^{+}\\tau^{-}) + \\sigma_{B}");
	data_gr->Draw("ap");
	sigma_gr->Draw("sl");
  theApp.Run();
  return 0;
}


