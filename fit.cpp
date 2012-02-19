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


unsigned DEBUG=0;

using namespace std;
int main(int argc, char ** argv)
{
  string filename="test.dat";
  if(argc>1) filename=argv[argc-1];
	ifstream file(filename.c_str());
	if(!file) { cerr << "cant open file " << filename << endl; exit(0);}
	cout << "Reading data from file " << filename << endl;
	FillData2(file,1.6);
	TGraphErrors * data_gr = DataGraph("r_{i} \\varepsilon \\sigma(e^{+}e^{-}\\rightarrow \\tau^{+}\\tau^{-}) + \\sigma_{B}");
	
	TMinuit *minuit = new TMinuit(3); 
	minuit->SetFCN(Lfcn);

	Double_t arglist[10];
	Int_t ierflg = 0;

 	minuit->SetErrorDef(0.5);
	minuit->mnparm( 0, "MTAU", 0,   0.2,   -1,    +1, ierflg);
	minuit->mnparm( 1, "EFF",  0.1,  0.1,  0,   1, ierflg);
	minuit->mnparm( 2, "BG",    0,  0.4,   0,   20 , ierflg);
	arglist[0]=2;
	minuit->mnexcm("SET STR", arglist ,1,ierflg);
	minuit->Migrad();
	//minuit->mnhess();
	//minuit->mnhess();
	//minuit->mnmnos();
  cout << "Fit is finished" << endl;
  //minuit->mnexcm("MINOS", arglist ,1,ierflg);
	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	minuit->mnprin(4,amin);
  
	double par[3],erpar[3];
	for( int i =0 ; i < 3; i++)
  {
		minuit->GetParameter(i,par[i], erpar[i]);
	}
	TApplication theApp("tau", &argc, argv);
	TCanvas * sigma_c = new TCanvas("sigma", "sigma", 800,900); 
	sigma_c->SetGrid();
	data_gr->Draw("ap");
	TGraph * sigma_gr = SigmaGraph(par[0],par[1],par[2],100);
	sigma_gr->Draw("c");

  theApp.Run();
  return 0;
}


