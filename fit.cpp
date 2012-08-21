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


#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/program_options.hpp>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TLatex.h>


#include <TROOT.h>
#include <TApplication.h>
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("Fit tau mass","Fit tau mass", initfuncs);

#include "sigma/Sigma.h"

// For measurement calculation speed.
#include <ibn/timer.h>


#include "fit/FitFunctions.h"


unsigned DEBUG=0;

using namespace std;
int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  opt_desc.add_options()
    ("help,h","Print this help")
    ("data",po::value<std::string>()->default_value("scan.txt"), "File with data")
    ("spread",po::value<double>()->default_value(1.43),"energy spread" )
    ("ivanos","Draw lfcn")
    ("novp", "No vacuum polarization")
//    ("fix","Fix parameter")
    ;
  po::positional_options_description pos;
  pos.add("data",-1);
  po::variables_map opt; //options container
  try
  {
    po::store(po::command_line_parser(argc, argv).options(opt_desc).positional(pos).run(), opt);
    //std::ifstream config_file("fit.cfg");
    //po::store(po::parse_config_file(config_file,opt_desc,true), opt);
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
  string filename=opt["data"].as<string>();
  if(argc>1) filename=argv[argc-1];
	ifstream file(filename.c_str());
	if(!file) { cerr << "cant open file " << filename << endl; exit(0);}
	cout << "Reading data from file " << filename << endl;

	FillData2(file,Sw);
	//FillData(file);
	TGraphErrors * data_gr = DataGraph("r_{i} \\varepsilon \\sigma(e^{+}e^{-}\\rightarrow \\tau^{+}\\tau^{-}) + \\sigma_{B}");
	
	TMinuit *minuit = new TMinuit(3); 
	minuit->SetFCN(Lfcn);
	TApplication theApp("tau", &argc, argv);
  if(opt.count("ivanos"))
  {
    TCanvas * ivc = new TCanvas;
    TGraph * g = ivanos3(minuit, -1,+1,0.1);
    g->Draw("al");
    theApp.Run();
  }

	Double_t arglist[10];
	Int_t ierflg = 0;

 	minuit->SetErrorDef(0.5);
	minuit->mnparm( 0, "MTAU",    0,   0.5,  0, 0, ierflg);
	minuit->mnparm( 1,  "EFF",  0.01,   0.2,   0,  1, ierflg);
	minuit->mnparm( 2,   "BG",    1,     0.5,   0,  100000, ierflg);
	arglist[0]=2;
	minuit->mnexcm("SET STR", arglist ,1,ierflg);
	minuit->Migrad();
	minuit->Migrad();
	minuit->mnhess();
	minuit->mnmnos();
  cout << "Fit is finished" << endl;
  //arglist[0]=1e10;
  //minuit->mnexcm("MNMNOS", arglist ,1,ierflg);
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
  double M = par[0] + MTAUSHIFT;
  double dM = erpar[0];
  double EPS = par[1]*100; //efficiency %
  double dEPS = erpar[1]*100; 
  double BG = abs(par[2]); //background
  double dBG = abs(erpar[2]); //background
  char buf[1024];
  sprintf(buf, "MTAU   = %7.2f +- %4.2f MeV", M , erpar[0]);
  cout << buf << endl;
  sprintf(buf, "M-MPDG = %7.2f +- %4.2f MeV", par[0], erpar[0]);
  cout << buf << endl;
  //cout << "M - MPDG = " <<  par[0] << " +- " <<  erpar[0] << " MeV" << endl;
	TCanvas * sigma_c = new TCanvas("sigma", "sigma", 1.2*640,1.2*480); 
	sigma_c->SetGrid();
	data_gr->Draw("ap");
	TGraph * fit_gr = SigmaGraph(par[0],par[1],par[2],200);
	fit_gr->Draw("c");
  data_gr->GetXaxis()->SetTitle("E, MeV");
  data_gr->GetYaxis()->SetTitle("#sigma_{obs}, pb");

  double x,y;
  x=M-5;
  y=15;
  std::vector <double> v(&data_gr->GetY()[0], &data_gr->GetY()[data_gr->GetN()]);
  std::sort(v.begin(),v.end());
  y=v.back();

  char texbuf[1024];
  //sprintf(texbuf,"M_{#tau} = %7.2f #pm %4.2f^{%+4.2f}_{%+4.2f} MeV", M, dM, minuit->fErp[0],minuit->fErn[0]);
  sprintf(texbuf,"M_{#tau} = %7.2f^{%+4.2f}_{%+4.2f} MeV", M, minuit->fErp[0],minuit->fErn[0]);
  TLatex Mtex(x,y,texbuf);
  Mtex.Draw();
  sprintf(texbuf,"#varepsilon = %3.1f #pm %3.1f %%", EPS, dEPS);
  TLatex EPStex(x,y*0.9,texbuf);
  EPStex.Draw();
  //sprintf(texbuf,"#sigma_{BG} = %3.1f #pm %3.1f pb", BG, dBG);
  sprintf(texbuf,"#sigma_{BG} = %3.1f^{%+4.2f}_{%+4.2f} pb", BG, minuit->fErp[2],minuit->fErn[2]);
  TLatex BGtex(x,y*0.8,texbuf);
  BGtex.Draw();

  sprintf(texbuf,"M_{#tau}-M_{PDG} =%4.2f #pm %4.2f MeV ",M-MTAUSHIFT , sqrt(dM*dM + DMTAU_PDG*DMTAU_PDG));
  TLatex DMtex(x,y*0.6,texbuf);
  DMtex.Draw();

  TGraphErrors * glum = new TGraphErrors(POINT_NUMBER);
  for(int i=0;i<POINT_NUMBER; i++)
  {
    double r = double(NEE[i])/double(NGG[i]);
    glum->SetPoint(i,ENERGY[i], r);
    double dr = r*sqrt(1./NGG[i] + 1./NEE[i]);
    glum->SetPointError(i, ENERGY_ERROR[i], dr);
  }
  TCanvas * clum = new TCanvas;
  glum->SetMarkerStyle(21);
  glum->Draw("ap");

  theApp.Run();

  return 0;
}


