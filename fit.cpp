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
#include <ibn/averager.h>


#include "fit/FitFunctions.h"


unsigned DEBUG=0;

double correct_energy(double dmjpsi, double dmpsi2s)
{
  double dm1 = dmjpsi/2;
  double dm2 = dmpsi2s/2;
  double m2=3686.109/2;
  double m1=3096.916/2;
  
  for(int i=0;i<POINT_NUMBER; i++)
  {
    double dE = dm1 + (ENERGY[i]-m1) *  (dm2-dm1) / (m2-m1);
    double E = ENERGY[i];
    ENERGY[i] = E - dE;
    cout << E << " - " << dE << " = " << ENERGY[i] << endl;
  }
}

double correct_luminocity(unsigned lumopt=LUM_GG)
{
  ibn::averager<double> gg_cross_section_averager;
  ibn::averager<double> bb_cross_section_averager;
  double E0 = ENERGY[0];
  for(int i=0;i<POINT_NUMBER;i++)
  {
    if(LUM[i]==0) continue;
    gg_cross_section_averager.add(NGG[i]*sq(ENERGY[i]/E0)/LUM[i]);
    bb_cross_section_averager.add(NEE[i]*sq(ENERGY[i]/E0)/LUM[i]);
  }
  double ggS0 = gg_cross_section_averager.average();
  double bbS0 = bb_cross_section_averager.average();
  cout << "Estimation gamma-gamma cross section: " << ggS0 << " nb" << endl;
  cout << "Estimation Bhabha cross section: " << bbS0 << " nb" << endl;
  //to change the luminosity
  cout << setw(15) << "energy" << setw(15) << "old lum" << setw(15) << "change" << setw(15) << " new lum " << endl;
  for(int i=0;i<POINT_NUMBER;i++)
  {
    double old_lum = LUM[i];
    switch(lumopt)
    {
      case LUM_GG:
        LUM[i] = NGG[i]*sq(ENERGY[i]/E0)/ggS0;
        break;
      case LUM_BB:
        LUM[i] = NEE[i]*sq(ENERGY[i]/E0)/bbS0;
        break;
    }
    cout << setw(15) << ENERGY[i] << setw(15) << old_lum << setw(15) << LUM[i] - old_lum << setw(15) << LUM[i] << endl;
  }
}

using namespace std;
int main(int argc, char ** argv)
{
  namespace po=boost::program_options;
  po::options_description opt_desc("Allowed options");
  opt_desc.add_options()
    ("help,h","Print this help")
    ("data",po::value<std::string>()->default_value("scan.txt"), "File with data")
    ("spread",po::value<double>()->default_value(1.45514),"energy spread" )
    ("mjpsi",po::value<double>()->default_value(0.116), "Difference between measurement and PDG mass of the JPSI" )
    ("mpsi2s",po::value<double>()->default_value(0.096),"Difference between measurement and PDG mass of the PSI2S" )
    ("dmjpsi",po::value<double>()->default_value(0.072), "Difference between measurement and PDG mass of the JPSI" )
    ("dmpsi2s",po::value<double>()->default_value(0.122),"Difference between measurement and PDG mass of the PSI2S" )
    ("lum",po::value<string>()->default_value("gg"), "luminosity  gg - gamma-gamma, bb-Bhabha")
    ("correct-energy","Apply energy correction to CBS energies")
    ("correct-efficiency","Apply efficiency corrections from file effcor.dat")
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
  if(opt.count("correct-energy"))
  {
    cout << "Correct CBS energies: " << endl;
    correct_energy(opt["mjpsi"].as<double>(), opt["mpsi2s"].as<double>());
  }

  if(opt.count("correct-efficiency"))
  {
    cout << "Correct efficiency from file: " << "effcor.dat" << endl;
    ifstream file("effcor.dat");
    if(!file)
    {
      cerr << " No efficiency correction file: " << "effcor.dat" << endl;
      exit(1);
    }
    for(int i=0; !file.eof() && i<POINT_NUMBER; i++)
    {
      file >> EFCOR[i];
      cout <<  i << " " << EFCOR[i] << endl;
    }
  }

  if(opt.count("lum"))
  {
    string s = opt["lum"].as<string>();
    if(s=="gg") LUMINOCITY = LUM_GG;
    if(s=="bb") LUMINOCITY = LUM_BB;
    correct_luminocity(LUMINOCITY);
  }
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
//	minuit->mnparm( 2,   "BG",    1,     0.5,   0,  100000, ierflg);
	minuit->mnparm( 2,   "BG",    1,     0.5,   0,  0, ierflg);
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

  //int npar=3;
  //double fcn=0;
  //int iflag=0;
  //Lfcn(npar,0,fcn, par,iflag);
  //cout << "2*L/ndf = " << fcn*2<<"/("<<POINT_NUMBER<<"-"<<npar<<")="<< 2*fcn/(POINT_NUMBER-npar) << endl;
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

