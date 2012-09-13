
#include <iomanip>

#include <TMinuit.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <algorithm>

#include <sigma/Sigma.h>


extern unsigned DEBUG;
const double MTAUSHIFT = MTAU_PDG2011;
//double PRECISION = 1e-3;
double DDelta=1;

Double_t sigma_total (Double_t *x, Double_t *p)	{
	return sigma_total(*x, p[0], p[1],0.001);
}

//double sigma_all (double *x, double *p)	{
//	return sigma_all(4 *x[0] * x[0], p[0], 0.001);
//}

double sigma_fit (double *x, double *p)	{
	//p[0] - m_tau
	//p[1] -- epsilon (efficiency)
	return p[1]*sigma_total(2*x[0], DDelta, p[0]+1800,0.001);
}


//double fit_function(double mtau, double eps, double bg)
//{
//  return eps * sigma + fabs(par[2]));
//}

double factorial(unsigned N)	{
	if(N<=0) return 1;
	double f=1.;
	for(unsigned i = 1; i <=N ; i++)	f*=i;
	return f;
}

double logfactorial(unsigned N)	{
	if( N <=0) return 0;
	double f=0;
	for( unsigned i = 1; i<=N; i++) f+=log(i);
	return f;
}
int POINT_NUMBER;
const int MAX_POINT_NUMBER=1024;

double DELTA[MAX_POINT_NUMBER]; // Energy distribution
int    EVENT[MAX_POINT_NUMBER]; //amount of registred event
double ENERGY[MAX_POINT_NUMBER]; //energy of beam
double ENERGY_ERROR[MAX_POINT_NUMBER]; //energy error
double EFCOR[MAX_POINT_NUMBER]; //corrections to efisiency
double LUM[MAX_POINT_NUMBER];

long  NGG[MAX_POINT_NUMBER]; //number of  gamma-gamma events
long  NEE[MAX_POINT_NUMBER]; //number of Bhabha events


enum 
{
  LUM_GG,
  LUM_BB
};

unsigned LUMINOCITY=LUM_GG;

bool ischar(istream &is,char Ch) 
{
    char c;
    is.get(c);
    is.putback(c);
    if( c == Ch)
    {
        return true;
    }
    return false;
} 

void FillData(istream & file)
{
	double tmp;
	int i;
	int colw=20;
	cout << setw(6) << "POINT"<<setw(colw)<<"E[MeV]"<<setw(colw) << "DELTA[MeV]" << setw(colw) << "LUM[1/pb]"  << setw(6) << "N"   << setw(5) << "R" << endl; 
	file.ignore(1024,'\n');
	for(i=0; !file.eof(); i++)
  {
		file >> tmp >> ENERGY[i] >> DELTA[i] >> LUM[i] >>EVENT[i] >> EFCOR[i];
		//EFCOR[i]=1;
		LUM[i]/=1000; //????????? ?? ????????? ?????????? ? ???????? ?????????
		if(file.eof()) break;
		cout << setw(6) << i+1 << setw(colw) << ENERGY[i] << setw(colw)  << DELTA[i] << setw(colw) << LUM[i] << setw(6) << EVENT[i] << setw(5) << EFCOR[i] << endl;
	}
	POINT_NUMBER = i;
	cout << "POINT_NUMBER=" << POINT_NUMBER << endl;
	double lum=0;
	for(int i = 0; i < POINT_NUMBER; i++)
  {
		 lum+= LUM[i];
	}
	cout << "TOTAL LUMINOCITY=" <<lum << endl;
}

void FillData2(istream & file, double sigmaW_psi2s /* energy spread at psi resonance */)
{
	double tmp;
  double Sw; 
  double dSw;
  double W;
  double dW;
  double lum;
  double lum_cor;
  unsigned Ntt; // Number of tau events
  unsigned Nee;
  unsigned Ngg;
	int colw=20;
  double efcor;
	cout << setw(6) << "POINT"<<setw(colw)<<"E[MeV]"<<setw(colw) << "DELTA[MeV]" << setw(colw) << "LUM[1/pb]"  << setw(6) << "N"   << setw(5) << "R" << endl; 
	int i;
	for(i=0; !file.eof(); i++)
  {
    char c = file.get();
    if(c=='#')
    {
      file.ignore(65535,'\n');
      i--;
      continue;
    }
		if(file.eof()) break;
    file.putback(c);
    file >> tmp >> lum >> W >> dW >> Sw >> dSw  >>  Ntt >> Nee >> Ngg >> efcor;
    EVENT[i] = Ntt;
    ENERGY[i] = W/2;
    ENERGY_ERROR[i]=dW/2;
    DELTA[i]  = sigmaW_psi2s*pow(W/3686.109,2);
		LUM[i] = lum/=1000; //????????? ?? ????????? ?????????? ? ???????? ?????????
    EFCOR[i] = efcor;
    NEE[i]=Nee;
    NGG[i]=Ngg;
    file.ignore(65535,'\n');
		if(file.eof()) break;
		cout << setw(6) << i+1 << setw(colw) << ENERGY[i] << setw(colw)  << DELTA[i] << setw(colw) << LUM[i] << setw(8) << EVENT[i] << setw(12) << EFCOR[i] << endl;
	}
	POINT_NUMBER = i;
	cout << "POINT_NUMBER=" << POINT_NUMBER << endl;
  lum=0;
	for(int i = 0; i < POINT_NUMBER; i++)
  {
		 lum+= LUM[i];
	}
	cout << "TOTAL LUMINOCITY=" <<lum << endl;
}

double InPol(int N, double *X, double *Y, double x )	{
	int imin=0;
	int imax=N;
	for( int i = 0; i < N; i++)	{
		if ( x >= X[i] && X[i] > X[imin]) imin=i; 
		if ( x <= X[i]  && X[i] < X[imax]) imax = i;
	}
	if(imin == imax && x < X[imin]) 	{
		if(imax <N-1) imax++;
		else imax --;
	}
	if(imin == imax && x > X[imin]) 	{
		if(imin <N-1) imin++;
		else imin --;
	}
  if(imin==imax) return Y[imin];
	return Y[imin] + (x - X[imin])/(X[imax] - X[imin] )*(Y[imax] - Y[imin]);
}

double InPol2(unsigned N, double *X, double *Y, double x )
{
  double x1=-std::numeric_limits<double>::max();
  double x2=+std::numeric_limits<double>::max();
  double y1,y2,y;
  bool f1=false, f2=false; //found something flag
  for(unsigned i=0;i<N;++i)
  {
    //find maximal value x1 that lower then x
    if(x>=X[i] && X[i]>x1)
    {
      x1=X[i];
      y1=Y[i];
      f1=true;
    }
    //find minimal value x2 that higher then x
    if(x<=X[i] && X[i]<x2)
    {
      x2=X[i];
      y2=Y[i];
      f2=true;
    }
  }
  //if found something good then calculate the result
  if(f1  && f2) y = x1==x2 ? y1 :y1+(x-x1)*(y2-y1)/(x2-x1);
  if(!f1 && f2) y=y2;
  if(!f2 && f1) y=y1;
  if(!f1 && !f2) 
  {
    cerr << "ERROR: Unable to found interpolation for: " << x << endl;
  }
  //cout << "x1="<< x1 << " y1="<<y1 << ", x2="<<x2 << ", y2="<<y2 << ", x=" << x << ", y=" << y << endl;
  return y;
}


double Delta( double E)
{
	/*
	int imin=0;
	int imax=POINT_NUMBER;
	for( int i = 0; i < POINT_NUMBER; i++)	{
		if ( E >= ENERGY[i] && ENERGY[i] > ENERGY[imin]) imin=i; 
		if ( E < ENERGY[i]  && ENERGY[i] < ENERGY[imax]) imax = i;
	}
	if(imin == imax && E < ENERGY[imin]) 	{
		if(imax <POINT_NUMBER-1) imax++;
		else imax --;
	}
	if(imin == imax && E > ENERGY[imin]) 	{
		if(imin <POINT_NUMBER-1) imin++;
		else imin --;
	}
	return DELTA[imin] + (E - ENERGY[imin])/(ENERGY[imax] - ENERGY[imin] )*(DELTA[imax] - DELTA[imin]);
	*/
	return InPol(POINT_NUMBER, ENERGY, DELTA, E);
}

inline double logfac(double n)
{
  return n*log(n) - n + (log(n*(1.+4.*n*(1.+2.*n))))/6. + LOGPI2;
}

static double FCN_PREVIOS_VALUE=1000;
void Lfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{  
 	 static unsigned FCN_CALL_NUMBER=0;	
	 double LL=0;
	 double mu;
	 double sigma;
	 double fac;
	 for ( int i = 0; i<POINT_NUMBER; i++)
   {
     sigma = sigma_total(2*ENERGY[i], DELTA[i], par[0]+MTAUSHIFT,PRECISION);
     mu = LUM[i] *  (fabs(par[1])*EFCOR[i] * fabs(sigma) + fabs(par[2]));
     if(mu <= 0 )
     {
       cerr << "Error: mu < 0\n";
       cerr << i  << " " << ENERGY[i] << " mu=" << mu << " sigma=" << sigma << endl;
       cerr << "mtau="<<par[0]<<"  eps=" << par[1] << "  bg = " << par[2] << endl;
       f =  1e300*( 1 - mu);
       return;
     }
     
     //if(sigma == 0)
     //{	
     //  cerr.precision(12);
     //  cerr << i << " " << ENERGY[i] <<  ", mu = " << mu ;
     //  cerr << ", Sigma = " << sigma  << endl;
     //  cerr << "trying to increez accuracy" << endl;
     //  sigma = sigma_total(2*ENERGY[i], DELTA[i], par[0]+MTAUSHIFT,PRECISION/1000.);
     //  mu = LUM[i] * ( par[1] * EFCOR[i] * sigma + par[2]);
     //  if( sigma <= 0)
     //  {
     //    cerr << " Sigma is steel below zeor";
     //  }
     //  exit(1);
     //}

     if ( EVENT[i] < 0)
     {
       cerr << "ERROR: number of event is low the zero" << endl;
       exit(1);
     }
     //if(EVENT[i] <= 50) fac = 	log( factorial(EVENT[i]) );
     //else fac = EVENT[i]*log(EVENT[i])-EVENT[i];
     if(EVENT[i] <= 5) fac = 	log( factorial(EVENT[i]) );
     else fac = logfac(EVENT[i]);
     LL+= EVENT[i] * log(mu) - mu -  fac ;
     cout.precision(6);
     //cout << ENERGY[i]<<": "<<sigma<<", "<<DELTA[i] <<", pars: "<<par[0]<<", "<<par[1]<<", "<<par[2]<< ": mu="<<mu << endl;
   }
	 f= -LL;
	 FCN_PREVIOS_VALUE = f;
}

/*
TGraph * SigmaGraph(double mtau, double  EFFECT, double BG, int NN = 40)	{
	if(NN > MAX_POINT_NUMBER) NN=MAX_POINT_NUMBER;
	double S[MAX_POINT_NUMBER];
	double EE[MAX_POINT_NUMBER];
	double EEmin=ENERGY[0];
	double EEmax=ENERGY[0];
	double efcor=1;
	double delta=1;
	for(int i = 0; i< POINT_NUMBER; i++) if(EEmin > ENERGY[i]) EEmin = ENERGY[i];
	for(int i = 0; i< POINT_NUMBER; i++) if(EEmax < ENERGY[i]) EEmax = ENERGY[i];
	for(int i = 0; i < NN; i++)	{
		EE[i] = EEmin + (EEmax - EEmin)/(NN-1)*i;
		efcor = InPol(POINT_NUMBER, ENERGY, EFCOR, EE[i]);
		delta = InPol(POINT_NUMBER, ENERGY, DELTA, EE[i]);
		S[i] =  EFFECT*sigma_total(2*EE[i],delta,mtau+MTAUSHIFT,0.001)*efcor +BG;
		//cout <<  i << " "<< EE[i] << "    " << S[i] << "  " << Delta(EE[i]) << endl;
		//cout << " mtau = " << mtau << " MTAUSHIFT = " << MTAUSHIFT << endl;
	}
	return new TGraph(NN,EE,S);
}
*/
TGraph * SigmaGraph(double MTAU, double  EFFECT, double BG, int NN )
{
  BG = fabs(BG);
  EFFECT = fabs(EFFECT);
	if(NN > MAX_POINT_NUMBER) NN=MAX_POINT_NUMBER;
	double S[MAX_POINT_NUMBER];
	double EE[MAX_POINT_NUMBER];
	double EEmin=ENERGY[0];
	double EEmax=ENERGY[0];
	double lum=1, efcor = 1,delta = 1;
	for(int i = 0; i< POINT_NUMBER; i++) if(EEmin > ENERGY[i]) EEmin = ENERGY[i];
	for(int i = 0; i< POINT_NUMBER; i++) if(EEmax < ENERGY[i]) EEmax = ENERGY[i];
  EEmin-=5;
  EEmax+=5;
  //calculate parameters for DELTA
  TF1 f("fun_spread","[0]*x*x");
  f.SetParameter(0,0);
  TGraph g(POINT_NUMBER, ENERGY,DELTA);
  g.Fit("fun_spread","goff");
  double Spar=f.GetParameter(0);
	for(int i = 0; i < NN; i++)
  {
		EE[i] = (EEmin + (EEmax - EEmin)/(NN-1)*i);
    delta = Spar*sq(EE[i]);
		S[i]=(EFFECT*sigma_total(2*EE[i],delta,MTAU+MTAUSHIFT,0.0001) + BG);
	}
	return new TGraph(NN,EE,S);
}

TGraph * SigmaGraph(TMinuit * myMin, int NN = 40)	{
	double MTAU, MTAU_error;
	double EFFECT, EFFECT_error;
	double BG, BG_error;
	myMin->GetParameter(0,MTAU,MTAU_error);
	myMin->GetParameter(1,EFFECT,EFFECT_error);
	myMin->GetParameter(2,BG,BG_error);
	return SigmaGraph(MTAU,EFFECT,BG,NN);
}


TGraph * SigmaGraphNew(double MTAU, double  EFFECT, double BG, int NN )	{
	if(NN > MAX_POINT_NUMBER) NN=MAX_POINT_NUMBER;
	double S[MAX_POINT_NUMBER];
	double EE[MAX_POINT_NUMBER];
	double EEmin=ENERGY[0];
	double EEmax=ENERGY[0];
	double lum=1, efcor = 1,delta = 1;
	for(int i = 0; i< POINT_NUMBER; i++) if(EEmin > ENERGY[i]) EEmin = ENERGY[i];
	for(int i = 0; i< POINT_NUMBER; i++) if(EEmax < ENERGY[i]) EEmax = ENERGY[i];
	for(int i = 0; i < NN; i++)	{
		EE[i] = (EEmin + (EEmax - EEmin)/(NN-1)*i);
		double delta = 1.07;
		S[i]=sigma_total(2*EE[i],delta,MTAU+MTAUSHIFT,0.0001);
	}
	return new TGraph(NN,EE,S);
}

//TGraph * SigmaGraphN(TMinuit * myMin, int NN = 40)	{
//	double MTAU, MTAU_error;
//	double EFFECT, EFFECT_error;
//	double BG, BG_error;
//	myMin->GetParameter(0,MTAU,MTAU_error);
//	myMin->GetParameter(1,EFFECT,EFFECT_error);
//	myMin->GetParameter(2,BG,BG_error);
//	return SigmaGraphN(MTAU,EFFECT,BG,NN);
//}


//TGraphErrors * DataGraph(const char * title )	{
//	double S[MAX_POINT_NUMBER];
//	double Ser[MAX_POINT_NUMBER];
//	double Wer[MAX_POINT_NUMBER];
//	for(int i = 0; i < POINT_NUMBER; i++ )	{
//		Wer[i] = 0;
//		S[i] = EVENT[i]/LUM[i]/EFCOR[i];
//		Ser[i] = sqrt(EVENT[i])/LUM[i]/EFCOR[i];
//	}
//	TGraphErrors * gr = new TGraphErrors(POINT_NUMBER, ENERGY,S,Wer,Ser);
//	gr->SetTitle(title);
//	gr->Draw("A*");	
//	gr->GetXaxis()->SetTitle("W/2 [MeV]");
//	gr->GetYaxis()->SetTitle("\\frac{N}{L \\varepsilon r } [pb]");
//	return gr;
//}

TGraphErrors * DataGraph(const char * title)
{
	double S[MAX_POINT_NUMBER];
	double Ser[MAX_POINT_NUMBER];
	for(int i = 0; i < POINT_NUMBER; i++ )
  {
		S[i] = EVENT[i]/LUM[i]/EFCOR[i];
		Ser[i] = sqrt(EVENT[i])/LUM[i]/EFCOR[i];
		//S[i] = EVENT[i]/LUM[i];
		//Ser[i] = sqrt(EVENT[i])/LUM[i];
    cout << i << " " << ENERGY[i] << " " << LUM[i] << " " << EVENT[i] << " " << S[i] << endl;
	}
	TGraphErrors * gr = new TGraphErrors(POINT_NUMBER, ENERGY,S,ENERGY_ERROR,Ser);
	gr->SetTitle(title);
	//gr->Draw("A*");	
	//gr->GetXaxis()->SetTitle("E, MeV");
	//gr->GetYaxis()->SetTitle("\\frac{N}{L} [pb]");
	return gr;
}


TGraphErrors * DataGraphN(const char * title )	{
	double Ner[MAX_POINT_NUMBER];
	double N[MAX_POINT_NUMBER];
	double Eer[MAX_POINT_NUMBER];
	double E[MAX_POINT_NUMBER];
	for(int i = 0; i < POINT_NUMBER; i++ )	{
		N[i]=EVENT[i];
		Ner[i] = sqrt(EVENT[i]);
		Eer[i] = 0;
		E[i] = 2*ENERGY[i];
	}
	TGraphErrors * gr = new TGraphErrors(POINT_NUMBER, E,N,Eer,Ner);
	gr->SetTitle(title);
	gr->Draw("A*");	
	gr->GetXaxis()->SetTitle("W/2, MeV");
	gr->GetYaxis()->SetTitle("N");
	return gr;
}


double F2(double *x, double *p)	{
	int npar;
	double gin;
	double f;
	int iflag;
	double par[3] = {x[0],p[0],x[1]};
	Lfcn(npar, &gin, f, par, iflag);
	return f;
}

double fcn_value(double * p)	{
	double fedm;
	int npari, nparx, istat;
  double fcn_minimum;
	Lfcn(npari,&fedm, fcn_minimum, p, istat);
	return fcn_minimum;
}

void get_parameters(TMinuit * minuit,double *par, double *par_err )	{
	for( int i =0 ; i < 3; i++)	{
		minuit->GetParameter(i,par[i], par_err[i]);
	}
}

double fcn_migrad_value(TMinuit * minuit)	{
		double par[3];
		double par_er[3];
		minuit->Migrad();
		//minuit->mnhess();
		//minuit->mnimpr();
		get_parameters(minuit,par,par_er);
		return fcn_value(par) ;
}

TGraph * ivanos(TMinuit * minuit, double min, double max, double step=0.01)	{
	double MTAU, MTAU_error,mt;
	double EFFECT, EFFECT_error;
	double BG, BG_error;
	minuit->GetParameter(0,MTAU,MTAU_error);
	minuit->GetParameter(1,EFFECT,EFFECT_error);
	minuit->GetParameter(2,BG,BG_error);
	if(max < 0) max = 0;
	if(min > 0) min = 0;
	if(step < 0) step = -step;
	if(step == 0) step = 0.01;
	int NL = -int(min/step);
	int NR = int (max/step);
	double * M = new double[NL+NR];
	double * FCN = new double[NL+NR];
	int ief;
	double fmin, fedm, errdef;
	int npari, nparx, istat;
	double par[3];
	double p_tmp;
  double fcn_minimum;
	minuit->GetParameter(0,par[0],MTAU_error);
	minuit->GetParameter(1,par[1],EFFECT_error);
	minuit->GetParameter(2,par[2],BG_error);
	Lfcn(npari,&fedm, fcn_minimum, par, istat);
	double left_error=0,  right_error=0;
	minuit->FixParameter(0);
	minuit->SetPrintLevel(0);
	ofstream res("ivanos.dat");
	for(int i = 0; i <  NL; i++)	{
		M[i] = MTAU - i*step;
		minuit->mnparm( 0, "MTAU", M[i],   0.2,   0,    0, ief);
		minuit->Migrad();
		minuit->GetParameter(0,par[0],MTAU_error);
		minuit->GetParameter(1,par[1],EFFECT_error);
		minuit->GetParameter(2,par[2],BG_error);
		Lfcn(npari,&fedm, fmin, par, istat);
		FCN[i] = fmin;
		if(FCN[i] >= fcn_minimum + 0.5) left_error = -i*step;
		res << M[i] << "  " << FCN[i] << endl;
	}
	for(int i = 0; i <  NR; i++)	{
		M[i+NL] = MTAU + i*step;
		minuit->mnparm( 0, "MTAU", M[i+NL],   0.2,   0,    0, ief);
		minuit->Migrad();
		minuit->GetParameter(0,par[0],MTAU_error);
		minuit->GetParameter(1,par[1],EFFECT_error);
		minuit->GetParameter(2,par[2],BG_error);
		Lfcn(npari,&fedm, fmin, par, istat);
		FCN[i+NL] = fmin;
		if(FCN[i+NL] >= fcn_minimum + 0.5) right_error = i*step;
		res << M[i+NL] << " " << FCN[i+NL] << endl;
	}
	if(left_error == 0)	{
		cout << endl;
		cout << "*********WARNING: Negative error for MTAU not calculated. You should expand negative region" << endl;
	}
	if(right_error == 0)	{
		cout << endl;
		cout << "*********WARNING: Positive error for MTAU not calculated. You should expand positive region" << endl;
	}
	cout << "*************************************************************** " << endl;
	cout << "     IVANOS ERRORS FOR MTAU\n";
	cout << "*************************************************************** " << endl;
	cout << " MTAU = " << MTAU << " + " << right_error << " - " << left_error << endl;
	//res  << "# MTAU = " << MTAU << " + " << right_error << " - " << left_error << endl;
	
	TGraph * gr = new TGraph(NR+NL, M, FCN);
	return gr;
}


struct Errors { double value, negative, positive;};
struct Point {
	double E; //energy
	double L; //part of luminocity in 
};

Errors ivanos(TMinuit * minuit,  int np=0, double step=1)	{
	const int MAX_ITERATIONS=100;
	int count=0;;
	double par[3], par_er[3];
	double p[3], p_er[3];
	double fcn_minimum;
	double neg, pos;
	int ierf;
	double prec=0.001;
	minuit->Migrad();
	minuit->mnhess();
	minuit->mnimpr();
	get_parameters(minuit,par,par_er);
	get_parameters(minuit,p,p_er);
	fcn_minimum = fcn_value(par);
	double fcnup = fcn_minimum + 0.5;
	minuit->FixParameter(np);
	double fcn_h=fcn_minimum, fcn_l=fcn_minimum, fcn_tmp, fcn_prev;
	double p_h, p_l, p_prev;
	minuit->SetPrintLevel(0);
	while( fcn_h < fcnup )	{
		fcn_l = fcn_h;
		p[np]=p[np]-step;
		minuit->mnparm(np, "MTAU", p[np],   step,   0,    0, ierf);
		fcn_h = fcn_migrad_value(minuit);
		cout << "find level for negative error " << p[np] << " fcn_h = " << fcn_h << " fcn_l = " << fcn_l << endl;
	}	
	p_h = p[np];
	p_l = p[np]+step;
	cout << "minimum = " << fcn_minimum << endl;
	cout << " p_h = " << p_h << "  p_l = " << p_l << endl;
	while (  fabs(fcn_h - fcn_l) > prec && count++ < MAX_ITERATIONS)	{
		p_prev=p[np];
		fcn_prev = fcn_tmp;
		p[np] = p_l + (p_h - p_l)/2.;
		minuit->mnparm(np, "MTAU", p[np],   step,   0,    0, ierf);
	 	fcn_tmp = fcn_migrad_value(minuit);	
		if(fcn_tmp  > fcnup) {
			fcn_h = fcn_tmp;
			p_h = p[np];
		}
		if(fcn_tmp < fcnup) {
			fcn_l = fcn_tmp;
			p_l = p[np];
		}
		if(fcn_tmp == fcnup)	break;
		cout << "******************************************************************\n\n";
		cout << setw(15) <<"****" <<setw(15) <<"high" << setw(15) << "low" << setw(15) << "prev" << endl;
		cout << setw(15)<<"param:"<<setw(15) << p_h  << setw(15) <<  p_l  << setw(15) <<  p_prev   << endl;
		cout << setw(15)<<"fcn:" <<setw(15) << fcn_h << setw(15) <<  fcn_l << setw(15) <<  fcn_prev << endl;
		cout << setw(15) << "Level: " << setw(15) << fcnup << endl;
		cout << "******************************************************************\n\n";
		if(fcn_h < fcnup || fcn_l > fcnup ) {
			cerr << " ERROR fcn_h or fcn_l out of range\n";
			exit(1);
		}
	}
	neg = par[np] - p_h;
	if(count >= MAX_ITERATIONS) {
		neg=0;
		cerr << "WARNIGN: cant determine negative errors\n";
	}


	
	fcn_h = fcn_minimum;
	fcn_l = fcn_minimum;
	p[0]=par[0];
	p[1]=par[1];
	p[2]=par[2];
	count = 0;
	while( fcn_h < fcnup)	{
		fcn_l = fcn_h;
		minuit->mnparm(np, "MTAU", p[np]+=step,   step,   0,    0, ierf);
		fcn_h = fcn_migrad_value(minuit);
		cout << p[np] << "finding levele for positive error:  fcn_h = " << fcn_h << " fcn_l = " << fcn_l << endl;
	}	
	p_h = p[np];
	p_l = p[np]-step;
	while (  fabs(fcn_h - fcn_l) > prec  && count++ < MAX_ITERATIONS)	{
		p_prev=p[np];
		fcn_prev = fcn_tmp;
		p[np] = p_l + (p_h - p_l)/2.;
    minuit->mnparm(np, "MTAU", p[np],   step,   0,    0, ierf);
	 	fcn_tmp = fcn_migrad_value(minuit);	
		if(fcn_tmp > fcnup) {
			fcn_h = fcn_tmp;
			p_h = p[np];
		}
		if(fcn_tmp < fcnup) {
			fcn_l = fcn_tmp;
			p_l = p[np];
		}
		if(fcn_tmp == fcnup)	break;
		cout << "******************************************************************\n\n";
		cout << setw(15) << "****" <<setw(15) <<"high" << setw(15) << "high" << setw(15) << "prev" << endl;
		cout << setw(15)<<"p:"<<setw(15) << p_h  << setw(15) <<  p_l  << setw(15) <<  p_prev   << endl;
		cout << setw(15)<<"fcn:" <<setw(15) << fcn_h << setw(15) <<  fcn_l << setw(15) <<  fcn_prev << endl;
		cout << setw(15) << "Level: " << setw(15) <<fcnup << endl;
		cout << "******************************************************************\n\n";
		//cout << "Positive: " << p[np] << " fcn_h = " << fcn_h << " fcn_l = " << fcn_l << " fcn_tmp = " << fcn_tmp<< endl;
		if(fcn_h < fcnup || fcn_l > fcnup ) {
			cerr << " ERROR fcn_h or fcn_l out of range\n";
			exit(1);
		}
	}
	pos = p_h-par[np];
	if(count >= MAX_ITERATIONS) {
		pos=0;
		cerr << "WARNIGN: cant determine positive errors\n";
	}
	cout << "*************************************************************** " << endl;
	cout << "     IVANOS ERRORS FOR PARAMETER " << np << endl;
	cout << "*************************************************************** " << endl;
	cout << " P" << np << " = " << par[np] << " + " << pos << " - " << neg << endl;
	minuit->Release(np);
	Errors err = {par[np], -neg, pos};
	return err;
}
 
struct Spread_t
{
  double E; //energy
  double D; //spread
};

void simulate_data ( TRandom *rnd, Point * SCAN, double lum, double mtau, double epsilon, double bg)	{
	double mu;
	double sigma;
	Spread_t spPsi = {3096.917/2., 0.9 };
	Spread_t spPsiPrime = { 3686.111/2.,1.33 };

	for(int i =0; i< POINT_NUMBER; i++)	{
		ENERGY[i] = SCAN[i].E + MTAU;
		DELTA[i] =spPsi.D + (spPsiPrime.D - spPsi.D)/(spPsiPrime.E - spPsi.E)*(ENERGY[i] - spPsi.E); 
		LUM[i] = SCAN[i].L*lum;
		sigma =  sigma_total(2*ENERGY[i],DELTA[i],mtau,PRECISION);
		mu = LUM[i]*( epsilon*sigma + bg ); 	
		EVENT[i] = rnd->Poisson(mu);
		EFCOR[i]=1;
		cout << i << "   " << ENERGY[i] << "    " << DELTA[i] << "    " << LUM[i] << "    " << EVENT[i] << "    " << EFCOR[i] <<  endl; //" mu = " << mu << ", sigma = " << sigma << "bp" << endl;
	}	
}

void simulate_data ( TRandom *rnd, double mtau, double epsilon, double bg)	{
    double mu;
	double sigma;
	cout << "POINT_NUMBER = " << POINT_NUMBER << endl;
	for(int i =0; i< POINT_NUMBER; i++)	{
		sigma =  sigma_total(2*ENERGY[i],DELTA[i],mtau,PRECISION);
		mu = LUM[i]*( EFCOR[i]*epsilon*sigma + bg ); 	
		EVENT[i] = rnd->Poisson(mu);
		cout << i << "   " << ENERGY[i] << "    " << DELTA[i] << "    " << LUM[i] << "    " << EVENT[i] << "    " << EFCOR[i] <<  endl; //" mu = " << mu << ", sigma = " << sigma << "bp" << endl;
	}	
}


TGraph * contur(TMinuit *minuit,int p1, int p2, int N,double level=0.5)	
{   double * x = new double[N];
    double * y = new double[N];
    int rc;
    minuit->SetErrorDef(level);
    minuit->mncont(p1,p2,N,x,y,rc);
    return new TGraph(N,x,y);
}

void contur(TMinuit * minuit, int p1, int p2, int N,double * levels, int lnum){   
    TCanvas * contur_c = new TCanvas("contur", "contur", 800,900); 
    TH2F * h = new TH2F("h","\\varepsilon - m\\tau",100,-1,+2,100,0.01,0.07);
    h->Draw("a");
    for (unsigned i = 0; i<lnum; i++)	{
	TGraph * gr = contur(minuit,p1,p2,N, levels[i]);
	if (levels[i] == 0.5) gr->SetLineWidth(3);
	gr->Draw("sl");
	contur_c->Modified();
	contur_c->Update();
    }
}
void contur(TMinuit * minuit){   
    double levels[] = {5,4,3,2,1,0.5,0.4,0.3,0.2,0.1,0.001};
    contur(minuit,0,1,100,levels,11);
}

TGraph * ivanos3(TMinuit * minuit, double mmin, double mmax, double dm)
{
	minuit->SetFCN(Lfcn);
	Double_t arglist[10];
	Int_t ierflg = 0;
 	minuit->SetErrorDef(0.5);
	minuit->mnparm( 0, "MTAU",    0,   0.1,  0, 0, ierflg);
	minuit->mnparm( 1,  "EFF", 1e-3,   0.005,   0,  1, ierflg);
	minuit->mnparm( 2,   "BG",    1,     0.5,   0,  100000, ierflg);
  TGraph * graph = new TGraph;
  TCanvas * c = new TCanvas;
  graph->Draw("a*");
  int i=0;
  for(double M = mmin; M<mmax; M+=dm)
  {
    minuit->mnparm( 0, "MTAU", M,   0.1,  0, 0, ierflg);
    minuit->FixParameter(0);
    minuit->Migrad();
    double par[3],er[3];
    double p_tmp;
    minuit->GetParameter(0,par[0],er[0]);
    minuit->GetParameter(1,par[1],er[1]);
    minuit->GetParameter(2,par[2],er[2]);
    double fcn, fedm, errdef;
    int npar, nparx, istat;
    Lfcn(npar, &fedm, fcn, par, istat);
    cout << "M=" << M << " "  << fcn << endl;
    graph->SetPoint(i, M, fcn);
    c->Modified();
    c->Update();
    i++;
  }
  return graph; 
};
