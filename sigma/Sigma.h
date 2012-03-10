#ifndef IBN_SIGMA0_H
#define IBN_SIGMA0_H

#include <iostream>

#include <ibn/integral.h>
#include "VacuumPolarization.h"
#include "FinRadCor.h"
#include "IniRadCor.h"


extern unsigned DEBUG;

using namespace std;

inline   double sigma_total(double W, double delta, double mt, double prec);

inline double sigma_final_radcor( double s, double mt);
inline double sigma_all_radcor( double s, double mt, double prec);
inline double sigma_tree( double s, double mt);


static bool IS_VP_COR=true;
static bool NO_INI_RADCOR=false;


inline double velocity(double s, double mt )
{ 
    return sqrt(1.-sq(mt*2)/s); 
}

inline double velocity( double gamma2)
{ 
    return sqrt(1.- 1./gamma2); 
}

inline bool is_below_threshold(double gamma2)
{
	return  (gamma2 < 1);
}

inline bool is_below_threshold(double s, double mt)
{
	return  (s < 4*mt*mt);
}
/* Name convention:
 *  tree --- tree level
 *  vp  ---- vacuum polarization
 *  vpl ---- lepton vacuum polarization 
 *  vph --- hadron vacuum polarization
 *  fsrc --- final state radiative corections
 *  fc  --- colomb interection in final state
 *  total --- cross section with all radiative corrections and with energy spread
 *  isrc  --- initias state radiative corrections
 */

// ***************  Crossection at tree level  ***************************************** //

inline double sigma_tree(double s,double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*v*(3.-v*v)/2./s;
}

inline double sigma_tree( double g)
{
	if( is_below_threshold(g) ) return 0;
	double v=velocity(g);
	return v*(3.-v*v)/g;
}


inline double sigma_tree_v( double s,double mt)
{
	double v=velocity(s,mt);
	return SIGMA_CONST*(3.-v*v)/2./s;
}


inline double sigma_tree_v( double gamma2)
{
	double v=velocity(gamma2);
	return (3.-v*v)/gamma2;
}

inline double sigma_tree_fc( double s, double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*(3.-v*v)/2./s*vFc(v);
}

inline double sigma_tree_fc_vpl(double s, double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*(3.-v*v)/2./s*vFc(v)*Vp_lepton(s);
}

inline double sigma_tree_fc_vpl_vph( double s, double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*(3.-v*v)/2./s*vFc(v)*Vp_lepton(s)*Vp_hadron(s);
}


inline double sigma_tree_fc_vpl_fsrc( double s, double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*(3.-v*v)/2./s*vFc(v)*Vp_lepton(s)*fsrc(v);
}



inline bool test_bad(double x)
{
	return std::isnan(x) || std::isinf(x);
};



inline double sigma_tree_fc_vp_fsrc( double s, double mt)	
{
  double sigma=0; //cross section
  double v=0; //velocity
  double vfc=1; //colomb corrections
  double vp=1; //vacuum polarization
  double frc=1; //final state radiative corrections
  if(DEBUG==1) cout << "s="<<s <<", mt=" << mt << std::endl;
	if(std::isnan(s) || std::isinf(s)) std::cout << "Bad s=" << s << std::endl;
  // Рабочая функция для интегрирования. Учтены все поправки в конечном сосстоянии и поляризация вакуума.
	if ( is_below_threshold(s,mt) ) return 0;
	v=velocity(s,mt);
	vfc=vFc(v);
	if(IS_VP_COR) vp=VP(s,mt);
	frc=fsrc(v);
	if(std::isnan(v) || std::isinf(v)) std::cout << "vFc("<<v<<") = " << vfc << endl;
	if(test_bad(vp)) std::cout << "vp("<<sqrt(s/4) <<") = " << vp << endl;
	if(test_bad(frc)) std::cout << "fsrc("<<v<<") = " << frc << endl;
	sigma= SIGMA_CONST*(3.-v*v)/(2*s)*vfc*vp*frc;
  if(std::isnan(sigma) || std::isinf(sigma)) sigma=0;
  return sigma;
}

// Bukin's sequence 

inline double sigma_tree_vpl(double s, double mt)
{
	if ( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*v*(3.-v*v)/2./s*Vp_lepton(s);;
}

inline double sigma_tree_vp(double s, double mt)
{
	if ( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*v*(3.-v*v)/2./s*VP(s,mt);;
}

inline double sigma_tree_vp_fsrc(double s, double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*v*(3.-v*v)/2./s*VP(s,mt)*fsrc(v);
}

inline double sigma_tree_vp_fsrc_fc(double s, double mt)
{
	if( is_below_threshold(s,mt) ) return 0;
	double v=velocity(s,mt);
	return SIGMA_CONST*(3.-v*v)/2./s*VP(s,mt)*fsrc(v)*vFc(v);
}




// Crossection only with final radiative corrections. 
inline double sigma_final_radcor(double gamma2)
{
	if ( is_below_threshold(gamma2) ) return 0;
	double v=velocity(gamma2);
	double fr = Fr(v);
	if(fr<0)	{
		cerr << "warning: Fr < 0!!!" << endl;
		exit(1);
	}
	//return sigma_tree_v(gamma2)*vFc(v)*fr;
	return sigma_tree_v(gamma2)*fsrc(v);
}

// Cross sectin with final state colomb interaction and final state radiative corrections //
inline double sigma_final_radcor(double s,double mt)
{
	double g2 = s/sq(2.*mt);
	return sigma_final_radcor(g2)*SIGMA_CONST/2./sq(2*mt);
}


// Cross section with final state radiative corrections and vacuum polarisation of intermediate photon
inline double sigma_final_vp_radcor ( double s, double mt )
{
	return sigma_tree_fc_vpl_fsrc(s,mt);
	double  gamma2 = s/sq(2.*mt);
	if(gamma2<=1) return 0;
	double sig = sigma_final_radcor(gamma2)*SIGMA_CONST/2./sq(2*mt);
	//if(IS_VP_COR) sig*=Vp(s);//sig*=VP_tau(s, mt); //ATTENTION here is approximate value off vacuum polarization
	if(IS_VP_COR) sig*=VP_tau2(s, mt); //ATTENTION here is approximate value off vacuum polarization
	return sig;
}


/*

class SigmaFinVP 	{
	const double mtau;
	public:
		SigmaFinVP(double mt) : mtau(mt) {};
		double operator()(double s)	{ return sigma_final_vp_radcor(s, mtau); }
};
*/
/*
double sigma0(double *x, double *p)	{
	if(is_below_threshold(sq(*x),p[0])) return 0;
	return sigma_tree(sq(*x),p[0]);
}

double sigma1(double *x, double *p)	{
	double mt = p[0];
	double s = *x;
	return sigma_final_radcor(s, mt);
}
*/
/*Подинтегральное выражение для использования в алгоритме интегрирования */


class F_sigma
{
	const double s;
	const double mt;
	public:
	F_sigma(double s_, double mt_) : s(s_), mt(mt_){}
	double operator()(double x)
  {
		return F(x,s)*sigma_final_vp_radcor(s*(1.-x),mt);
	}
};




// Сечение с учетом всех радпоправок.
inline double sigma_all_radcor_sing(double s, double mt, double prec)
{
	double eps = 1. - sq(2*mt)/s;	
	if (eps < 0) return 0;
	double xmax = min(eps,2*60/sqrt(s));
	double I = 0,integ;
	integ=ibn::dgaus(F_sigma(s,mt),0,xmax,prec);
	I+=integ;
	return I;
}

class F_reg1_sigma
{
	const double mt;
	const double s;
	double L;
	double b;
public:
	F_reg1_sigma(double s_, double mt_)	: mt(mt_), s(s_)
  {
		L = log(s/ME/ME);
		b = 2*ALPHAPI*(L-1.);
	}
	double operator()(double x)
  {
		return F_reg1(b,x)*sigma_final_vp_radcor(s*(1.-x),mt);
	}
};



class F_reg_x_b_sigma
{
	const double mt;
	const double s;
	double L;
	double b;
	double rb;
public:
	F_reg_x_b_sigma(double s_, double mt_)	: mt(mt_), s(s_)
  {
		L = log(s/ME/ME);
		b = 2*ALPHAPI*(L-1.);
		rb = 1./b;
	}
	double operator()(double y)
  {
		return F_reg_x_b(b,L)*sigma_final_vp_radcor(s*(1.-pow(y,rb)),mt);
	}
};



class F_reg_ln_x_sigma
{
	const double mt;
	const double s;
public:
	F_reg_ln_x_sigma(double s_, double mt_)	: mt(mt_), s(s_)
  {
	}
	double operator()(double y)
  {
		double x = y*y;	
		return F_reg_ln_x(x)*sigma_final_vp_radcor (s*(1.-x),mt);
	}
};


class F2_reg_x_sigma
{
	const double s;
	const double b;
	const double L;
	const double mt;
public:
	F2_reg_x_sigma(double s_, double b_, double L_,double mt_)	: s(s_), b(b_), L(L_),mt(mt_) {}
	double operator()(double y)
  {
		double x=exp(y);
		return F2_reg_x(x,s,b,L)*sigma_final_vp_radcor(s*(1.-x), mt);
	}
};
/*
class GFSreg2	{
	const double W;
	const double mtau;
	const double delta;
public:
	GFSreg2(double w, double mt, double dt)	: W(w), mtau(mt), delta(dt) {}
	double operator()(double *z)	{
		double eps = (1. - 4.* sq(mtau / z[1]));
//	double L = 2*log(z[1]/mtau);
		double b = 2*ALPHAPI*( 2*log(z[1]/mtau) - 1.);
		double x = z[0]*z[0]*eps;
		return F_reg_ln_x(x)*sigma_final_vp_radcor(z[1]*z[1]*(1.-x),mtau)*Gauss(W-z[1],delta)*sqrt(eps);
	}
};
*/
// Сечение с учетом всех радпоправок.
inline double sigma_all_radcor(double s, double mt, double prec)
{
	if ( NO_INI_RADCOR )
  {
		return sigma_final_vp_radcor(s,mt);
	}
	double eps = 1. - sq(2*mt)/s;	
	if (eps <= 0) return 0;
	//double xmax = min(eps,2*60/sqrt(s));
	double xmax = eps;
	//double xmax = 2*60/sqrt(s);
	double ymax;
	//double ymin;
	double integ=0, I = 0;
	double L = log(s/ME/ME);
	double b = 2*ALPHAPI*(L-1.);
		
	ymax = xmax;
	integ=b*ibn::dgaus(F_reg1_sigma(s,mt),0,ymax,prec*10);
	I+=integ;
	//cout << "f_reg1="<< integ << endl;
	
	ymax = pow(xmax,b);
	integ=ibn::dgaus(F_reg_x_b_sigma(s,mt),0,ymax,prec);
	I+=integ;
	//cout << "f_reg_x_b="<< integ << endl;
	
	ymax = sqrt(xmax);
	integ= b*b *ibn::dgaus(F_reg_ln_x_sigma(s,mt),0,ymax,prec)/4.;
	I+=integ;
	//cout << "f_reg_lnx="<< integ << endl;
	//return I;

	//Излучение электрон-позитронных пар с начального состояния.
	/*
	ymax = log(xmax);
	ymin = log(4*ME/sqrt(s));
	integ = sq(ALPHAPI)*ibn::dgaus( F2_reg_x_sigma(s,b,L,mt),ymin,ymax,prec*100);
	I+=integ;
	*/
	return I;
}


class XB
{
	const double rb;
	public:
		XB(double b): rb(b) {};
		double operator() (double y)	{ return pow(y,rb); }
};

class Y2
{
	public:
		Y2(void) {};
		double operator() (double y)	{ return y*y; }
};

class EXP
{
	public:
		EXP(void) {};
		double operator() (double y)	{ return exp(y); }
};
/*
inline double sigma_all_radcor2(double s, double mt, double prec)	{
	double eps = 1. - sq(2*mt)/s;	
	if (eps <= 0) return 0;
	double xmax = eps;
	double ymax;
	double ymin;
	double integ=0, I = 0;
	double L = log(s/ME/ME);
	double b = 2*ALPHAPI*(L-1.);
	SigmaFinVP sigm(mt);	
	FReg1 freg1(b);
	FRegxb fregxb(b,L);   XB xb(b);
	FRegLn fregln; 			  Y2 y2;
	FReg2x freg2x(s,b,L); EXP ex;
	
	integ=b*svertka(freg1,sigm,0,xmax,prec*10);
	I+=integ;

	integ = svertka_reg(fregxb,sigm,xb,0,pow(xmax,b),prec);
	I+=integ;

	integ = b*b/4. * svertka_reg(fregln,sigm,y2,0,sqrt(xmax),prec) ;
	I+=integ;

	//Излучение электрон-позитронных пар с начального состояния.
	integ = sq(ALPHAPI) * svertka_reg(freg2x,sigm,ex,log(4*ME/sqrt(s)),log(xmax),prec*100);
	I+=integ;

	return I;
}
*/
// Вычисление сечения с учетом разброса энергии в пучке.

static double INTEGRAL_PRECISION;
class Spread
{
	const double W;
	const double delta;
	const double mt;
	double prec;
	public:
	
	Spread(double W_, double delta_, double mt_, double p): W(W_), delta(delta_), mt(mt_)
  {
		prec = p;
	}
	double operator()(double w)
  {
			return exp(-sq(w-W)/(2*delta*delta))*sigma_all_radcor(w*w,mt,INTEGRAL_PRECISION);
		//if(w>=1777) return exp(-sq(w-W)/(2*delta*delta))*100 ; 
		//else return 0;
	}
};

inline  double sigma_total2(double W, double delta, double mt, double prec)
{
	double nsig=3;
	Spread theSpread(W,delta,mt,prec);

		 INTEGRAL_PRECISION = prec;
	//	double P[7] = { 0, 2.14, 3.03, 3.71, 4.29 , 4.8 , 5.25};
	double P[7] = { 0, 3.0, 4.8 , 5.25, 7};
	double result=0;
  if(W < 2*mt)
  {
    for(int i = 0; i < 2 ; i++)
    {
      result+=ibn::dgaus(theSpread,max(2*mt, W + delta * P[i]),   max(2*mt, W + delta * P[i+1]), INTEGRAL_PRECISION);
      //cerr << "result = " << result << ", min=" << W+delta*P[i] << ", max = " << W+delta*P[i+1] << ", P[i]=" <<P[i] <<  endl;
      //cerr <<  "W = "<< W << " E = " << W/2. << ", mt=" << mt << endl ; 
      if(INTEGRAL_PRECISION < 0.1) INTEGRAL_PRECISION*=100;
    }
    if (result ==0)
    { 
      result+=ibn::dgaus(theSpread,2*mt,   2*mt+(nsig-1)*delta, INTEGRAL_PRECISION);
      //cerr << " result = " << result << endl;
    } ;
    return result/delta/sqrt(2*M_PI);
  }
  else
  {
    INTEGRAL_PRECISION = prec;
    if(W > 2*mt+nsig*delta)
    {
      result=0;
      result += ibn::dgaus(theSpread,2*mt+nsig*delta, W+2*nsig*delta, prec);
      INTEGRAL_PRECISION*=100;
      result += ibn::dgaus(theSpread,2*mt, 2*mt + nsig*delta, INTEGRAL_PRECISION);
      //if (result ==0) { cerr << " result = 0!!!" << endl ; } ;
      return result/delta/sqrt(2*M_PI);
    }

    /*
       if( W > 2*mt + 2*nsig*delta)	{
       result += ibn::dgaus(theSpread,2*mt+2*nsig*delta, W+6*delta, INTEGRAL_PRECISION);
       INTEGRAL_PRECISION*=100;
       result += ibn::dgaus(theSpread,2*mt+nsig*delta, 2*mt + 2*nsig*delta , INTEGRAL_PRECISION);
       INTEGRAL_PRECISION*=100;
    //result += ibn::dgaus(theSpread,2*mt, 2*mt + nsig*delta, INTEGRAL_PRECISION);
    return result/delta/sqrt(2*M_PI);
    } */

    return 1./(sqrt(2*M_PI)*delta)*ibn::dgaus(theSpread,2*mt, W+2*nsig*delta, prec);
  }
  return 1./(sqrt(2*M_PI)*delta)*ibn::dgaus(theSpread,2*mt, W+2*nsig*delta, prec);
}







// Сечение с вычислением двойного интеграла 
// Мнемоника:
// GFS --- Gauss , F (фадинские поправки) , Sigma
// reg --- часть функции без особенностей
// xb  --- особенность x^{\betta -1}
// ln  --- логарифмическая особенность


inline double Sigma ( double s, double mt)
{
  double sig = sigma_tree_fc_vp_fsrc(s,mt);
  if ( std::isnan(sig) || std::isinf(sig) )
  {
    sig = -1;
    cout << "ERROR: sigma = nan || inf \n";
    exit(1);
  }
  return sig;
} 

double Gauss( double w, double delta)
{
	return exp(- sq(w/delta)/2. )/delta/sqrt( 2*M_PI );
}

class GFSreg
{
	const double W;
	const double mtau;
	const double delta;
public:
	GFSreg( double w, double mt, double dt)	: W(w), mtau(mt), delta(dt) {}
	double operator()(double *y)
  {
		double eps =  1. - 4.* sq(mtau / y[1]);
    if(eps<=0) return 0;
		double x = y[0]*eps;
		double b = 2*ALPHAPI*( 2*log(y[1]/ME) - 1.);
		double result =  eps*b*F_reg1(b,x)*Sigma(y[1]*y[1]*( 1.-x),mtau)*Gauss(W-y[1],delta);
    if(std::isnan(result)|| std::isinf(result))
    {
      cerr << "GFSreg: " << result << " W=" << W << " b=" << b << " x=" <<  x << " eps = " << eps << endl;
    }
    return result;
	}
	int dimension(void) { return 2; }
};

class GFSxb	{
	const double W;
	const double mtau;
	const double delta;
public:
	GFSxb( double w, double mt, double dt)	: W(w), mtau(mt), delta(dt) {}
	double operator()( double *z )
  {
		double eps = (1. - 4.* sq(mtau / z[1]));
    if(eps<=0) return 0;
		double L = 2*log(z[1]/ME);
		double b = 2*ALPHAPI*(L-1.);
		double x = pow(z[0],1./b)*eps;
	  double result  = pow(eps,b)*F_reg_x_b(b,L)*Sigma(z[1]*z[1]*( 1.-x),mtau)*Gauss(W-z[1],delta);
    if(std::isnan(result)|| std::isinf(result))
    {
      cerr << "GFSxb: " << result << " W=" << W << " b=" << b << " x=" <<  x << " eps = " << eps <<  " L=" << L << endl;
    }
    return result;
	}
	int dimension(void) { return 2; }
};

class GFSln	{
	const double W;
	const double mtau;
	const double delta;
public:
	GFSln(double w, double mt, double dt)	: W(w), mtau(mt), delta(dt) {}
	double operator()(double *z)
  {
		double eps = ( 1. - 4.* sq(mtau / z[1]));
    if(eps<=0) return 0;
		double b = 2*ALPHAPI*( 2*log(z[1]/ME) - 1.);
		double x = z[0]*z[0]*eps;
		double result= 0.75*b*b*F_reg_ln_x(x)*Sigma(z[1]*z[1]*(1.-x),mtau)*Gauss(W-z[1],delta)*sqrt(eps);
    if(std::isnan(result) || std::isinf(result))
    {
      cerr << "GFSln: " << result << " W=" << W << " b=" << b << " x=" <<  x << " eps = " << eps << endl;
    }
    return result;
	}
	int dimension(void) { return 2; }
};



inline double partial_sigma(double W, double  delta,double mt,double wmin,double wmax,double prec)	{
	if(wmin==wmax) return 0;
	double a[2] = { 0, wmin};
	double b[2] = { 1, wmax};
	double I=0;
	double relerr1,relerr2,relerr3;
	double integ=0;
	integ= ibn::dgaus( GFSreg(W,mt,delta), a, b, 100*prec, relerr1);
	I+=integ;

	integ= ibn::dgaus( GFSxb(W,mt,delta), a, b, prec, relerr2);
	I+=integ;
	
	integ= ibn::dgaus( GFSln(W,mt,delta), a, b, 100*prec, relerr3);
	I+=integ;
	return I;
}

inline double sigma_total(double W, double delta, double mt, double prec)	{
	double I=0;//,integ=0;
//	int nd=5;
	double range = 10*delta;
	//double range = 3*delta;
	I+=partial_sigma(W,delta, mt, max(2*mt,W-range), max(2*mt,W+range),prec);
	return I;
}

#endif
