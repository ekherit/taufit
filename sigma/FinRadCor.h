#ifndef IBN_FIN_RAD_COR_H
#define IBN_FIN_RAD_COR_H
#include "Const.h"
#include <iostream>
#include <complex>
#include <vector>
//#include"polylog.h"
#include "Integral.h"
/* Final states corrections... 
 * \sigma_1 = SIGMA_CONST b(3-b)/2 * Fc*Fr
 * 
 * M.B.Voloshin The onset of e^+e^- \rightarrow \tau^+ \tau^- at threshold revisitted
 * hep-ph/012207
 */


inline double Fc(double v); // Кулоновское взаимодейсвтие в конечном состоянии 
inline double vFc(double v); // Fc * v
inline double Fr(double v);  //радиационные поправки (петлевые и всякая хрень)

// Некоторые вспомогательные функции 

static const double PRECISION=1e-10;
static const int LI2MAXORDER = 2000;
inline double Li2(double x);
inline double Li2_array(double x);
inline double Li2_int(double x);
inline double Li2_int_reg(double x);
inline double Li2_int2( double x1, double x2);
inline double Li2_int2_reg( double x1, double x2);
inline double S(double b);
inline double h(double v, double mt = MTAU);


class Polylog;
using namespace std;
// Colomb intaraction correction in final states (Voloshin) 
inline double Fc(double v)	{
    double z=PIALPHA/v;
    return z/(1-exp(-z));
}

// Fc*v function
inline double vFc(double v)	{
    if(v==0) return PIALPHA;
    return PIALPHA/(1. - exp(-PIALPHA/v));
}


inline double Fr(double v)	{
    double fr = 1. + ALPHAPI* (S(v) + 2./3. * h(v) ) ;
    return fr;
}



inline double S(double v)	{
    if(v < 0.0001) return -4;
    return 1./v * ( - PI2/2. + (1.+v*v)*( PI2/6.0 + log((1.+v)/2.0)*log((1.+v)/(1.-v)) + 
		2. * ( Li2((1.-v)/(1.+v)) + Li2((1.+v)/2.) - Li2((1.-v)/2.) -
		    2.*Li2(v) ) + Li2(v*v) ) +
	    log((1.+v)/(1.-v)) * ( 11./8.*(1.+v*v) - 3.*v + sq(v*v)/(3-v*v)/2.) +
	    6.*v*log((1.+v)/2.) - 4.*v*log(v) + 3./4.*v*(5.-3.*v*v)/(3.-v*v) );
}

/*inline double S(double v)	{
  return 1./v * ( (1.+v*v)*( -PI2/3.0 + log((1.+v)/2.0)*log((1.+v)/(1.-v)) + 
  2. * ( Li2((1.-v)/(1.+v)) + 
  Li2_int2((1.-v)/2.,(1.+v)/2.)
  - 2.*Li2(v) ) + Li2(v*v) ) +
  log((1.+v)/(1.-v)) * ( 11./8.*(1.+v*v) - 3.*v + sq(v*v)/(3-v*v)/2.) +
  6.*v*log((1.+v)/2.) - 4.*v*log(v) + 3./4.*v*(5.-3.*v*v)/(3.-v*v) );
  }*/


inline double Li2(double x)	
{
    /*
       double sum=0;
       for(int i = 1; i<max; i++)	{
       sum+=pow(x,i)/(i*i);
       }
       return sum;
     */
    //return Li2_array(x);
    return Li2_int_reg(x);
}
/*
   inline double Li2_array(double x)	{
   if(x==1.) return PI2/6.;
   int n1 = (int)floor(x*POLYLOG_ARRAY_SIZE);
   int n2 = n1+1;
   return PolylogArray[n1]+(PolylogArray[n2]-PolylogArray[n1])*(x-double(n1)/POLYLOG_ARRAY_SIZE);
   }
 */
class Polylog	{
    public:
	Polylog(void){}
	double operator() (double x)	{
	    return -log(1.-x)/x;
	}	
};

class Polylog_reg	{
    public:
	Polylog_reg(void){}
	double operator() (double y)	{
	    return -4.*y/(1.-y*y)*log(y);
	}	
};

inline double Li2_int(double x)	{
    if(x==1.) return PI2/6.;
    if(x==0) return 0;
    return dgaus(Polylog(),0,x,PRECISION);
}

inline double Li2_int_reg(double x)	{
    if(x==1.) return PI2/6.;
    if(x==0) return 0;
    return dgaus(Polylog_reg(),sqrt(1.-x),1.,PRECISION);
}

inline double Li2_int2_reg( double x1, double x2)	{
    return dgaus(Polylog_reg(),sqrt(1.-x2),sqrt(1.-x1),PRECISION);
}

inline double Li2_int2( double x1, double x2)	{
    return dgaus(Polylog(),x1,x2,PRECISION);
}

class h_sub1 {
    typedef complex<double> comp_t;
    const double zv;
    const double t;
    const double lambda;
    public:
    h_sub1( double zv_, double t_, double l_) : zv(zv_), t(t_), lambda(l_) {}
    comp_t operator()( double x ) 	{
	return pow(t + zv/x , comp_t(-1., lambda)) / pow( t+1. + zv/x, comp_t(1.,lambda) ) *
	    ( 1. +0.5 * sq(x)) * sqrt( 1. - x*x) /x/x;
    }
};

class h_sub2 {
    typedef complex<double> comp_t;
    const double zv;
    const double lambda;
    public:	
    h_sub2( double zv_, double l_) : zv(zv_), lambda(l_) {}
    comp_t operator()( double etta)	{
	double t = etta /(1.-etta);
	comp_t tmp = pow( (1+t)/t, comp_t(0,lambda)) * dgaus_comp( h_sub1(zv,t,lambda) , 0, 1,PRECISION)/ sq(1.-etta);
	//cerr << tmp << endl;
	return tmp;
    }
};

inline double h_int(double v, double mt=MTAU)	{
    //Попытка вычислить h прямым интегрированием.
    double lambda = ALPHA/2./v;
    double zv = ME/mt/v;
    complex<double> tmp=dgaus_comp(h_sub2(zv,lambda) , 0, 1, PRECISION); 
    return -2.*lambda * tmp.imag();
}
inline double h(double v, double mt )	{
    double H0 =  log(mt*ALPHA/ME) + GAMMA_E + 1./6.;
    if(v==0) return H0;
    double z = PIALPHA /v;
    double h0 = (1. - (1.+z)*exp(-z))/(1. - exp(-z));
    double b;

    if( v <= 0.01)  {
	b = 0.998929  - 7.71335 *v + 1033.43*v*v;
	return b*h0*H0;
    } 
    if( v <= 0.04)	{
	b = 0.927085 + 9.46920 * v  - 31.2432*v*v;
	return b*h0*H0;
    }
    if(v>0.04) return h0*(log(2*mt*v/ME) - 5./6.);
    //				else return h0*(log(mt*ALPHA/ME) + GAMMA_E +1./6.); 
    return 0;
}



/* this part creates tables for final state radiative corrections */


class FSRC	{
    double precision;
    /* 
     * 	maximum and minimal meaning of tau velocity
     * 	belong this value direct radiation corrections will be performed
     */
    double vmax, vmin; 
    double dv; //elementary step in velocity (it depends on precision )
    vector <double> rad_cor_array;
    double vel(unsigned i )	{	return (vmin  + i*dv); }
    public:
    FSRC(void){};
    FSRC( double p ,double  mx = 0.4, double mn = 0 ) 	{
	Init(p,mx,mn);		
    }	
    void Init( double p,double  mx = 0.3,double  mn = 0) {
	cout << "Init FSRC table v=("<<mn<<","<<mx<<") ..." << flush ;
	vmin = mn;
	vmax = mx;
	if(vmin <  0 ) { 
	    cerr << "error: FinRadCor class: vmin below zero." << endl;
	    exit(1);
	}
	if(vmax <= vmin) {
	    cerr << "error: FinRadCor class: vmax <= vmin in Fin." << endl;
	    exit(1);
	}
	if(p==0) {
	    cerr << "error: FinRadCor class: precision equal zero." << endl;
	    exit(1);
	}
	precision = p;
	unsigned N = unsigned (1./precision);
	double v;
	dv = (vmax-vmin)/N;
	rad_cor_array.resize(N+1);
	for ( unsigned i = 0 ; i< rad_cor_array.size() ; i++)	{
	    v= vel(i);
	    rad_cor_array[i] = Fr(v);
	    //cout << i << ": " << v << ": " << rad_cor_array[i] << endl;

	}
	std::cout << " OK\n";
    }
    double operator()(double v)	{
	if(v >= vmin && v <=vmax) {
	    double x=(v-vmin)/dv;
	    size_t i=size_t(x);
	    double r;
	    if(i<rad_cor_array.size()-1)    {
		r=rad_cor_array[i] + (x - double(i))*(rad_cor_array[i+1]-rad_cor_array[i]);
	    }
	    else {
		r=rad_cor_array[i];
	    }
	    return  r;
	} else {
	    return vFc(v)*Fr(v);
	}
    }
};

class FinRadCor	{
    double precision;
    /* 
     * 	maximum and minimal meaning of tau velocity
     * 	belong this value direct radiation corrections will be performed
     */
    double vmax, vmin; 
    double dv; //elementary step in velocity (it depends on precision )
    vector <double> rad_cor_array;
    double vel(unsigned i )	{	return (vmin  + i*dv); }
    public:
    FinRadCor(void){};
    FinRadCor( double p ,double  mx = 0.4, double mn = 0 ) 	{
	Init(p,mx,mn);		
    }	
    void Init( double p,double  mx = 0.3,double  mn = 0) {
	vmin = mn;
	vmax = mx;
	if(vmin <  0 ) { 
	    cerr << "error: FinRadCor class: vmin below zero." << endl;
	    exit(1);
	}
	if(vmax <= vmin) {
	    cerr << "error: FinRadCor class: vmax <= vmin in Fin." << endl;
	    exit(1);
	}
	if(p==0) {
	    cerr << "error: FinRadCor class: precision equal zero." << endl;
	    exit(1);
	}
	precision = p;
	unsigned N = unsigned (1./precision);
	double v;
	dv = (vmax-vmin)/N;
	rad_cor_array.resize(N+1);
	for ( unsigned i = 0 ; i< rad_cor_array.size() ; i++)	{
	    v= vel(i);
	    rad_cor_array[i] = vFc(v)*Fr(v);
	    //cout << i << ": " << v << ": " << rad_cor_array[i] << endl;

	}
    }
    double operator()(double v)	{
	if(v >= vmin && v <=vmax) {
	    return  rad_cor_array[  unsigned ( (v - vmin)/dv ) ];
	} else {
	    return vFc(v)*Fr(v);
	}
    }
};



//static FinRadCor FRC(1e-3,0.2,0.0); 

//static FSRC fsrc(1e-5,0.2,0.0); 
//static FSRC fsrc(1e-5,0.3,0.0);  убрал, так как точка 1888 не попадает.
static FSRC fsrc(1e-6,0.4,0.0);  //расширил диапазон и добавил точности.
// 0.4 --- до 1940 Mev
// 0.3 - 1862
// 0.2   1813
// 0.1   1785

// this functions are used by root TF1
double S(double * x, double *p)	{
    return S(*x);
} 

double h(double *x, double *p)	{
    return h(*x);
}

double vFc(double *x, double *p)	{
    return vFc(*x);
}
double Fr(double *x, double *p)	{
    return Fr(*x)-1.;
}
double Fc(double *x, double *p)	{
    return Fc(*x);
}
double Li2(double *x, double *p)	{
    //return Li2_array(*x);
    return Li2_int_reg(*x);
}
/*
double frc(double *x, double *p)	{
    return FRC(*x);
}
*/
#endif
