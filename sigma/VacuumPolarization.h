#ifndef IBN_VACUUM_POLARIZATION_H
#define IBN_VACUUM_POLARIZATION_H

extern "C"	
{
#include <vp/vp.h>
}

#include "Const.h"
#include <iostream>
#include <vector>
#include <cstdlib>

#include "vplepton.h"
using namespace std;

static bool VP_LIB_INIT=false;

inline double Vp(double s, double mtau = MTAU)
{
  double S=s*1e-6;
  double ReH=0,ImH=0,ReL=0,ImL=0;
  double Re, Im;
  vp::lepton(s,&ReL,&ImL,mtau);  // Считаем лептонную часть поляризации
  if(!VP_LIB_INIT)
  {
    vp_init();
    VP_LIB_INIT=true;
  }
  vp_hadron(&S,&ReH,&ImH);  // Считаем адронную часть поляризации
  Re=ReL+ReH;
  Im=ImL+ImH;
  cout << "S = " << s << " ReL=" << ReL << " ImL=" << ImL << "  ReH=" << ReH << " ImH=" << ImH << endl;
  return 1./((1.-Re)*(1.-Re) + Im*Im);
}

inline double Vp_old(double S)
{
  S=S*1e-6;
  double ReH=0,ImH=0,ReL=0,ImL=0;
  double Re, Im;
  vp_lepton(&S,&ReL,&ImL);  // Считаем лептонную часть поляризации
  vp_hadron(&S,&ReH,&ImH);  // Считаем адронную часть поляризации
  Re=ReL+ReH;
  Im=ImL+ImH;
  return 1./((1.-Re)*(1.-Re) + Im*Im);
}

inline double Vp_hadron(double S)
{
  S=S*1e-6;
  double ReH=0,ImH=0;
  double Re, Im;
  vp_hadron(&S,&ReH,&ImH);  // Считаем адронную часть поляризации
  //std::cout << sqrt(S) << "  " << ReH << "    " << ImH << std::endl;
  Re=ReH;
  Im=ImH;
  //return 1./((1.-Re)*(1.-Re) + Im*Im);
  return 1./((1.-Re)*(1.-Re) + Im*Im);
}


inline double ReHadron( double S )
{
  S=S*1e-6;
  double Re, Im;
  vp_hadron(&S, &Re,&Im); 
  return Re;
}

inline double ImHadron( double S )
{
  S=S*1e-6;
  double Re, Im;
  vp_hadron(&S, &Re,&Im); 
  return Im;
}

inline double Vp_lepton(double S)
{
  S=S*1e-6;
  double ReL=0,ImL=0;
  double Re, Im;
  vp_lepton(&S,&ReL,&ImL);  // Считаем лептонную часть поляризации
  Re=ReL;
  Im=ImL;
  return 1./((1.-Re)*(1.-Re) + Im*Im);
}

inline double VP_tau(double S,double mt)
{
  double x = sqrt(S)/2.0 - mt;
  if(x<0) x = 0;
  return 1.0 + 3.31349e-2 - 1.21436e-4*x - 2.01719e-6*x*x - 3.64106e-8*x*x*x;
}

inline double VP_tau2(double S,double mt)
{
  double x = sqrt(S)/2.0 - mt;
  if(x<0) x = 0;
  //if( x > 23 ) return Vp(S);
  return 1.0 + 3.31163e-2 - 1.06988e-4*x - 4.56869e-6*x*x + 1.22004e-7*x*x*x - 3.21668e-9*x*x*x*x;
}

/* Подгонка в диапазоне от 0 до + 50 MeV полиномом девятой степени
   здесь все в куче и лептонная и адронная часть. 
   FCN=64.7929 FROM MIGRAD    STATUS=FAILED        822 CALLS         823 TOTAL
   EDM=2.24267e-12    STRATEGY= 1      ERR MATRIX NOT POS-DEF
   EXT PARAMETER                APPROXIMATE        STEP         FIRST   
   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           3.53393e-02   3.06498e-07   1.68511e-08   5.00381e-01
   2  p1          -1.33187e-04   2.35831e-08   6.35085e-11  -2.28261e+00
   3  p2           3.83315e-06   7.22536e-10   1.82779e-12  -8.64102e+02
   4  p3          -1.19355e-06   1.77034e-11   5.69130e-13  -4.48394e+04
   5  p4           1.13451e-07   4.13114e-13   5.40977e-14  -3.47515e+06
   6  p5          -6.23188e-09   9.45323e-15   2.97159e-15  -2.55502e+07
   7  p6           2.03424e-10   2.11637e-16   9.70002e-17  -1.12574e+10
   8  p7          -3.97272e-12   4.59654e-18   1.89434e-18  -3.63634e+10
   9  p8           4.29167e-14   9.57794e-20   2.04643e-20  -1.58101e+13
   10  p9          -2.00078e-16   1.87715e-21   9.54044e-23  -4.81777e+14
   */

/*
   Подгонка в области пси штриха в подложке 1836 до 1839 МэВ
   EXT PARAMETER                APPROXIMATE        STEP         FIRST   
   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0          -9.70098e+03   2.63572e+00   2.91029e-01  -1.98432e-10
   2  p1           1.05745e+01   1.43446e-03   3.17234e-04   6.31486e-07
   3  p2          -2.88166e-03   7.80587e-07   8.64498e-08  -6.11931e-04
   Картинка лежит в result/tau06/vp-psi2-zoom-fit.eps

*/
inline double VP_pol9( double S, double mt)
{
  double x = sqrt(S)/2.0 - mt;
  if ( x < 0 ) x = 0;
  if ( x > 47 )
  {
    if( x+mt > 1836 && x+mt < 1839 )
    {
      double p[3] =
      { 
        -9.70098e+03,
        1.05745e+01,
        -2.88166e-03
      };
      double sum=0,pw=1;
      for(int i=0;i<3;i++)
      {
        sum+=pw*p[i];
        pw*=(x+mt);
      }
      return sum;
    }
    if( x+mt > 1887 && x+mt < 1890)
    {
      return 1+0.54438;
    }
    std::cout << "Using vp lib for  E=" << x+mt << std::endl;
    return Vp(S);
  }
  double sum=1, pw=1;
  double p[10] = 
  {  3.53393e-02, 
    -1.33187e-04, 
    3.83315e-06, 
    -1.19355e-06, 
    1.13451e-07, 
    -6.23188e-09, 
    2.03424e-10, 
    -3.97272e-12, 
    4.29167e-14, 
    -2.00078e-16 };
  for ( int i = 0; i< 10; i++)
  {
    sum+=pw*p[i];
    pw*=x;
  }
  return sum;
}


inline double Vpol(double g2)
{
  double S=g2*1e-6*MTAU2;
  double ReH,ImH,ReL,ImL;
  double Re, Im;
  vp_lepton(&S,&ReL,&ImL);  // Считаем лептонную часть поляризации
  vp_hadron(&S,&ReH,&ImH);  // Считаем адронную часть поляризации
  Re=ReL+ReH;
  Im=ImL+ImH;
  return 1./((1.-Re)*(1.-Re) + Im*Im);
}


inline double vp_ROOT(double *x, double *p)
{
  double E=*x,S;
  double ReH,ImH,ReL,ImL;
  //double E_Re,E_Im;
  double Re,Im;
  double P;

  //double D_Re,D_Im,D_ReIm;

  S=E*E;
  if(E<0) S=-S;

  vp_lepton(&S,&ReL,&ImL);  // Считаем лептонную часть поляризации
  vp_hadron(&S,&ReH,&ImH);  // Считаем адронную часть поляризации
  //vp_hadron_err(&S,&D_Re,&D_Im,&D_ReIm);
  //E_Re=sqrt(D_Re);
  //E_Im=sqrt(D_Im);

  Re=ReL+ReH;
  Im=ImL+ImH;
  P= (1-Re)*(1-Re) + Im*Im;
  return 1/P-1;
}

inline double VP1777( double S, double mt)
{
  double x = sqrt(S)/2.0 - mt;
  if ( x > 47 && x <0 )
  { 
    std::cerr << "ERROR: VP1777  bad energy E=" << x+mt << " is out of range\n";
    exit(1);
  }
  double sum=1, pw=1;
  double p[10] = {  
    3.53393e-02, 
    -1.33187e-04, 
    3.83315e-06, 
    -1.19355e-06, 
    1.13451e-07, 
    -6.23188e-09, 
    2.03424e-10, 
    -3.97272e-12, 
    4.29167e-14, 
    -2.00078e-16 
  };
  for ( int i = 0; i< 10; i++)
  {
    sum+=pw*p[i];
    pw*=x;
  }
  return sum;
}

class VPtab
{
  std::vector <double> buf;
  bool isinit;
  public:
  double Emin, Emax,dE;
  VPtab(double emin, double emax,double de)
  {
    Emin=emin;
    Emax=emax;
    dE=de;
    isinit =false;
  }
  void Init(void)
  {
    std::cout << "Init vp library ... ";
    vp_init();
    cout << " OK\n";
    std::cout << "Init VPtab table " << "("<<Emin<<","<<Emax<<") dE="<<dE<<" ... " << flush;
    buf.resize(size_t((Emax-Emin)/dE));
    double E,vp;
    for(size_t i=0;i<buf.size();i++)
    {
      E=Emin+i*dE;
      vp=Vp(4*E*E);
      if(std::isnan(vp) || std::isinf(vp))
      {
        cerr << " ERROR  VPtab init bad vp value: " << vp << " i="<< i << " E=" << E << endl;
        exit(1);
      }
      buf[i]=vp;
    }
    std::cout << " OK\n";
    isinit = true;
  }

  double operator()(double S,double mt)
  {
    if(!isinit) Init();
    double E=sqrt(S/4);
    if(E>=Emin && E<Emax)
    {
      double x=(E-Emin)/dE;
      size_t i=size_t(x);
      double vp;
      /*
         std::cout << "Using new fucntion\n";
         std::cout << "E="<<E << " S=" << S << std::endl;
         std::cout << "bufsize="<<buf.size() << std::endl;
         i=
         std::cout << "i="<< i<<std::endl;
         */
      if(i<buf.size()-1)
      {
        vp=buf[i]+ (x-i)*(buf[i+1]-buf[i]);
      }
      else
      {
        return buf[i];
      }
      return vp;
    }
    std::cerr << "ERROR: VPtab energy is out of range E=" << E << "\n";
    exit(1);
  }
};

// hadron tab
class VPHtab
{
  std::vector <double> RE;
  std::vector <double> IM;
  bool isinit;
  public:
  double Emin, Emax,dE;
  double Smin, Smax,dS;
  unsigned N;

  VPHtab(double E, double Erange, double de)
  {
    Emin=E-Erange/2.;
    Emax=E+Erange/2.;
    dE=de;
    Smin = 4*Emin*Emin;
    Smax = 4*Emax*Emax;
    dS   = 4*(Emin+Emax)*dE;
    isinit =false;
  }
  void Init(void)
  {
    std::cout << "Init vp hadron tabular ... ";
    std::cout << "("<<Emin<<","<<Emax<<") dE="<<dE<<" ... " << flush;
    if(!VP_LIB_INIT) vp_init();
    VP_LIB_INIT=true;
    cout << " OK\n";
    N =  (Smax-Smin)/dS;
    RE.resize(N);
    IM.resize(N);
    double S,s;
    double re,im;
    for(size_t i=0;i<N;i++)
    {
      S = Smin + i*dS;
      s = S*1e-6;
      vp_hadron(&s,&re,&im);
      RE[i]=re;
      IM[i]=im;
      cout << sqrt(S)/2. << " " << re << " " << im << endl;
    }
    std::cout << " OK\n";
    isinit = true;
  }

  double vp(double S, double & re, double & im)
  {
    if(!isinit) Init();
    double x=(S-Smin)/dS;
    int i=int(x);
    if(i>=0 && i<N-1)
    {
      re=RE[i]+ (x-i)*(RE[i+1]-RE[i]);
      im=IM[i]+ (x-i)*(IM[i+1]-IM[i]);
    }
    else 
    {
      double s=S*1e-6;
      cout << "out of range" << endl;
      vp_hadron(&s,&re,&im);
    }
  }
};

static VPtab VPT(1820,1910,0.001);  //расширил диапазон и добавил точности.
static VPHtab VPHT(MTAU,100,0.01);  //расширил диапазон и добавил точности.

//Лептонную часть считаю сам, а адронную по таблице.
inline double Vp2(double s, double mtau = MTAU)
{
  double S=s*1e-6;
  double ReH=0,ImH=0,ReL=0,ImL=0;
  double Re, Im;
  vp::lepton(s,&ReL,&ImL,mtau);  // Считаем лептонную часть поляризации
  VPHT.vp(s, ReH, ImH);
  Re=ReL+ReH;
  Im=ImL+ImH;
  //cout << "S = " << s << " ReL=" << ReL << " ImL=" << ImL << "  ReH=" << ReH << " ImH=" << ImH << endl;
  return 1./((1.-Re)*(1.-Re) + Im*Im);
}

static inline double VP(double S, double mt=MTAU)
{
  if(S<0 || std::isnan(S) || std::isinf(S))
  {
    return 1e300;
  }
  double E = sqrt(S)/2.0;
  if(E<1777+47) return VP1777(S,mt);
  if(E>=VPT.Emin && E<VPT.Emax)    return VPT(S,mt);
  std::cerr << "Using VP\n";
  std::cerr << S << " " << E << " " << "mt=" << mt << std::endl;
  return Vp(S);
}


#endif
