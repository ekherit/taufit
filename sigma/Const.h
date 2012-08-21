#ifndef IBN_CONST_H
#define IBN_CONST_H

#include <cmath>

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

//mass of electron
//static const double ME=0.510998902;// mass of electron
static const double ME_PDG2011=0.510998910; //+-0.000000013 
static const double ME_PDG=ME_PDG2011;// mass of electron
static const double ME=ME_PDG;// mass of electron

//muon
static const double MMU_PDG2011=105.6583668; //+-0.0000038
static const double MMU_PDG = MMU_PDG2011;
static const double MMU = MMU_PDG;

//tau mass
static const double MTAU_PDG2002 = 1776.99; // mass of tau lepton PDG-2002
static const double MTAU_PDG2011=1776.82; //error 0.16; 
static const double DMTAU_PDG2011=0.16;
static const double MTAU = MTAU_PDG2011; 
static const double MTAU_PDG = MTAU_PDG2011; 
static const double DMTAU_PDG = DMTAU_PDG2011; 
static const double MTAU2 = MTAU*MTAU; // m_\tau^2

static const double ALPHA=1./137.03599976; //fine structure constant e^2/4pi
//static const double ALPHA=1./137.036; //fine structure constant e^2/4pi
static const double ALPHAPI=ALPHA/M_PI; //\alpha / \pi
static const double PIALPHA=M_PI*ALPHA; 


static const double GAMMA_E = 0.577215664901532861; //eiler number
static const double PI2 = M_PI * M_PI;  // \pi^2
static const double LOGPI = log(M_PI);
static const double LOGPI2 = log(M_PI)*0.5;

static const double LAMBDAE = 3.861592642e-11;//cm
static const double SIGMA_TOMSON = 0.665245854; //barn


//const double NBARN = sq(LAMBDAE*ME)*1e33; // ??????? 1/Mev^2 ? ?????????

const double SIGMA_CONST = 1e12*SIGMA_TOMSON * ME*ME/2.; // ??????? 1/Mev^2 ? ????????? ????? ???????????? ???????
const double SIGMA_CONST2 = 1e12*SIGMA_TOMSON*ME*ME/4.; // ??????? 1/Mev^2 ? ????????? ????? ???????????? ???????

// ????????? ?????????????? ???????
inline double sq(double x)	{ return x*x; }
inline double cb(double x)  { return x*x*x;}
inline double max(double x1, double x2)	{ return x1> x2 ? x1 : x2; }
inline double min(double x1, double x2)	{ return x1 < x2 ? x1 : x2; }
#endif
