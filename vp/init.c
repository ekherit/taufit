#include "inc.h"


static Point CreatePoint2_R(double e,double r,double er)
 {
  double s=e*e;
  double aux=s*eemm(s);
  return CreatePoint(s,aux*r,aux*er);
 }

static Point CreatePoint2_MeVR(double e,double r,double er)
 {
  return CreatePoint2_R(0.001*e,r,er);
 }

static Point CreatePoint2_nb(double e,double sigma,double esigma)
 {
  double s =e*e;
  double y =s*(nb2GeV*sigma);
  double ey=s*(nb2GeV*esigma);

  return CreatePoint(s,y,ey);
 }

static Point CreatePoint2_MeVnb(double e,double sigma,double esigma)
 {
  return CreatePoint2_nb(0.001*e,sigma,esigma);
 }


static void AddToChain_R(Point V,double e,double r,double er)
 {
  double s=e*e;
  double aux=s*eemm(s);
  AddPoint(V, s, r*aux, er*aux);
 }

static void AddToChain_R_percent(Point V,double e,double r,double er)
 {
  double s=e*e;
  double aux=s*eemm(s);
  AddPoint(V, s, r*aux, (r*er/100)*aux);
 }

static void AddToChain_MeVR(Point V,double e,double r,double er)
 {
  AddToChain_R(V, 0.001*e, r, er);
 }

static void AddToChain_nb(Point V,double e,double sigma,double esigma)
 {
  double s =e*e;
  double y =s*(nb2GeV*sigma);
  double ey=s*(nb2GeV*esigma);
  
  AddPoint(V, s, y, ey);
 }

static void AddToChain_MeVnb(Point V,double e,double sigma,double esigma)
 {
  AddToChain_nb(V,0.001*e,sigma,esigma);
 }

//--------------------------------- Самая больша и главная функция 
void VP_Init()
 {
   Point Vrho,Vhigh,Vkpkm,Vkskl,V4pic,V4pin,Vphi,V5pi,V3pi_low,V3pi_high;

//------------------------------------------------------------------

/*-------------------------------------------------------------------
 Narrow resonances -> omega phi J/psi psi(2S) Y(1S) Y(2S) Y(3S) Y(4S)
-------------------------------------------------------------------*/
   Mr[0]   = 0.78194;
   Gr[0]   = 8.41e-3;
   Beer[0] = 7.07e-5;

   Mr[1]   = 1.019483;
   Gr[1]   = 4.28e-3;
   Beer[1] = 2.91e-4;

   Mr[2]   = 3.09688;
   Gr[2]   = 0.087e-3;
   Beer[2] = 6.02e-2;

   Mr[3]   = 3.686;
   Gr[3]   = 0.277e-3;
   Beer[3] = 8.5e-3;

   Mr[4]   = 9.4604;
   Gr[4]   = 0.0525e-3;
   Beer[4] = 2.52e-2;

   Mr[5]   = 10.02330;
   Gr[5]   = 0.044e-3;
   Beer[5] = 1.18e-2;

   Mr[6]   = 10.3553;
   Gr[6]   = 0.0263e-3;
   Beer[6] = 0.476e-6;

   Mr[7]   = 10.58;
   Gr[7]   = 10.e-3;
   Beer[7] = 2.8e-5;


/*-----------------------------------------------------------------------
 Континуум
-----------------------------------------------------------------------*/


//----------------------------------------------------------- \pi\pi
        Vrho=CreatePoint2_nb(treshold2,0.,0.);


//-------------------------------------------------------- Оля
        AddToChain_MeVnb(Vrho,	400.,	102.,	12.);
//-------------------------------------------------------- CMD
        AddToChain_MeVnb(Vrho,	360.,	78.9,	6.8);
        AddToChain_MeVnb(Vrho,	380.,	100.2,	9.4);
        AddToChain_MeVnb(Vrho,	410.,	102.1,	8.1);
        AddToChain_MeVnb(Vrho,	430.,	126.3,	12.4);
        AddToChain_MeVnb(Vrho,	438.,	143.3,	11.9);
        AddToChain_MeVnb(Vrho,	470.,	142.4,	9.2);
        AddToChain_MeVnb(Vrho,	540.,	208.1,	14.0);
        AddToChain_MeVnb(Vrho,	580.,	245.1,	16.5);
        AddToChain_MeVnb(Vrho,	620.,	347.4,	20.1);
//-------------------------------------------------------- \pi\pi (Ваня Логашенко)
        AddToChain_MeVnb(Vrho, 610.50,  327.8, 46.2 );
        AddToChain_MeVnb(Vrho, 620.50,  390.6, 29.1 );
        AddToChain_MeVnb(Vrho, 630.50,  426.9, 28.6 );
        AddToChain_MeVnb(Vrho, 640.51,  426.2, 27.9 ); 
        AddToChain_MeVnb(Vrho, 650.49,  480.3, 31.6 ); 
        AddToChain_MeVnb(Vrho, 660.50,  498.8, 27.8 ); 
        AddToChain_MeVnb(Vrho, 670.50,  560.5, 30.0 ); 
        AddToChain_MeVnb(Vrho, 680.59,  697.2, 32.3 ); 
        AddToChain_MeVnb(Vrho, 690.43,  718.6, 25.5 ); 
        AddToChain_MeVnb(Vrho, 700.52,  818.1, 21.2 ); 
        AddToChain_MeVnb(Vrho, 710.47,  926.1, 33.4 ); 
        AddToChain_MeVnb(Vrho, 720.25, 1048.1, 27.0 ); 
        AddToChain_MeVnb(Vrho, 730.24, 1107.2, 38.1 ); 
        AddToChain_MeVnb(Vrho, 740.20, 1191.9, 37.1 ); 
        AddToChain_MeVnb(Vrho, 750.28, 1302.9, 36.2 ); 
        AddToChain_MeVnb(Vrho, 760.18, 1306.9, 36.5 ); 
        AddToChain_MeVnb(Vrho, 764.17, 1289.6, 30.6 ); 
        AddToChain_MeVnb(Vrho, 770.11, 1304.6, 33.5 ); 
        AddToChain_MeVnb(Vrho, 774.38, 1265.7, 33.7 ); 
        AddToChain_MeVnb(Vrho, 778.17, 1321.6, 37.7 ); 
        AddToChain_MeVnb(Vrho, 780.17, 1234.9, 34.8 ); 
        AddToChain_MeVnb(Vrho, 782.23, 1026.6, 18.5 ); 
        AddToChain_MeVnb(Vrho, 784.24,  935.9, 26.5 );
        AddToChain_MeVnb(Vrho, 786.04,  809.5, 28.9 );
        AddToChain_MeVnb(Vrho, 790.10,  859.3, 30.0 );
        AddToChain_MeVnb(Vrho, 794.14,  799.6, 22.7 );
        AddToChain_MeVnb(Vrho, 800.02,  782.7, 18.1 );
        AddToChain_MeVnb(Vrho, 810.14,  669.7, 14.4 );
        AddToChain_MeVnb(Vrho, 820.02,  619.0, 20.0 );
        AddToChain_MeVnb(Vrho, 829.97,  517.5, 18.5 );
        AddToChain_MeVnb(Vrho, 839.10,  415.4, 17.3 );
        AddToChain_MeVnb(Vrho, 849.24,  344.8, 17.0 );
        AddToChain_MeVnb(Vrho, 859.60,  345.7, 15.9 );
        AddToChain_MeVnb(Vrho, 869.50,  260.7, 10.7 );
        AddToChain_MeVnb(Vrho, 879.84,  233.6, 18.1 );
        AddToChain_MeVnb(Vrho, 889.72,  191.1,  7.7 );
        AddToChain_MeVnb(Vrho, 900.04,  172.1,  6.4 );
        AddToChain_MeVnb(Vrho, 910.02,  149.1,  6.8 );
        AddToChain_MeVnb(Vrho, 919.56,  129.7,  6.5 );
        AddToChain_MeVnb(Vrho, 930.11,  120.3,  7.7 );
        AddToChain_MeVnb(Vrho, 942.19,  106.8,  5.1 );
        AddToChain_MeVnb(Vrho, 951.84,   92.9,  4.8 );
        AddToChain_MeVnb(Vrho, 961.52,   86.2,  4.6 );

   
//-------------------------------------------------------- \pi\pi (Федя Игнатов)
        AddToChain_MeVnb(Vrho, 2*490.0,	73.32,	3.97);
        AddToChain_MeVnb(Vrho, 2*520.0,	45.35,	3.01);
        AddToChain_MeVnb(Vrho, 2*525.0,	39.67,	2.31);
        AddToChain_MeVnb(Vrho, 2*530.0,	37.41,	3.34);
        AddToChain_MeVnb(Vrho, 2*535.0,	39.72,	2.66);
        AddToChain_MeVnb(Vrho, 2*540.0,	30.59,	2.44);
        AddToChain_MeVnb(Vrho, 2*545.0,	34.51,	2.25);
        AddToChain_MeVnb(Vrho, 2*550.0,	27.81,	2.22);
        AddToChain_MeVnb(Vrho, 2*555.0,	26.31,	2.08);
        AddToChain_MeVnb(Vrho, 2*560.0,	25.50,	3.26);
        AddToChain_MeVnb(Vrho, 2*565.0,	26.08,	1.85);
        AddToChain_MeVnb(Vrho, 2*570.0,	22.76,	1.55);
        AddToChain_MeVnb(Vrho, 2*575.0,	21.55,	2.38);
        AddToChain_MeVnb(Vrho, 2*580.0,	19.98,	1.53);
        AddToChain_MeVnb(Vrho, 2*585.0,	17.75,	1.81);
        AddToChain_MeVnb(Vrho, 2*590.0,	16.97,	1.46);
        AddToChain_MeVnb(Vrho, 2*595.0,	15.24,	1.43);
        AddToChain_MeVnb(Vrho, 2*600.0,	15.47,	1.24);
        AddToChain_MeVnb(Vrho, 2*605.7,	15.30,	1.32);
        AddToChain_MeVnb(Vrho, 2*610.0,	15.12,	1.18);
        AddToChain_MeVnb(Vrho, 2*615.0,	14.33,	1.36);
        AddToChain_MeVnb(Vrho, 2*620.0,	13.77,	1.26);
        AddToChain_MeVnb(Vrho, 2*625.0,	11.09,	0.95);
        AddToChain_MeVnb(Vrho, 2*630.0,	10.45,	0.93);
        AddToChain_MeVnb(Vrho, 2*635.0,	9.41,	0.72);
        AddToChain_MeVnb(Vrho, 2*640.0,	9.71,	0.80);
        AddToChain_MeVnb(Vrho, 2*645.0,	8.34,	0.65);
        AddToChain_MeVnb(Vrho, 2*650.0,	9.29,	0.73);
        AddToChain_MeVnb(Vrho, 2*655.0,	7.79,	0.83);
        AddToChain_MeVnb(Vrho, 2*660.0,	6.24,	0.75);
        AddToChain_MeVnb(Vrho, 2*665.0,	4.93,	0.46);
        AddToChain_MeVnb(Vrho, 2*670.0,	6.53,	0.84);
        AddToChain_MeVnb(Vrho, 2*675.0,	4.07,	0.57);
        AddToChain_MeVnb(Vrho, 2*680.0,	4.83,	0.56);
        AddToChain_MeVnb(Vrho, 2*685.0,	5.23,	1.15);
        AddToChain_MeVnb(Vrho, 2*690.0,	4.06,	0.64);


//-------------------------------------------------------- R выше 1.4
   	Vhigh=CreatePoint2_R(1.3801,1.68,	0.03);
   	
	AddToChain_R(Vhigh,1.4025,	1.68,	0.03);
	AddToChain_R(Vhigh,1.417,	1.66,	0.07);
	AddToChain_R(Vhigh,1.4385,2.01,	0.12);
	AddToChain_R(Vhigh,1.4585,1.78,	0.09);
	AddToChain_R(Vhigh,1.4795,1.95,	0.03);
	AddToChain_R(Vhigh,1.5045,1.87,	0.03);
	AddToChain_R(Vhigh,1.5345,1.66,	0.11);
	AddToChain_R(Vhigh,1.5605,2.23,	0.11);
	AddToChain_R(Vhigh,1.581,	2.04,	0.03);
	AddToChain_R(Vhigh,1.6005,2.2,	0.03);
	AddToChain_R(Vhigh,1.62,	2.25,	0.05);
	AddToChain_R(Vhigh,1.639,	2.16,	0.03);
	AddToChain_R(Vhigh,1.659,	2.28,	0.03);
	AddToChain_R(Vhigh,1.6805,2.25,	0.04);
	AddToChain_R(Vhigh,1.7005,2.17,	0.03);
	AddToChain_R(Vhigh,1.7195,2.16,	0.11);
	AddToChain_R(Vhigh,1.7385,2,	0.03);
	AddToChain_R(Vhigh,1.759,	2.09,	0.03);
	AddToChain_R(Vhigh,1.7805,2.12,	0.05);
	AddToChain_R(Vhigh,1.8005,2.34,	0.08);
	AddToChain_R(Vhigh,1.8195,2.09,	0.03);
	AddToChain_R(Vhigh,1.8385,1.89,	0.05);
	AddToChain_R(Vhigh,1.859,	1.97,	0.07);
	AddToChain_R(Vhigh,1.8805,2.09,	0.09);
	AddToChain_R(Vhigh,1.9005,2.06,	0.07);
	AddToChain_R(Vhigh,1.9195,2.17,	0.09);
	AddToChain_R(Vhigh,1.9395,2.19,	0.11);
	AddToChain_R(Vhigh,1.961,	2.54,	0.26);
	AddToChain_R(Vhigh,1.986,	2.44,	0.12);
	AddToChain_R(Vhigh,2.05,	2.18,	0.07);
	AddToChain_R(Vhigh,2.2,	2.38,	0.07);
	AddToChain_R(Vhigh,2.375,	2.38,	0.07);
	AddToChain_R(Vhigh,2.495,	2.39,	0.08);
	AddToChain_R(Vhigh,2.595,	2.53,	0.04);
	AddToChain_R(Vhigh,2.7,	2.3,	0.07);
	AddToChain_R(Vhigh,2.8,	2.17,	0.06);
	AddToChain_R(Vhigh,2.9,	2.22,	0.07);
	AddToChain_R(Vhigh,3.025,	2.21,	0.05);
	AddToChain_R(Vhigh,2.62,	2.85,	0.74);
	AddToChain_R(Vhigh,2.8,	2.54,	0.46);
	AddToChain_R(Vhigh,3,	2.59,	0.15);
	AddToChain_R(Vhigh,3.2,	2.21,	0.07);
	AddToChain_R(Vhigh,3.4,	2.38,	0.07);
	AddToChain_R(Vhigh,3.575,	2.23,	0.06);
	AddToChain_R(Vhigh,3.685,	2.23,	0.08);
	AddToChain_R(Vhigh,3.73,	2.1,	0.08);
	AddToChain_R(Vhigh,3.7475,2.47,	0.09);
	AddToChain_R(Vhigh,3.7585,2.77,	0.11);
	AddToChain_R(Vhigh,3.764,	3.29,	0.27);
	AddToChain_R(Vhigh,3.7675,3.8,	0.33);
	AddToChain_R(Vhigh,3.77,	3.55,	0.14);
	AddToChain_R(Vhigh,3.7725,3.12,	0.24);
	AddToChain_R(Vhigh,3.776,	3.26,	0.26);
	AddToChain_R(Vhigh,3.7815,3.28,	0.12);
	AddToChain_R(Vhigh,3.791,	2.62,	0.11);
	AddToChain_R(Vhigh,3.8135,2.38,	0.1);
	AddToChain_R(Vhigh,3.85,	2.47,	0.11);
	AddToChain_R(Vhigh,3.895,	2.64,	0.11);
	AddToChain_R(Vhigh,3.9275,3.18,	0.14);
	AddToChain_R(Vhigh,3.94,	2.94,	0.13);
	AddToChain_R(Vhigh,3.95,	2.97,	0.13);
	AddToChain_R(Vhigh,3.96,	2.79,	0.12);
	AddToChain_R(Vhigh,3.97,	3.29,	0.13);
	AddToChain_R(Vhigh,3.98,	3.13,	0.14);
	AddToChain_R(Vhigh,3.9925,3.06,	0.15);
	AddToChain_R(Vhigh,4.0025,3.16,	0.14);
	AddToChain_R(Vhigh,4.01,	3.53,	0.16);
	AddToChain_R(Vhigh,4.02,	4.43,	0.16);
	AddToChain_R(Vhigh,4.0265,4.58,	0.18);
	AddToChain_R(Vhigh,4.03,	4.58,	0.2);
	AddToChain_R(Vhigh,4.034,	4.32,	0.17);
	AddToChain_R(Vhigh,4.0405,4.4,	0.17);
	AddToChain_R(Vhigh,4.05,	4.23,	0.17);
	AddToChain_R(Vhigh,4.06,	4.65,	0.19);
	AddToChain_R(Vhigh,4.07,	4.14,	0.2);
	AddToChain_R(Vhigh,4.08,	4.24,	0.21);
	AddToChain_R(Vhigh,4.09,	4.06,	0.17);
	AddToChain_R(Vhigh,4.1,	3.97,	0.16);
	AddToChain_R(Vhigh,4.11,	3.92,	0.16);
	AddToChain_R(Vhigh,4.12,	4.11,	0.24);
	AddToChain_R(Vhigh,4.13,	3.99,	0.15);
	AddToChain_R(Vhigh,4.14,	3.83,	0.15);
	AddToChain_R(Vhigh,4.15,	4.21,	0.18);
	AddToChain_R(Vhigh,4.16,	4.12,	0.15);
	AddToChain_R(Vhigh,4.17,	4.12,	0.15);
	AddToChain_R(Vhigh,4.18,	4.18,	0.17);
	AddToChain_R(Vhigh,4.19,	4.01,	0.14);
	AddToChain_R(Vhigh,4.2,	3.87,	0.16);
	AddToChain_R(Vhigh,4.21,	3.2,	0.16);
	AddToChain_R(Vhigh,4.22,	3.62,	0.15);
	AddToChain_R(Vhigh,4.23,	3.21,	0.13);
	AddToChain_R(Vhigh,4.2395,3.24,	0.12);
	AddToChain_R(Vhigh,4.246,	2.97,	0.11);
	AddToChain_R(Vhigh,4.25,	2.71,	0.12);
	AddToChain_R(Vhigh,4.255,	2.88,	0.11);
	AddToChain_R(Vhigh,4.26,	2.97,	0.11);
	AddToChain_R(Vhigh,4.265,	3.04,	0.13);
	AddToChain_R(Vhigh,4.27,	3.26,	0.12);
	AddToChain_R(Vhigh,4.281,	3.08,	0.12);
	AddToChain_R(Vhigh,4.3,	3.11,	0.12);
	AddToChain_R(Vhigh,4.32,	2.96,	0.12);
	AddToChain_R(Vhigh,4.3375,3.27,	0.15);
	AddToChain_R(Vhigh,4.35,	3.49,	0.14);
	AddToChain_R(Vhigh,4.3625,3.47,	0.13);
	AddToChain_R(Vhigh,4.3775,3.5,	0.15);
	AddToChain_R(Vhigh,4.39,	3.48,	0.16);
	AddToChain_R(Vhigh,4.4,	3.91,	0.16);
	AddToChain_R(Vhigh,4.41,	3.79,	0.15);
	AddToChain_R(Vhigh,4.42,	3.68,	0.14);
	AddToChain_R(Vhigh,4.43,	4.02,	0.16);
	AddToChain_R(Vhigh,4.44,	3.85,	0.17);
	AddToChain_R(Vhigh,4.45,	3.75,	0.15);
	AddToChain_R(Vhigh,4.46,	3.66,	0.17);
	AddToChain_R(Vhigh,4.475,	3.54,	0.17);
	AddToChain_R(Vhigh,4.4975,3.49,	0.14);
	AddToChain_R(Vhigh,4.52,	3.25,	0.13);
	AddToChain_R(Vhigh,4.54,	3.23,	0.14);
	AddToChain_R(Vhigh,4.565,	3.62,	0.13);
	AddToChain_R(Vhigh,4.64,	3.37,	0.1);
	AddToChain_R(Vhigh,4.775,	3.66,	0.14);
	AddToChain_R(Vhigh,4.925,	3.47,	0.32);
	AddToChain_R(Vhigh,5.05,	3.42,	0.12);
	AddToChain_R(Vhigh,5.25,	3.57,	0.11);
	AddToChain_R(Vhigh,5.5,	3.41,	0.1);
	AddToChain_R(Vhigh,5.75,	3.44,	0.11);
	AddToChain_R(Vhigh,6,	3.5,	0.1);
	AddToChain_R(Vhigh,6.25,	3.31,	0.1);
	AddToChain_R(Vhigh,6.5,	3.37,	0.1);
	AddToChain_R(Vhigh,6.75,	3.42,	0.09);
	AddToChain_R(Vhigh,7,	3.35,	0.1);
	AddToChain_R(Vhigh,7.15,	3.57,	0.11);
	AddToChain_R(Vhigh,7.3,	3.35,	0.14);
	AddToChain_R(Vhigh,7.3,	3.51,	0.11);
	AddToChain_R(Vhigh,7.425,	3.86,	0.19);
	AddToChain_R(Vhigh,7.465,	3.89,	0.19);
	AddToChain_R(Vhigh,7.515,	3.89,	0.19);
	AddToChain_R(Vhigh,7.6,	3.57,	0.19);
	AddToChain_R(Vhigh,7.7,	3.94,	0.2);
	AddToChain_R(Vhigh,7.8,	3.59,	0.2);
	AddToChain_R(Vhigh,7.875,	3.66,	0.2);
	AddToChain_R(Vhigh,7.975,	3.67,	0.2);
	AddToChain_R(Vhigh,8.1,	3.73,	0.17);
	AddToChain_R(Vhigh,8.2,	3.64,	0.17);
	AddToChain_R(Vhigh,8.3,	3.4,	0.17);
	AddToChain_R(Vhigh,8.4,	3.63,	0.17);
	AddToChain_R(Vhigh,8.5,	3.54,	0.18);
	AddToChain_R(Vhigh,8.608,	3.66,	0.18);
	AddToChain_R(Vhigh,8.708,	3.44,	0.09);
	AddToChain_R(Vhigh,8.8,	3.53,	0.09);
	AddToChain_R(Vhigh,8.9,	3.71,	0.09);
	AddToChain_R(Vhigh,9,	3.67,	0.12);
	AddToChain_R(Vhigh,9.1,	3.63,	0.12);
	AddToChain_R(Vhigh,9.2,	3.57,	0.12);
	AddToChain_R(Vhigh,10.295,3.48, 0.24);
	AddToChain_R(Vhigh,9.3,	3.6,	0.12);
	AddToChain_R(Vhigh,9.395,	3.41,	0.09);
	AddToChain_R(Vhigh,9.495,	3.63,	0.09);
	AddToChain_R(Vhigh,9.6,	3.55,	0.1);
	AddToChain_R(Vhigh,9.7,	3.52,	0.1);
	AddToChain_R(Vhigh,9.8,	3.58,	0.1);
	AddToChain_R(Vhigh,9.9,	3.49,	0.11);
	AddToChain_R(Vhigh,10,	3.68,	0.14);
	AddToChain_R(Vhigh,10.1,	3.56,	0.14);
	AddToChain_R(Vhigh,10.2,	3.49,	0.14);
	AddToChain_R(Vhigh,10.295,3.48,	0.24);

//--------- Взято из Behrend et al., Phys. Lett. 183B, 400(1987) 

//------------------------------------------------------------ CELLO
	AddToChain_R_percent(Vhigh,	14.04,	4.10,	2.6);
	AddToChain_R_percent(Vhigh,	22.00,	3.86,	3.0);
	AddToChain_R_percent(Vhigh,	33.80,	3.74,	2.6);
	AddToChain_R_percent(Vhigh,	38.28,	3.89,	2.6);
	AddToChain_R_percent(Vhigh,	41.50,	4.03,	4.1);
	AddToChain_R_percent(Vhigh,	44.20,	4.01,	2.5);
	AddToChain_R_percent(Vhigh,	46.60,	4.20,	8.5); 
//------------------------------------------------------------ JADE
	AddToChain_R_percent(Vhigh,	14.04,	3.94,	3.6);
	AddToChain_R_percent(Vhigh,	22.00,	4.11,	3.2);
	AddToChain_R_percent(Vhigh,	25.01,	4.24,	6.8);
	AddToChain_R_percent(Vhigh,	27.66,	3.85,	12.5);
	AddToChain_R_percent(Vhigh,	30.38,	3.85,	4.9);
	AddToChain_R_percent(Vhigh,	35.01,	3.94,	2.5);
	AddToChain_R_percent(Vhigh,	40.32,	4.07,	4.7);
	AddToChain_R_percent(Vhigh,	43.53,	4.05,	5.0);
	AddToChain_R_percent(Vhigh,	46.47,	4.11,	5.9); 
//------------------------------------------------------------ MARKJ
	AddToChain_R_percent(Vhigh,	22.00,	3.66,	2.2);
	AddToChain_R_percent(Vhigh,	25.0,	3.89,	5.4);
	AddToChain_R_percent(Vhigh,	30.60,	4.09,	3.4);
	AddToChain_R_percent(Vhigh,	33.82,	3.71,	1.6);
	AddToChain_R_percent(Vhigh,	36.36,	3.78,	4.0);
	AddToChain_R_percent(Vhigh,	40.36,	3.75,	4.0);
	AddToChain_R_percent(Vhigh,	43.58,	3.91,	1.5);
	AddToChain_R_percent(Vhigh,	45.48,	4.17,	4.8); 
//------------------------------------------------------------ PLUTO
	AddToChain_R_percent(Vhigh,	27.60,	4.07,	7.1);
	AddToChain_R_percent(Vhigh,	30.80,	4.11,	3.2); 
//------------------------------------------------------------ TASSO
	AddToChain_R_percent(Vhigh,	14.00,	4.14,	7.3);
	AddToChain_R_percent(Vhigh,	22.00,	3.89,	4.4);
	AddToChain_R_percent(Vhigh,	25.00,	3.72,	10.2);
	AddToChain_R_percent(Vhigh,	33.00,	3.74,	7.2);
	AddToChain_R_percent(Vhigh,	34.00,	4.14,	3.1);
	AddToChain_R_percent(Vhigh,	41.50,	4.11,	2.9);
	AddToChain_R_percent(Vhigh,	44.20,	4.28,	3.8); 
//------------------------------------------------------------

//----------------------------------- Разница между Брейт-Вигнером и истинным значением сечения для \phi

        Vphi=CreatePoint2_MeVnb(1010.,		-56.0913,	0.);

        AddToChain_MeVnb(Vphi,	1011,	-62.1797,	0.);
        AddToChain_MeVnb(Vphi,	1012,	-68.8979,	0.);
        AddToChain_MeVnb(Vphi,	1013,	-75.8778,	0.);
        AddToChain_MeVnb(Vphi,	1014,	-81.9922,	0.);
        AddToChain_MeVnb(Vphi,	1015,	-84.0626,	0.);
        AddToChain_MeVnb(Vphi,	1015.05,	-83.9423,	0.);
        AddToChain_MeVnb(Vphi,	1015.1,	-83.7907,	0.);
        AddToChain_MeVnb(Vphi,	1015.15,	-83.606,	0.);
        AddToChain_MeVnb(Vphi,	1015.2,	-83.388,	0.);
        AddToChain_MeVnb(Vphi,	1015.25,	-83.1713,	0.);
        AddToChain_MeVnb(Vphi,	1015.3,	-82.8802,	0.);
        AddToChain_MeVnb(Vphi,	1015.35,	-82.5498,	0.);
        AddToChain_MeVnb(Vphi,	1015.4,	-82.1763,	0.);
        AddToChain_MeVnb(Vphi,	1015.45,	-81.8014,	0.);
        AddToChain_MeVnb(Vphi,	1015.5,	-81.3408,	0.);
        AddToChain_MeVnb(Vphi,	1015.55,	-80.8317,	0.);
        AddToChain_MeVnb(Vphi,	1015.6,	-80.2724,	0.);
        AddToChain_MeVnb(Vphi,	1015.65,	-79.7043,	0.);
        AddToChain_MeVnb(Vphi,	1015.7,	-79.0848,	0.);
        AddToChain_MeVnb(Vphi,	1015.75,	-78.3618,	0.);
        AddToChain_MeVnb(Vphi,	1015.8,	-77.5785,	0.);
        AddToChain_MeVnb(Vphi,	1015.85,	-76.7828,	0.);
        AddToChain_MeVnb(Vphi,	1015.9,	-75.8695,	0.);
        AddToChain_MeVnb(Vphi,	1015.95,	-74.9432,	0.);
        AddToChain_MeVnb(Vphi,	1016,	-73.8931,	0.);
        AddToChain_MeVnb(Vphi,	1016.05,	-72.7966,	0.);
        AddToChain_MeVnb(Vphi,	1016.1,	-71.6531,	0.);
        AddToChain_MeVnb(Vphi,	1016.15,	-70.4294,	0.);
        AddToChain_MeVnb(Vphi,	1016.2,	-69.097,	0.);
        AddToChain_MeVnb(Vphi,	1016.25,	-67.6481,	0.);
        AddToChain_MeVnb(Vphi,	1016.3,	-66.2068,	0.);
        AddToChain_MeVnb(Vphi,	1016.35,	-64.5774,	0.);
        AddToChain_MeVnb(Vphi,	1016.4,	-62.9174,	0.);
        AddToChain_MeVnb(Vphi,	1016.45,	-61.13,	0.);
        AddToChain_MeVnb(Vphi,	1016.5,	-59.243,	0.);
        AddToChain_MeVnb(Vphi,	1016.55,	-57.3278,	0.);
        AddToChain_MeVnb(Vphi,	1016.6,	-55.2351,	0.);
        AddToChain_MeVnb(Vphi,	1016.65,	-53.0698,	0.);
        AddToChain_MeVnb(Vphi,	1016.7,	-50.7611,	0.);
        AddToChain_MeVnb(Vphi,	1016.75,	-48.4242,	0.);
        AddToChain_MeVnb(Vphi,	1016.8,	-45.8931,	0.);
        AddToChain_MeVnb(Vphi,	1016.85,	-43.2488,	0.);
        AddToChain_MeVnb(Vphi,	1016.9,	-40.5297,	0.);
        AddToChain_MeVnb(Vphi,	1016.95,	-37.6575,	0.);
        AddToChain_MeVnb(Vphi,	1017,	-34.766,	0.);
        AddToChain_MeVnb(Vphi,	1017.05,	-31.7185,	0.);
        AddToChain_MeVnb(Vphi,	1017.1,	-28.5131,	0.);
        AddToChain_MeVnb(Vphi,	1017.15,	-25.2456,	0.);
        AddToChain_MeVnb(Vphi,	1017.2,	-21.9886,	0.);
        AddToChain_MeVnb(Vphi,	1017.25,	-18.4795,	0.);
        AddToChain_MeVnb(Vphi,	1017.3,	-14.934,	0.);
        AddToChain_MeVnb(Vphi,	1017.35,	-11.364,	0.);
        AddToChain_MeVnb(Vphi,	1017.4,	-7.60278,	0.);
        AddToChain_MeVnb(Vphi,	1017.45,	-3.90259,	0.);
        AddToChain_MeVnb(Vphi,	1017.5,	-0.0969238,	0.);
        AddToChain_MeVnb(Vphi,	1017.55,	3.7417,	0.);
        AddToChain_MeVnb(Vphi,	1017.6,	7.53418,	0.);
        AddToChain_MeVnb(Vphi,	1017.65,	11.3293,	0.);
        AddToChain_MeVnb(Vphi,	1017.7,	15.156,	0.);
        AddToChain_MeVnb(Vphi,	1017.75,	18.9365,	0.);
        AddToChain_MeVnb(Vphi,	1017.8,	22.5085,	0.);
        AddToChain_MeVnb(Vphi,	1017.85,	26.1125,	0.);
        AddToChain_MeVnb(Vphi,	1017.9,	29.5222,	0.);
        AddToChain_MeVnb(Vphi,	1017.95,	32.8279,	0.);
        AddToChain_MeVnb(Vphi,	1018,	35.9333,	0.);
        AddToChain_MeVnb(Vphi,	1018.05,	38.657,	0.);
        AddToChain_MeVnb(Vphi,	1018.1,	41.2444,	0.);
        AddToChain_MeVnb(Vphi,	1018.15,	43.4451,	0.);
        AddToChain_MeVnb(Vphi,	1018.2,	45.2773,	0.);
        AddToChain_MeVnb(Vphi,	1018.25,	46.6301,	0.);
        AddToChain_MeVnb(Vphi,	1018.3,	47.6794,	0.);
        AddToChain_MeVnb(Vphi,	1018.35,	48.1616,	0.);
        AddToChain_MeVnb(Vphi,	1018.4,	48.1101,	0.);
        AddToChain_MeVnb(Vphi,	1018.45,	47.478,	0.);
        AddToChain_MeVnb(Vphi,	1018.5,	46.1597,	0.);
        AddToChain_MeVnb(Vphi,	1018.55,	44.2625,	0.);
        AddToChain_MeVnb(Vphi,	1018.6,	41.7505,	0.);
        AddToChain_MeVnb(Vphi,	1018.65,	38.3953,	0.);
        AddToChain_MeVnb(Vphi,	1018.7,	34.4456,	0.);
        AddToChain_MeVnb(Vphi,	1018.75,	29.5684,	0.);
        AddToChain_MeVnb(Vphi,	1018.8,	24.075,	0.);
        AddToChain_MeVnb(Vphi,	1018.85,	17.8435,	0.);
        AddToChain_MeVnb(Vphi,	1018.9,	10.8279,	0.);
        AddToChain_MeVnb(Vphi,	1018.95,	3.1228,	0.);
        AddToChain_MeVnb(Vphi,	1019,	-5.19702,	0.);
        AddToChain_MeVnb(Vphi,	1019.05,	-14.1399,	0.);
        AddToChain_MeVnb(Vphi,	1019.1,	-23.7292,	0.);
        AddToChain_MeVnb(Vphi,	1019.15,	-33.7192,	0.);
        AddToChain_MeVnb(Vphi,	1019.2,	-44.2312,	0.);
        AddToChain_MeVnb(Vphi,	1019.25,	-55.0583,	0.);
        AddToChain_MeVnb(Vphi,	1019.3,	-66.1265,	0.);
        AddToChain_MeVnb(Vphi,	1019.35,	-77.4021,	0.);
        AddToChain_MeVnb(Vphi,	1019.4,	-88.7495,	0.);
        AddToChain_MeVnb(Vphi,	1019.45,	-100.06,	0.);
        AddToChain_MeVnb(Vphi,	1019.5,	-111.26,	0.);
        AddToChain_MeVnb(Vphi,	1019.55,	-122.24,	0.);
        AddToChain_MeVnb(Vphi,	1019.6,	-132.917,	0.);
        AddToChain_MeVnb(Vphi,	1019.65,	-143.213,	0.);
        AddToChain_MeVnb(Vphi,	1019.7,	-153.022,	0.);
        AddToChain_MeVnb(Vphi,	1019.75,	-162.31,	0.);
        AddToChain_MeVnb(Vphi,	1019.8,	-170.897,	0.);
        AddToChain_MeVnb(Vphi,	1019.85,	-178.906,	0.);
        AddToChain_MeVnb(Vphi,	1019.9,	-186.135,	0.);
        AddToChain_MeVnb(Vphi,	1019.95,	-192.634,	0.);
        AddToChain_MeVnb(Vphi,	1020,	-198.305,	0.);
        AddToChain_MeVnb(Vphi,	1020.05,	-203.222,	0.);
        AddToChain_MeVnb(Vphi,	1020.1,	-207.335,	0.);
        AddToChain_MeVnb(Vphi,	1020.15,	-210.591,	0.);
        AddToChain_MeVnb(Vphi,	1020.2,	-213.102,	0.);
        AddToChain_MeVnb(Vphi,	1020.25,	-214.711,	0.);
        AddToChain_MeVnb(Vphi,	1020.3,	-215.687,	0.);
        AddToChain_MeVnb(Vphi,	1020.35,	-215.942,	0.);
        AddToChain_MeVnb(Vphi,	1020.4,	-215.442,	0.);
        AddToChain_MeVnb(Vphi,	1020.45,	-214.216,	0.);
        AddToChain_MeVnb(Vphi,	1020.5,	-212.518,	0.);
        AddToChain_MeVnb(Vphi,	1020.55,	-210.184,	0.);
        AddToChain_MeVnb(Vphi,	1020.6,	-207.255,	0.);
        AddToChain_MeVnb(Vphi,	1020.65,	-204.067,	0.);
        AddToChain_MeVnb(Vphi,	1020.7,	-200.232,	0.);
        AddToChain_MeVnb(Vphi,	1020.75,	-196.086,	0.);
        AddToChain_MeVnb(Vphi,	1020.8,	-191.672,	0.);
        AddToChain_MeVnb(Vphi,	1020.85,	-186.744,	0.);
        AddToChain_MeVnb(Vphi,	1020.9,	-181.702,	0.);
        AddToChain_MeVnb(Vphi,	1020.95,	-176.446,	0.);
        AddToChain_MeVnb(Vphi,	1021,	-171.081,	0.);
        AddToChain_MeVnb(Vphi,	1021.05,	-165.361,	0.);
        AddToChain_MeVnb(Vphi,	1021.1,	-159.669,	0.);
        AddToChain_MeVnb(Vphi,	1021.15,	-153.89,	0.);
        AddToChain_MeVnb(Vphi,	1021.2,	-148.06,	0.);
        AddToChain_MeVnb(Vphi,	1021.25,	-142.2,	0.);
        AddToChain_MeVnb(Vphi,	1021.3,	-136.205,	0.);
        AddToChain_MeVnb(Vphi,	1021.35,	-130.353,	0.);
        AddToChain_MeVnb(Vphi,	1021.4,	-124.468,	0.);
        AddToChain_MeVnb(Vphi,	1021.45,	-118.698,	0.);
        AddToChain_MeVnb(Vphi,	1021.5,	-112.99,	0.);
        AddToChain_MeVnb(Vphi,	1021.55,	-107.241,	0.);
        AddToChain_MeVnb(Vphi,	1021.6,	-101.699,	0.);
        AddToChain_MeVnb(Vphi,	1021.65,	-96.1914,	0.);
        AddToChain_MeVnb(Vphi,	1021.7,	-90.8492,	0.);
        AddToChain_MeVnb(Vphi,	1021.75,	-85.5114,	0.);
        AddToChain_MeVnb(Vphi,	1021.8,	-80.3466,	0.);
        AddToChain_MeVnb(Vphi,	1021.85,	-75.3549,	0.);
        AddToChain_MeVnb(Vphi,	1021.9,	-70.4344,	0.);
        AddToChain_MeVnb(Vphi,	1021.95,	-65.6469,	0.);
        AddToChain_MeVnb(Vphi,	1022,	-60.9435,	0.);
        AddToChain_MeVnb(Vphi,	1022.05,	-56.4193,	0.);
        AddToChain_MeVnb(Vphi,	1022.1,	-52.0715,	0.);
        AddToChain_MeVnb(Vphi,	1022.15,	-47.7644,	0.);
        AddToChain_MeVnb(Vphi,	1022.2,	-43.5941,	0.);
        AddToChain_MeVnb(Vphi,	1022.25,	-39.5973,	0.);
        AddToChain_MeVnb(Vphi,	1022.3,	-35.729,	0.);
        AddToChain_MeVnb(Vphi,	1022.35,	-31.9495,	0.);
        AddToChain_MeVnb(Vphi,	1022.4,	-28.3693,	0.);
        AddToChain_MeVnb(Vphi,	1022.45,	-24.8033,	0.);
        AddToChain_MeVnb(Vphi,	1022.5,	-21.4329,	0.);
        AddToChain_MeVnb(Vphi,	1022.55,	-18.1801,	0.);
        AddToChain_MeVnb(Vphi,	1022.6,	-14.9753,	0.);
        AddToChain_MeVnb(Vphi,	1022.65,	-11.9171,	0.);
        AddToChain_MeVnb(Vphi,	1022.7,	-9.00415,	0.);
        AddToChain_MeVnb(Vphi,	1022.75,	-6.19739,	0.);
        AddToChain_MeVnb(Vphi,	1022.8,	-3.43457,	0.);
        AddToChain_MeVnb(Vphi,	1022.85,	-0.833618,	0.);
        AddToChain_MeVnb(Vphi,	1022.9,	1.69995,	0.);
        AddToChain_MeVnb(Vphi,	1022.95,	4.10767,	0.);
        AddToChain_MeVnb(Vphi,	1023,	6.45068,	0.);
        AddToChain_MeVnb(Vphi,	1023.05,	8.72974,	0.);
        AddToChain_MeVnb(Vphi,	1023.1,	10.869,	0.);
        AddToChain_MeVnb(Vphi,	1023.15,	12.9762,	0.);
        AddToChain_MeVnb(Vphi,	1023.2,	14.9509,	0.);
        AddToChain_MeVnb(Vphi,	1023.25,	16.9185,	0.);
        AddToChain_MeVnb(Vphi,	1023.3,	18.7389,	0.);
        AddToChain_MeVnb(Vphi,	1023.35,	20.5093,	0.);
        AddToChain_MeVnb(Vphi,	1023.4,	22.2097,	0.);
        AddToChain_MeVnb(Vphi,	1023.45,	23.8408,	0.);
        AddToChain_MeVnb(Vphi,	1023.5,	25.4471,	0.);
        AddToChain_MeVnb(Vphi,	1023.55,	26.9276,	0.);
        AddToChain_MeVnb(Vphi,	1023.6,	28.3672,	0.);
        AddToChain_MeVnb(Vphi,	1023.65,	29.7679,	0.);
        AddToChain_MeVnb(Vphi,	1023.7,	31.1282,	0.);
        AddToChain_MeVnb(Vphi,	1023.75,	32.3959,	0.);
        AddToChain_MeVnb(Vphi,	1023.8,	33.6107,	0.);
        AddToChain_MeVnb(Vphi,	1023.85,	34.8092,	0.);
        AddToChain_MeVnb(Vphi,	1023.9,	35.9075,	0.);
        AddToChain_MeVnb(Vphi,	1023.95,	37.0077,	0.);
        AddToChain_MeVnb(Vphi,	1024,	38.0296,	0.);
        AddToChain_MeVnb(Vphi,	1024.05,	39.0236,	0.);
        AddToChain_MeVnb(Vphi,	1024.1,	39.9742,	0.);
        AddToChain_MeVnb(Vphi,	1024.15,	40.855,	0.);
        AddToChain_MeVnb(Vphi,	1024.2,	41.7401,	0.);
        AddToChain_MeVnb(Vphi,	1024.25,	42.5733,	0.);
        AddToChain_MeVnb(Vphi,	1024.3,	43.3696,	0.);
        AddToChain_MeVnb(Vphi,	1024.35,	44.1435,	0.);
        AddToChain_MeVnb(Vphi,	1024.4,	44.8705,	0.);
        AddToChain_MeVnb(Vphi,	1024.45,	45.5518,	0.);
        AddToChain_MeVnb(Vphi,	1024.5,	46.2151,	0.);
        AddToChain_MeVnb(Vphi,	1024.55,	46.8605,	0.);
        AddToChain_MeVnb(Vphi,	1024.6,	47.4641,	0.);
        AddToChain_MeVnb(Vphi,	1024.65,	48.0516,	0.);
        AddToChain_MeVnb(Vphi,	1024.7,	48.6117,	0.);
        AddToChain_MeVnb(Vphi,	1024.75,	49.1235,	0.);
        AddToChain_MeVnb(Vphi,	1024.8,	49.6218,	0.);
        AddToChain_MeVnb(Vphi,	1024.85,	50.1065,	0.);
        AddToChain_MeVnb(Vphi,	1024.9,	50.568,	0.);
        AddToChain_MeVnb(Vphi,	1024.95,	51.007,	0.);
        AddToChain_MeVnb(Vphi,	1025,	51.4242,	0.);
        AddToChain_MeVnb(Vphi,	1025.05,	51.8013,	0.);
        AddToChain_MeVnb(Vphi,	1025.1,	52.1782,	0.);
        AddToChain_MeVnb(Vphi,	1025.15,	52.5357,	0.);
        AddToChain_MeVnb(Vphi,	1025.2,	52.8747,	0.);
        AddToChain_MeVnb(Vphi,	1025.25,	53.1957,	0.);
        AddToChain_MeVnb(Vphi,	1025.3,	53.5081,	0.);
        AddToChain_MeVnb(Vphi,	1025.35,	53.7786,	0.);
        AddToChain_MeVnb(Vphi,	1025.4,	54.0503,	0.);
        AddToChain_MeVnb(Vphi,	1025.45,	54.3067,	0.);
        AddToChain_MeVnb(Vphi,	1025.5,	54.5562,	0.);
        AddToChain_MeVnb(Vphi,	1025.55,	54.7834,	0.);
        AddToChain_MeVnb(Vphi,	1025.6,	55.0049,	0.);
        AddToChain_MeVnb(Vphi,	1025.65,	55.1905,	0.);
        AddToChain_MeVnb(Vphi,	1025.7,	55.3861,	0.);
        AddToChain_MeVnb(Vphi,	1025.75,	55.5619,	0.);
        AddToChain_MeVnb(Vphi,	1025.8,	55.7329,	0.);
        AddToChain_MeVnb(Vphi,	1025.85,	55.8859,	0.);
        AddToChain_MeVnb(Vphi,	1025.9,	56.0346,	0.);
        AddToChain_MeVnb(Vphi,	1025.95,	56.1598,	0.);
        AddToChain_MeVnb(Vphi,	1026,	56.2879,	0.);
        AddToChain_MeVnb(Vphi,	1026.05,	56.4061,	0.);
        AddToChain_MeVnb(Vphi,	1026.1,	56.5089,	0.);
        AddToChain_MeVnb(Vphi,	1026.15,	56.6085,	0.);
        AddToChain_MeVnb(Vphi,	1026.2,	56.6997,	0.);
        AddToChain_MeVnb(Vphi,	1026.25,	56.7705,	0.);
        AddToChain_MeVnb(Vphi,	1026.3,	56.845,	0.);
        AddToChain_MeVnb(Vphi,	1026.35,	56.9121,	0.);
        AddToChain_MeVnb(Vphi,	1026.4,	56.9766,	0.);
        AddToChain_MeVnb(Vphi,	1026.45,	57.0289,	0.);
        AddToChain_MeVnb(Vphi,	1026.5,	57.0736,	0.);
        AddToChain_MeVnb(Vphi,	1026.55,	57.1014,	0.);
        AddToChain_MeVnb(Vphi,	1026.6,	57.1337,	0.);
        AddToChain_MeVnb(Vphi,	1026.65,	57.1644,	0.);
        AddToChain_MeVnb(Vphi,	1026.7,	57.1844,	0.);
        AddToChain_MeVnb(Vphi,	1026.75,	57.2031,	0.);
        AddToChain_MeVnb(Vphi,	1026.8,	57.2113,	0.);
        AddToChain_MeVnb(Vphi,	1026.85,	57.2101,	0.);
        AddToChain_MeVnb(Vphi,	1026.9,	57.208,	0.);
        AddToChain_MeVnb(Vphi,	1026.95,	57.206,	0.);
        AddToChain_MeVnb(Vphi,	1027,	57.1942,	0.);
        AddToChain_MeVnb(Vphi,	1027.05,	57.1822,	0.);
        AddToChain_MeVnb(Vphi,	1027.1,	57.1577,	0.);
        AddToChain_MeVnb(Vphi,	1027.15,	57.1328,	0.);
        AddToChain_MeVnb(Vphi,	1027.2,	57.1086,	0.);
        AddToChain_MeVnb(Vphi,	1027.25,	57.0798,	0.);
        AddToChain_MeVnb(Vphi,	1027.3,	57.0471,	0.);
        AddToChain_MeVnb(Vphi,	1027.35,	57.0113,	0.);
        AddToChain_MeVnb(Vphi,	1027.4,	56.9637,	0.);
        AddToChain_MeVnb(Vphi,	1027.45,	56.9211,	0.);
        AddToChain_MeVnb(Vphi,	1027.5,	56.8746,	0.);
        AddToChain_MeVnb(Vphi,	1027.55,	56.825,	0.);
        AddToChain_MeVnb(Vphi,	1027.6,	56.7728,	0.);
        AddToChain_MeVnb(Vphi,	1027.65,	56.7173,	0.);
        AddToChain_MeVnb(Vphi,	1027.7,	56.6526,	0.);
        AddToChain_MeVnb(Vphi,	1027.75,	56.595,	0.);
        AddToChain_MeVnb(Vphi,	1027.8,	56.5315,	0.);
        AddToChain_MeVnb(Vphi,	1027.85,	56.466,	0.);
        AddToChain_MeVnb(Vphi,	1027.9,	56.4008,	0.);
        AddToChain_MeVnb(Vphi,	1027.95,	56.3307,	0.);
        AddToChain_MeVnb(Vphi,	1028,	56.2547,	0.);
        AddToChain_MeVnb(Vphi,	1028.05,	56.1799,	0.);
        AddToChain_MeVnb(Vphi,	1028.1,	56.1065,	0.);
        AddToChain_MeVnb(Vphi,	1028.15,	56.0276,	0.);
        AddToChain_MeVnb(Vphi,	1028.2,	55.9504,	0.);
        AddToChain_MeVnb(Vphi,	1028.25,	55.8708,	0.);
        AddToChain_MeVnb(Vphi,	1028.3,	55.781,	0.);
        AddToChain_MeVnb(Vphi,	1028.35,	55.6986,	0.);
        AddToChain_MeVnb(Vphi,	1028.4,	55.6141,	0.);
        AddToChain_MeVnb(Vphi,	1028.45,	55.5286,	0.);
        AddToChain_MeVnb(Vphi,	1028.5,	55.441,	0.);
        AddToChain_MeVnb(Vphi,	1028.55,	55.3521,	0.);
        AddToChain_MeVnb(Vphi,	1028.6,	55.2572,	0.);
        AddToChain_MeVnb(Vphi,	1028.65,	55.1656,	0.);
        AddToChain_MeVnb(Vphi,	1028.7,	55.0734,	0.);
        AddToChain_MeVnb(Vphi,	1028.75,	54.9794,	0.);
        AddToChain_MeVnb(Vphi,	1028.8,	54.8845,	0.);
        AddToChain_MeVnb(Vphi,	1028.85,	54.7889,	0.);
        AddToChain_MeVnb(Vphi,	1028.9,	54.6871,	0.);
        AddToChain_MeVnb(Vphi,	1028.95,	54.5919,	0.);
        AddToChain_MeVnb(Vphi,	1029,	54.493,	0.);
        AddToChain_MeVnb(Vphi,	1029.05,	54.3932,	0.);
        AddToChain_MeVnb(Vphi,	1029.1,	54.2953,	0.);
        AddToChain_MeVnb(Vphi,	1029.15,	54.1939,	0.);
        AddToChain_MeVnb(Vphi,	1029.2,	54.0878,	0.);
        AddToChain_MeVnb(Vphi,	1029.25,	53.9871,	0.);
        AddToChain_MeVnb(Vphi,	1029.3,	53.8858,	0.);
        AddToChain_MeVnb(Vphi,	1029.35,	53.7821,	0.);
        AddToChain_MeVnb(Vphi,	1029.4,	53.6794,	0.);
        AddToChain_MeVnb(Vphi,	1029.45,	53.5745,	0.);
        AddToChain_MeVnb(Vphi,	1029.5,	53.4667,	0.);
        AddToChain_MeVnb(Vphi,	1029.55,	53.3624,	0.);
        AddToChain_MeVnb(Vphi,	1029.6,	53.2579,	0.);
        AddToChain_MeVnb(Vphi,	1029.65,	53.1526,	0.);
        AddToChain_MeVnb(Vphi,	1029.7,	53.0474,	0.);
        AddToChain_MeVnb(Vphi,	1029.75,	52.9356,	0.);
        AddToChain_MeVnb(Vphi,	1029.8,	52.8292,	0.);
        AddToChain_MeVnb(Vphi,	1029.85,	52.7228,	0.);
        AddToChain_MeVnb(Vphi,	1029.9,	52.6174,	0.);
        AddToChain_MeVnb(Vphi,	1029.95,	52.5103,	0.);
        AddToChain_MeVnb(Vphi,	1030,	52.4024,	0.);
        AddToChain_MeVnb(Vphi,	1031,	50.2261,	0.);
        AddToChain_MeVnb(Vphi,	1032,	48.0944,	0.);
        AddToChain_MeVnb(Vphi,	1033,	46.0572,	0.);
        AddToChain_MeVnb(Vphi,	1034,	44.1436,	0.);
        AddToChain_MeVnb(Vphi,	1035,	42.3584,	0.);
        AddToChain_MeVnb(Vphi,	1036,	40.6985,	0.);
        AddToChain_MeVnb(Vphi,	1037,	39.159,	0.);
        AddToChain_MeVnb(Vphi,	1038,	37.731,	0.);
        AddToChain_MeVnb(Vphi,	1039,	36.4049,	0.);
        AddToChain_MeVnb(Vphi,	1040,	35.1729,	0.);
        AddToChain_MeVnb(Vphi,	1041,	34.0256,	0.);
        AddToChain_MeVnb(Vphi,	1042,	32.9564,	0.);
        AddToChain_MeVnb(Vphi,	1043,	31.9574,	0.);
        AddToChain_MeVnb(Vphi,	1044,	31.022,	0.);
        AddToChain_MeVnb(Vphi,	1045,	30.1453,	0.);
        AddToChain_MeVnb(Vphi,	1046,	29.3214,	0.);
        AddToChain_MeVnb(Vphi,	1047,	28.5462,	0.);
        AddToChain_MeVnb(Vphi,	1048,	27.8154,	0.);
        AddToChain_MeVnb(Vphi,	1049,	27.1251,	0.);
        AddToChain_MeVnb(Vphi,	1050,	26.4722,	0.);
        AddToChain_MeVnb(Vphi,	1051,	25.8534,	0.);
        AddToChain_MeVnb(Vphi,	1052,	25.2663,	0.);
        AddToChain_MeVnb(Vphi,	1053,	24.7084,	0.);
        AddToChain_MeVnb(Vphi,	1054,	24.1773,	0.);
        AddToChain_MeVnb(Vphi,	1055,	23.6714,	0.);
        AddToChain_MeVnb(Vphi,	1056,	23.1885,	0.);
        AddToChain_MeVnb(Vphi,	1057,	22.7272,	0.);
        AddToChain_MeVnb(Vphi,	1058,	22.286,	0.);
        AddToChain_MeVnb(Vphi,	1059,	21.8635,	0.);
        AddToChain_MeVnb(Vphi,	1060,	21.4584,	0.);
        AddToChain_MeVnb(Vphi,	1061,	21.0697,	0.);
        AddToChain_MeVnb(Vphi,	1062,	20.6962,	0.);
        AddToChain_MeVnb(Vphi,	1063,	20.3372,	0.);
        AddToChain_MeVnb(Vphi,	1064,	19.9915,	0.);
        AddToChain_MeVnb(Vphi,	1065,	19.6586,	0.);
        AddToChain_MeVnb(Vphi,	1066,	19.3376,	0.);
        AddToChain_MeVnb(Vphi,	1067,	19.0277,	0.);
        AddToChain_MeVnb(Vphi,	1068,	18.7286,	0.);
        AddToChain_MeVnb(Vphi,	1069,	18.4394,	0.);
        AddToChain_MeVnb(Vphi,	1070,	18.1597,	0.);
        AddToChain_MeVnb(Vphi,	1071,	17.889,	0.);
        AddToChain_MeVnb(Vphi,	1072,	17.6268,	0.);
        AddToChain_MeVnb(Vphi,	1073,	17.3727,	0.);
        AddToChain_MeVnb(Vphi,	1074,	17.1262,	0.);
        AddToChain_MeVnb(Vphi,	1075,	16.887,	0.);
        AddToChain_MeVnb(Vphi,	1076,	16.6547,	0.);
        AddToChain_MeVnb(Vphi,	1077,	16.429,	0.);
        AddToChain_MeVnb(Vphi,	1078,	16.2096,	0.);
        AddToChain_MeVnb(Vphi,	1079,	15.9962,	0.);
        AddToChain_MeVnb(Vphi,	1080,	15.7885,	0.);
        AddToChain_MeVnb(Vphi,	1081,	15.5863,	0.);
        AddToChain_MeVnb(Vphi,	1082,	15.3893,	0.);
        AddToChain_MeVnb(Vphi,	1083,	15.1973,	0.);
        AddToChain_MeVnb(Vphi,	1084,	15.0101,	0.);
        AddToChain_MeVnb(Vphi,	1085,	14.8275,	0.);
        AddToChain_MeVnb(Vphi,	1086,	14.6492,	0.);
        AddToChain_MeVnb(Vphi,	1087,	14.4752,	0.);
        AddToChain_MeVnb(Vphi,	1088,	14.3053,	0.);
        AddToChain_MeVnb(Vphi,	1089,	14.1392,	0.);
        AddToChain_MeVnb(Vphi,	1090,	13.9768,	0.);
        AddToChain_MeVnb(Vphi,	1091,	13.8181,	0.);
        AddToChain_MeVnb(Vphi,	1092,	13.6628,	0.);
        AddToChain_MeVnb(Vphi,	1093,	13.5108,	0.);
        AddToChain_MeVnb(Vphi,	1094,	13.3621,	0.);
        AddToChain_MeVnb(Vphi,	1095,	13.2164,	0.);
        AddToChain_MeVnb(Vphi,	1096,	13.0737,	0.);
        AddToChain_MeVnb(Vphi,	1097,	12.9339,	0.);
        AddToChain_MeVnb(Vphi,	1098,	12.7969,	0.);
        AddToChain_MeVnb(Vphi,	1099,	12.6626,	0.);
        AddToChain_MeVnb(Vphi,	1100,	12.5308,	0.);
        AddToChain_MeVnb(Vphi,	1101,	12.4019,	0.);
        AddToChain_MeVnb(Vphi,	1102,	12.2753,	0.);
        AddToChain_MeVnb(Vphi,	1103,	12.1511,	0.);
        AddToChain_MeVnb(Vphi,	1104,	12.0291,	0.);
        AddToChain_MeVnb(Vphi,	1105,	11.9094,	0.);
        AddToChain_MeVnb(Vphi,	1106,	11.7918,	0.);
        AddToChain_MeVnb(Vphi,	1107,	11.6763,	0.);
        AddToChain_MeVnb(Vphi,	1108,	11.5627,	0.);
        AddToChain_MeVnb(Vphi,	1109,	11.4512,	0.);
        AddToChain_MeVnb(Vphi,	1110,	11.3415,	0.);
/*
        AddToChain_MeVnb(Vphi,	1111,	11.2336,	0.);
        AddToChain_MeVnb(Vphi,	1112,	11.1276,	0.);
        AddToChain_MeVnb(Vphi,	1113,	11.0232,	0.);
        AddToChain_MeVnb(Vphi,	1114,	10.9206,	0.);
        AddToChain_MeVnb(Vphi,	1115,	10.8196,	0.);
        AddToChain_MeVnb(Vphi,	1116,	10.7201,	0.);
        AddToChain_MeVnb(Vphi,	1117,	10.6223,	0.);
        AddToChain_MeVnb(Vphi,	1118,	10.5259,	0.);
        AddToChain_MeVnb(Vphi,	1119,	10.431,	0.);
        AddToChain_MeVnb(Vphi,	1120,	10.3376,	0.);
        AddToChain_MeVnb(Vphi,	1121,	10.2455,	0.);
        AddToChain_MeVnb(Vphi,	1122,	10.1548,	0.);
        AddToChain_MeVnb(Vphi,	1123,	10.0654,	0.);
        AddToChain_MeVnb(Vphi,	1124,	9.97731,	0.);
        AddToChain_MeVnb(Vphi,	1125,	9.89047,	0.);
        AddToChain_MeVnb(Vphi,	1126,	9.80487,	0.);
        AddToChain_MeVnb(Vphi,	1127,	9.72047,	0.);
        AddToChain_MeVnb(Vphi,	1128,	9.63724,	0.);
        AddToChain_MeVnb(Vphi,	1129,	9.55516,	0.);
        AddToChain_MeVnb(Vphi,	1130,	9.47419,	0.);
*/
//------------------------------------------------------------------ K^+ K^- (Olya)
        Vkpkm=CreatePoint2_nb(1.110000,	8.200000,	6.000000);

        AddToChain_nb(Vkpkm,	1.130000,	10.00000,	2.600000);
        AddToChain_nb(Vkpkm,	1.150000,	8.600000,	1.300000);
        AddToChain_nb(Vkpkm,	1.170000,	8.000000,	1.100000);
        AddToChain_nb(Vkpkm,	1.190000,	8.100000,	0.9000000);
        AddToChain_nb(Vkpkm,	1.210000,	5.200000,	0.7000000);
        AddToChain_nb(Vkpkm,	1.230000,	6.400000,	0.8000000);
        AddToChain_nb(Vkpkm,	1.250000,	5.900000,	0.7000000);
        AddToChain_nb(Vkpkm,	1.270000,	6.900000,	0.7000000);
        AddToChain_nb(Vkpkm,	1.290000,	5.800000,	0.7000000);
        AddToChain_nb(Vkpkm,	1.310000,	6.600000,	1.100000);
        AddToChain_nb(Vkpkm,	1.330000,	5.800000,	0.8000000);
        AddToChain_nb(Vkpkm,	1.350000,	6.100000,	0.6000000);
        AddToChain_nb(Vkpkm,	1.370000,	6.000000,	0.6000000);
        AddToChain_nb(Vkpkm,	1.380000,	5.500000,	0.6000000);

//------------------------------------------------------------------ K_S K_L (Петя Лукин)

        Vkskl=CreatePoint2_MeVnb( 1110.000,	3.424731,	 0.6242384);
        
        AddToChain_MeVnb(Vkskl,	1120.000,	2.674277,	 0.6743487);
        AddToChain_MeVnb(Vkskl,	1130.000,	2.153516,	 0.5738519);
        AddToChain_MeVnb(Vkskl,	1140.000,	2.676342,	 0.5633617);
        AddToChain_MeVnb(Vkskl,	1150.000,	2.645177,	 0.6632734);
        AddToChain_MeVnb(Vkskl,	1160.000,	1.526746,	 0.3955862);
        AddToChain_MeVnb(Vkskl,	1170.000,	1.928362,	 0.5413375);
        AddToChain_MeVnb(Vkskl,	1180.000,	1.368014,	 0.3877933);
        AddToChain_MeVnb(Vkskl,	1190.000,	1.658339,	 0.3967198);
        AddToChain_MeVnb(Vkskl,	1204.600,	0.8258442,	 0.2198748);
        AddToChain_MeVnb(Vkskl,	1225.000,	1.083510,	 0.2320116);
        AddToChain_MeVnb(Vkskl,	1250.600,	0.5480630,	 0.1400206);
        AddToChain_MeVnb(Vkskl,	1275.000,	0.5322962,	 0.1564946);
        AddToChain_MeVnb(Vkskl,	1295.800,	0.4784661,	 0.1342008);
        AddToChain_MeVnb(Vkskl,	1325.300,	0.1834043,	 0.9547332E-01);
        AddToChain_MeVnb(Vkskl,	1368.300,	0.2478264,	 0.7589923E-01);
               
//------------------------------------------------------------------- \pi^+\pi^-\pi^o

//---------- ND  (S.I.Dolinsky et al., Phys. Rep. 202 (1991) 99)

        V3pi_low=CreatePoint2_MeVnb(  661.,   0.0,2.4);

        AddToChain_MeVnb(V3pi_low,  671,   0.0,2.0);
        AddToChain_MeVnb(V3pi_low,  681,   1.5,4.1);
        AddToChain_MeVnb(V3pi_low,  691,   1.4,4.1);
        AddToChain_MeVnb(V3pi_low,  701,   1.6,2.7);
        AddToChain_MeVnb(V3pi_low,  711,   1.1,3.6);
        AddToChain_MeVnb(V3pi_low,  724,   4.3,1.8);
        AddToChain_MeVnb(V3pi_low,  735,   6.7,2.0);
        AddToChain_MeVnb(V3pi_low,  745,   14.1,1.9);
        AddToChain_MeVnb(V3pi_low,  755,   23.4,3.2);
        AddToChain_MeVnb(V3pi_low,  765,   63.9,5.2);
        AddToChain_MeVnb(V3pi_low,  805,   74.9,3.7);
        AddToChain_MeVnb(V3pi_low,  815,   41.7,2.1);
        AddToChain_MeVnb(V3pi_low,  824,   28.2,1.7);
        AddToChain_MeVnb(V3pi_low,  832,   24.8,3.5);
        AddToChain_MeVnb(V3pi_low,  841,   18.9,2.8);
        AddToChain_MeVnb(V3pi_low,  851,   13.5,2.6);
        AddToChain_MeVnb(V3pi_low,  861,   9.8,2.2);
        AddToChain_MeVnb(V3pi_low,  871,   12.0,2.3);
        AddToChain_MeVnb(V3pi_low,  881,   11.9,2.3);
        AddToChain_MeVnb(V3pi_low,  891,   8.8,2.1);
        AddToChain_MeVnb(V3pi_low,  901,   6.6,1.8);
        AddToChain_MeVnb(V3pi_low,  911,   9.2,2.3);
        AddToChain_MeVnb(V3pi_low,  921,   11.6,2.4);
        AddToChain_MeVnb(V3pi_low,  938,   4.2,1.8);
        AddToChain_MeVnb(V3pi_low,  954,   8.1,0.8);
        AddToChain_MeVnb(V3pi_low,  963,   7.7,2.3);
        AddToChain_MeVnb(V3pi_low,  973,   11.3,3.5);
        AddToChain_MeVnb(V3pi_low,  983,   10.1,2.7);
        AddToChain_MeVnb(V3pi_low,  993,   13.8,5.3);
        AddToChain_MeVnb(V3pi_low,  1003,  25.6,3.7);
//------- (L.M.Barkov et al., Preprint INP 89-15, Novosibirsk, 1989)

        AddToChain_MeVnb(V3pi_low,  840.,13.9,4.3);
        AddToChain_MeVnb(V3pi_low,  880.,6.5,2.8);
        AddToChain_MeVnb(V3pi_low,  940.,9.8,3.3);
        AddToChain_MeVnb(V3pi_low,  984.,10.0,2.5);
        AddToChain_MeVnb(V3pi_low,  1000.,38.,11.);
        AddToChain_MeVnb(V3pi_low,  1010.,53.,18.);


//---------- ND  (S.I.Dolinsky et al., Phys. Rep. 202 (1991) 99)
        
        V3pi_high=CreatePoint2_MeVnb( 1139,  2.2,1.0);
        
        AddToChain_MeVnb(V3pi_high,  1179,  2.3,1.0);
        AddToChain_MeVnb(V3pi_high,  1219,  3.4,0.9);
        AddToChain_MeVnb(V3pi_high,  1259,  2.2,1.1);
        AddToChain_MeVnb(V3pi_high,  1299,  3.3,1.0);
        AddToChain_MeVnb(V3pi_high,  1339,  4.3,1.1);
        AddToChain_MeVnb(V3pi_high,  1379.5,3.8,1.0);

//------- SND

        AddToChain_MeVnb(V3pi_high,  1110,  2.6,0.4);
        AddToChain_MeVnb(V3pi_high,  1120,  3.1,0.5);
        AddToChain_MeVnb(V3pi_high,  1130,  2.8,0.4);
        AddToChain_MeVnb(V3pi_high,  1137.5,2.9,0.5);
        AddToChain_MeVnb(V3pi_high,  1152.5,3.8,0.6);
        AddToChain_MeVnb(V3pi_high,  1162.5,3.2,0.4);
        AddToChain_MeVnb(V3pi_high,  1177.5,4.2,0.5);
        AddToChain_MeVnb(V3pi_high,  1190,  3.7,0.4);
        AddToChain_MeVnb(V3pi_high,  1200,  4.0,0.3);
        AddToChain_MeVnb(V3pi_high,  1210,  4.5,0.4);
        AddToChain_MeVnb(V3pi_high,  1220,  4.7,0.5);
        AddToChain_MeVnb(V3pi_high,  1230,  4.5,0.5);
        AddToChain_MeVnb(V3pi_high,  1240,  4.1,0.4);
        AddToChain_MeVnb(V3pi_high,  1250,  4.0,0.4);
        AddToChain_MeVnb(V3pi_high,  1260,  4.5,0.4);
        AddToChain_MeVnb(V3pi_high,  1270,  3.3,0.3);
        AddToChain_MeVnb(V3pi_high,  1280,  4.1,0.3);
        AddToChain_MeVnb(V3pi_high,  1290,  4.0,0.3);
        AddToChain_MeVnb(V3pi_high,  1300,  3.4,0.3);
        AddToChain_MeVnb(V3pi_high,  1310,  3.5,0.3);
        AddToChain_MeVnb(V3pi_high,  1320,  3.8,0.4);
        AddToChain_MeVnb(V3pi_high,  1330,  3.6,0.3);
        AddToChain_MeVnb(V3pi_high,  1340,  3.7,0.3);
        AddToChain_MeVnb(V3pi_high,  1350,  3.2,0.3);
        AddToChain_MeVnb(V3pi_high,  1360,  3.2,0.3);
        AddToChain_MeVnb(V3pi_high,  1370,  3.5,0.4);
        AddToChain_MeVnb(V3pi_high,  1380,  3.3,0.2);


//------- DM2 (A.Antonelli et al., Z.Phys. C56 (1992) 15.)

        AddToChain_MeVnb(V3pi_high,  1425,  2.59,0.57);
        AddToChain_MeVnb(V3pi_high,  1475,  2.19,0.51);
        AddToChain_MeVnb(V3pi_high,  1525,  1.26,0.47);
        AddToChain_MeVnb(V3pi_high,  1575,  1.94,0.43);
        AddToChain_MeVnb(V3pi_high,  1625,  2.24,0.30);
        AddToChain_MeVnb(V3pi_high,  1675,  1.26,0.26);
        AddToChain_MeVnb(V3pi_high,  1725,  0.64,0.21);
        AddToChain_MeVnb(V3pi_high,  1775,  0.71,0.30);
        AddToChain_MeVnb(V3pi_high,  1825,  0.38,0.17);
        AddToChain_MeVnb(V3pi_high,  1875,  0.50,0.21);
        AddToChain_MeVnb(V3pi_high,  1925,  0.39,0.20);
        AddToChain_MeVnb(V3pi_high,  1975,  0.54,0.30);
        AddToChain_MeVnb(V3pi_high,  2025,  0.16,0.11);
        AddToChain_MeVnb(V3pi_high,  2075,  0.28,0.15);
        AddToChain_MeVnb(V3pi_high,  2125,  0.30,0.14);
        AddToChain_MeVnb(V3pi_high,  2175,  0.23,0.16);
        AddToChain_MeVnb(V3pi_high,  2222.5,0.00,0.10);
        AddToChain_MeVnb(V3pi_high,  2400,  0.25,0.09);

//----------------------------------------------------- \pi^+\pi^-\pi^+\pi^-
      V4pic=CreatePoint2_nb( 0.7 , 0. , 0. );

//-------------------------------------- Саша Суханов

      AddToChain_nb(V4pic,	0.775,	0.036,	0.026);
      AddToChain_nb(V4pic,	0.79,	0.015,	0.011);
      AddToChain_nb(V4pic,	0.81,	0.15,	0.11);
      AddToChain_nb(V4pic,	0.82,	0.14,	0.08);
      AddToChain_nb(V4pic,	0.84,	0.14,	0.07);
      AddToChain_nb(V4pic,	0.88,	0.23,	0.08);
      AddToChain_nb(V4pic,	0.92,	0.37,	0.07);
      AddToChain_nb(V4pic,	0.94,	0.30,	0.09);
      AddToChain_nb(V4pic,	0.95,	0.31,	0.07);
      AddToChain_nb(V4pic,	0.958,	0.53,	0.09);
      AddToChain_nb(V4pic,	0.97,	0.70,	0.11);

//-------------------------------------- Коля Роот

      AddToChain_MeVnb(V4pic,	980,		0.2133997,	0.8624446E-01);
      AddToChain_MeVnb(V4pic,	1045.937,	1.204082,	0.1088124);
      AddToChain_MeVnb(V4pic,	1065.158,	2.132582,	0.1619177);
      AddToChain_MeVnb(V4pic,	1085.832,	2.185026,	0.1847182);
      AddToChain_MeVnb(V4pic,	1106.291,	2.856088,	0.2073246);
      AddToChain_MeVnb(V4pic,	1126.885,	3.183685,	0.2297862);
      AddToChain_MeVnb(V4pic,	1143.212,	4.270751,	0.2490733);
      AddToChain_MeVnb(V4pic,	1164.341,	4.621148,	0.2379854);
      AddToChain_MeVnb(V4pic,	1185.511,	6.132823,	0.2511422);
      AddToChain_MeVnb(V4pic,	1204.142,	8.310534,	0.2575552);
      AddToChain_MeVnb(V4pic,	1225.358,	9.061060,	0.3069048);
      AddToChain_MeVnb(V4pic,	1246.029,	10.16305,	0.2686263);
      AddToChain_MeVnb(V4pic,	1265.784,	11.55784,	0.2644564);
      AddToChain_MeVnb(V4pic,	1285.847,	12.43822,	0.2406759);
      AddToChain_MeVnb(V4pic,	1304.676,	14.54420,	0.2661371);
      AddToChain_MeVnb(V4pic,	1326.509,	14.83099,	0.2584724);
      AddToChain_MeVnb(V4pic,	1345.809,	16.92074,	0.3103572);
      AddToChain_MeVnb(V4pic,	1362.406,	17.10555,	0.3210002);
      AddToChain_MeVnb(V4pic,	1380,		19.79809,	0.3020165);

//--------------------------------------- Слава Шарый
 
      AddToChain_MeVnb(V4pic,	980,		1.07,	0.21);
      AddToChain_MeVnb(V4pic,	984.10,		1.21,	0.09);
      AddToChain_MeVnb(V4pic,	1003.80,	1.26,	0.12);
      AddToChain_MeVnb(V4pic,	1009.68,	1.43,	0.15);
      AddToChain_MeVnb(V4pic,	1015.66,	1.68,	0.16);
      AddToChain_MeVnb(V4pic,	1016.68,	1.60,	0.12);
      AddToChain_MeVnb(V4pic,	1017.66,	1.78,	0.12);
      AddToChain_MeVnb(V4pic,	1018.64,	1.72,	0.11);
      AddToChain_MeVnb(V4pic,	1019.62,	1.96,	0.10);
      AddToChain_MeVnb(V4pic,	1020.58,	1.76,	0.12);
      AddToChain_MeVnb(V4pic,	1021.64,	1.92,	0.17);
      AddToChain_MeVnb(V4pic,	1022.76,	1.72,	0.15);
      AddToChain_MeVnb(V4pic,	1027.74,	1.91,	0.15);
      AddToChain_MeVnb(V4pic,	1033.72,	1.88,	0.17);
      AddToChain_MeVnb(V4pic,	1040.00,	2.09,	0.15);
      AddToChain_MeVnb(V4pic,	1050.00,	2.37,	0.15);
      AddToChain_MeVnb(V4pic,	1060.00,	2.52,	0.13);
      AddToChain_MeVnb(V4pic,	1070.00,	3.11,	0.35);
      AddToChain_MeVnb(V4pic,	1080.00,	4.67,	0.18);
      AddToChain_MeVnb(V4pic,	1090.00,	3.52,	0.38);
      AddToChain_MeVnb(V4pic,	1100.00,	3.92,	0.19);
      AddToChain_MeVnb(V4pic,	1110.00,	4.48,	0.43);
      AddToChain_MeVnb(V4pic,	1120.00,	4.93,	0.25);
      AddToChain_MeVnb(V4pic,	1130.00,	5.45,	0.43);
      AddToChain_MeVnb(V4pic,	1140.00,	6.00,	0.28);
      AddToChain_MeVnb(V4pic,	1150.00,	7.32,	0.62);
      AddToChain_MeVnb(V4pic,	1160.00,	7.33,	0.30);
      AddToChain_MeVnb(V4pic,	1180.00,	8.60,	0.28);
      AddToChain_MeVnb(V4pic,	1190.00,	8.68,	0.42);
      AddToChain_MeVnb(V4pic,	1200.00,	11.57,	0.32);
      AddToChain_MeVnb(V4pic,	1210.00,	10.58,	0.51);
      AddToChain_MeVnb(V4pic,	1220.00,	11.66,	0.36);
      AddToChain_MeVnb(V4pic,	1230.00,	10.84,	0.52);
      AddToChain_MeVnb(V4pic,	1240.00,	13.36,	0.36);
      AddToChain_MeVnb(V4pic,	1250.00,	15.03,	0.50);
      AddToChain_MeVnb(V4pic,	1260.00,	13.79,	0.54);
      AddToChain_MeVnb(V4pic,	1270.00,	16.74,	0.49);
      AddToChain_MeVnb(V4pic,	1280.00,	16.66,	0.50);
      AddToChain_MeVnb(V4pic,	1290.00,	16.21,	0.45);
      AddToChain_MeVnb(V4pic,	1300.00,	18.32,	0.50);
      AddToChain_MeVnb(V4pic,	1310.00,	17.38,	0.53);
      AddToChain_MeVnb(V4pic,	1320.00,	17.28,	0.49);
      AddToChain_MeVnb(V4pic,	1330.00,	18.55,	0.45);
      AddToChain_MeVnb(V4pic,	1340.00,	22.19,	0.46);
      AddToChain_MeVnb(V4pic,	1350.00,	21.66,	0.53);
      AddToChain_MeVnb(V4pic,	1360.00,	23.99,	0.36);
      AddToChain_MeVnb(V4pic,	1370.00,	23.86,	0.55);
      AddToChain_MeVnb(V4pic,	1380,		24.09,	0.40);

//----------------------------------------------------- \pi^+\pi^-\pi^o\pi^o

      V4pin=CreatePoint2_MeVnb( 840, 0., 0.);

//--------------------------------------

      AddToChain_MeVnb(V4pin,	880,	0.58,	0.10);
      AddToChain_MeVnb(V4pin,	920,	0.81,	0.11);
      AddToChain_MeVnb(V4pin,	940,	1.44,	0.17);
      AddToChain_MeVnb(V4pin,	950,	2.12,	0.26);
      AddToChain_MeVnb(V4pin,	958,	2.11,	0.24);
      AddToChain_MeVnb(V4pin,	970,	3.46,	0.31);
      AddToChain_MeVnb(V4pin,	980,	5.61,	0.72);
      AddToChain_MeVnb(V4pin,	1040,	10.7,	0.9);
      AddToChain_MeVnb(V4pin,	1050,	9.7,	0.7);
      AddToChain_MeVnb(V4pin,	1060,	11.4,	0.9);
      AddToChain_MeVnb(V4pin,	1070,	11.4,	0.9);
      AddToChain_MeVnb(V4pin,	1080,	12.0,	1.0);
      AddToChain_MeVnb(V4pin,	1090,	11.8,	0.9);
      AddToChain_MeVnb(V4pin,	1100,	12.4,	1.0);

//------------------------------------- Коля Роот

      AddToChain_MeVnb(V4pin,	980.0000,	3.599568,	0.7674304);
      AddToChain_MeVnb(V4pin,	1045.937,	8.982689,	0.6320198);
      AddToChain_MeVnb(V4pin,	1065.158,	9.041348,	0.7081730);
      AddToChain_MeVnb(V4pin,	1085.832,	10.13851,	0.7965567);
      AddToChain_MeVnb(V4pin,	1106.291,	11.16306,	0.7953351);
      AddToChain_MeVnb(V4pin,	1126.885,	11.44387,	0.8368586);
      AddToChain_MeVnb(V4pin,	1143.212,	12.77563,	0.8246634);
      AddToChain_MeVnb(V4pin,	1164.341,	14.13116,	0.7887244);
      AddToChain_MeVnb(V4pin,	1185.511,	13.74229,	0.7105969);
      AddToChain_MeVnb(V4pin,	1204.142,	15.18844,	0.6566535);
      AddToChain_MeVnb(V4pin,	1225.358,	17.34144,	0.8033268);
      AddToChain_MeVnb(V4pin,	1246.029,	17.93040,	0.6710283);
      AddToChain_MeVnb(V4pin,	1265.784,	19.35792,	0.6438347);
      AddToChain_MeVnb(V4pin,	1285.847,	20.21490,	0.5763935);
      AddToChain_MeVnb(V4pin,	1304.676,	22.26375,	0.6186762);
      AddToChain_MeVnb(V4pin,	1326.509,	20.30975,	0.5678967);
      AddToChain_MeVnb(V4pin,	1345.809,	20.35806,	0.6409006);
      AddToChain_MeVnb(V4pin,	1362.406,	21.43022,	0.6753236);
      AddToChain_MeVnb(V4pin,	1380.000,	25.88614,	0.6487775);

//------------------------------------- Слава Шарый

      AddToChain_MeVnb(V4pin,	980.00,		5.56,	1.58);
      AddToChain_MeVnb(V4pin,	1019.20,	8.07,	1.65);
      AddToChain_MeVnb(V4pin,	1040.00,	9.17,	1.47);
      AddToChain_MeVnb(V4pin,	1050.00,	10.66,	1.35);
      AddToChain_MeVnb(V4pin,	1060.00,	11.45,	0.75);
      AddToChain_MeVnb(V4pin,	1070.00,	10.26,	1.25);
      AddToChain_MeVnb(V4pin,	1080.00,	13.67,	0.63);
      AddToChain_MeVnb(V4pin,	1090.00,	13.01,	1.37);
      AddToChain_MeVnb(V4pin,	1100.00,	14.41,	0.42);
      AddToChain_MeVnb(V4pin,	1110.00,	12.02,	0.86);
      AddToChain_MeVnb(V4pin,	1120.00,	15.67,	0.54);
      AddToChain_MeVnb(V4pin,	1130.00,	15.41,	0.88);
      AddToChain_MeVnb(V4pin,	1140.00,	16.50,	0.58);
      AddToChain_MeVnb(V4pin,	1150.00,	15.78,	1.16);
      AddToChain_MeVnb(V4pin,	1160.00,	17.34,	0.58);
      AddToChain_MeVnb(V4pin,	1180.00,	19.14,	0.54);
      AddToChain_MeVnb(V4pin,	1200.00,	19.54,	0.55);
      AddToChain_MeVnb(V4pin,	1210.00,	19.36,	0.92);
      AddToChain_MeVnb(V4pin,	1220.00,	20.10,	0.64);
      AddToChain_MeVnb(V4pin,	1230.00,	20.34,	1.00);
      AddToChain_MeVnb(V4pin,	1240.00,	22.03,	0.63);
      AddToChain_MeVnb(V4pin,	1250.00,	22.35,	0.85);
      AddToChain_MeVnb(V4pin,	1260.00,	22.33,	0.94);
      AddToChain_MeVnb(V4pin,	1270.00,	24.69,	0.82);
      AddToChain_MeVnb(V4pin,	1280.00,	24.03,	0.83);
      AddToChain_MeVnb(V4pin,	1290.00,	24.46,	0.76);
      AddToChain_MeVnb(V4pin,	1300.00,	26.57,	0.83);
      AddToChain_MeVnb(V4pin,	1310.00,	24.29,	0.87);
      AddToChain_MeVnb(V4pin,	1320.00,	25.92,	0.84);
      AddToChain_MeVnb(V4pin,	1330.00,	26.31,	0.75);
      AddToChain_MeVnb(V4pin,	1340.00,	27.19,	0.63);
      AddToChain_MeVnb(V4pin,	1350.00,	28.81,	0.83);
      AddToChain_MeVnb(V4pin,	1360.00,	28.35,	0.54);
      AddToChain_MeVnb(V4pin,	1370.00,	28.99,	0.81);
      AddToChain_MeVnb(V4pin,	1380.00,	28.80,	0.58);

//-----------------------------------------------------------------------------

      TabR=CreateTab(Vrho);
      AddToTab(TabR,Vhigh);
      AddToTab(TabR,Vkpkm);
      AddToTab(TabR,Vkskl);
      AddToTab(TabR,V3pi_low);
      AddToTab(TabR,V3pi_high);
      AddToTab(TabR,V4pic);
      AddToTab(TabR,V4pin);
      AddToTab(TabR,Vphi);
 //    AddToTab(TabR,V5pi);

//      PrintTab(TabR);
//	PrintChain(Vhigh);
 }


void VP_Done()
 {
  FreeTab(TabR);
 }
 
