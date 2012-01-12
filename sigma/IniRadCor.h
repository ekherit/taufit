#ifndef IBN_INI_RAD_COR_H
#define IBN_INI_RAD_COR_H

/*
 * Initial states radiative corections.
 * E.A.Kuraev, V.S.Fadin "On radiative corrections to the cross section for signle-photon
 * annihilation of an e^+e^- pair at high energy". Sov.J.Nucl.Phys. 41(3), March 1985
 *
 *  sigma(s) =  int_0^\epsilon dx sigma2(s(1-x)) F(x,s)
 */


/*  Для вычисления интерграла необходимо регуляризовать особенности при  x = 0
 *   и  x = 1, причем посленяя в нашей задачи при s \sim (2m_\tau)^2 не достигается
 *  однако для общности выделил и ее.
 */

inline double F_reg1( double b, double x); // Здесь нет сингулярности
inline double F_reg_x_b( double b, double L); //  особенность \beta x^\beta
inline double F_reg_ln_x( double x);  // особенность ln(x)
inline double F2_reg_x( double x, double s, double b, double L); // поправки второго порядка по альфа

#include "Const.h"

using namespace std;

/* Сингулярная функция, которая и представляет собой  поправки */
inline double F(double x, double s)	
{
				double L=log(s/ME/ME);
				double b=2.*ALPHAPI*(L-1);
				return b * F_reg1(b,x) +
												b * pow(x, b-1) * F_reg_x_b(b,L) +
												b * b * F_reg_ln_x(x) +
												F2_reg_x(x,s,b,L)* sq(ALPHAPI);
}


inline double F_reg1( double b, double x)	{
	return (-1. + 0.5*x + b/8.*(-6. + x));
}

class FReg1	{
	const double b;
	public:
		FReg1(double b_) : b(b_) {}
		double operator()(double x)		{ return F_reg1(b,x); }
};

/* Здесь можноп попытаться соптимизировать -- заменить константы на статические переменные, которые
	 вычисляются лишь один раз */
const double FXB1 = 1. + ALPHAPI*( PI2/3.-0.5);
const double FXB2 = 2*PI2 - 37./4.;
inline double F_reg_x_b( double b, double L )	{
				return (FXB1 + 0.75*b - b*b/24.*(L/3.+ FXB2));
}

class FRegxb	{
	const double b;
	const double L;
	public:
		FRegxb(double b_, double L_) : b(b_), L(L_) {}
		double operator()(double x)		{ return F_reg_x_b(b,L); }
}; 

inline double F_reg_ln_x(double x)	{
		if(x==0) return 0;
		return (-sqrt(x)*( 4*(2.-x)*log(x) + (1.+3.*sq(1.-x)) * log(1.-x)/x));
}

class FRegLn	{
	public:
		FRegLn(void) {}
		double operator()(double x)		{ return F_reg_ln_x(x); }
};

inline double F2_reg_x ( double x, double s, double b, double L)	{
	double cond = x - 4.*ME/sqrt(s);
	if(cond < 0) return 0;
	double L53 = log(s*sq(x/ME)) - 5./3.;
	return (pow(cond,b)*L53*L53*(2.-2.*x+x*x + b*L53/3.)/6. 
		+ 0.5*x*L*L*(2.*(1.-cb(1.-x))/(3.*(1.-x)) + (2.-x)*log(1.-x)+ x*0.5));
}


class FReg2x {
	const double s,b,L;
	public:
		FReg2x(double s_, double b_, double L_ ) :s(s_), b(b_) , L(L_) {}
		double operator()(double x)		{ return F2_reg_x(x,s,b,L); }
};

// Функции для руута
double F_reg1(double *x, double *p)	{
	return F_reg1(*p,*x);
}

double F_reg_ln_x(double *x, double *p)	{
	return F_reg_ln_x(*x);
}

double F2_reg_x(double *x, double *p)	{
	return F2_reg_x(*x,p[0],p[1],p[2]);
} 
#endif 
