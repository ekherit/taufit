/*
 * =====================================================================================
 *
 *       Filename:  interpolate.h
 *
 *    Description:  Interpolate
 *
 *        Version:  1.0
 *        Created:  20.04.2018 16:05:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once

#include <map>
#include <array>
#include <iostream>
#include <fmt/printf.h>
#include <TSpline.h>

//#define PRINT_SPLINE_RESULT 

template<class Function>
TSpline3 * interpolate(const Function & f, double xmin, double xmax, double precision=1e-5)
{
  std::map<double, double> cache_table;
  cache_table[xmin] = f(xmin);
  cache_table[xmax] = f(xmax);
  constexpr int N2=8;
  constexpr int N1=N2/2;
  std::array<double, N1+1> x1,y1;
  std::array<double, N2+1> x2,y2;
  auto splineit = [&](double a, double b) -> double
  {
    for(int i=0;i<N2+1;i++) 
    {
      x2[i] =  a+(b-a)/N2*i;
      auto it = cache_table.find(x2[i]);
      if( it == cache_table.end()) 
      {
        //if no cached value then calculate it and add
        y2[i] =  f(x2[i]);
        cache_table[x2[i]] = y2[i];
      }
      else 
      {
        y2[i] = it->second; //just use cache value
      }
    }
    for(int i=0;i<N1+1;i++) 
    {
      x1[i] = x2[2*i];
      y1[i] = y2[2*i];
    }
    double max_diff=0;
    TSpline3 spline("", &x1[0], &y1[0], x1.size()); 
    for(int i=0;i<=N2;i++)
    {
      double diff = fabs(y2[i] - spline.Eval(x2[i]));
      if(diff>max_diff) max_diff = diff;
    }
    return max_diff;
  };

  auto ita = cache_table.begin();
  auto itb = --cache_table.end();
#ifdef PRINT_SPLINE_RESULT
  fmt::printf("%18s %18s %18s %20s %8s\n","a","b","f(a)","Δmax","b-a");
#endif
  while(itb!=cache_table.end())
  {
    double a = ita->first;
    double b = itb->first;
    double max_diff = splineit(a,b);
#ifdef PRINT_SPLINE_RESULT
    fmt::printf("%18.9f %18.9f %18.9f %20.9f %8.3f\n",
        a, 
        b, 
        cache_table[a], 
        max_diff,
        b-a
        );
#endif
    if ( max_diff > precision )
    {
      int i = 0;
      for(int i=0;i<N1;i++) --itb;//move rigt edge to left
    }
    else
    {
      ita = itb;
      for(int i=0;i<N1 && itb!=cache_table.end();i++) ++itb;
    }
  }

  vector<double> X(cache_table.size());
  vector<double> Y(cache_table.size());
  int i = 0;
  for(auto & p : cache_table)
  {
    X[i] = p.first;
    Y[i] = p.second;
    i++;
  }
  auto spline = new TSpline3("spline", &X[0],&Y[0], X.size());

#ifdef PRINT_SPLINE_RESULT
    cout << "Spline result:" << endl;
    fmt::printf("%5s %18s %18s %18s %18s %18s\n", 
        "#",
        "x", 
        "y",
        "spline",
        "Δ",
        "δ");

    for(int i=0;i<X.size();++i)
    {
      fmt::printf("%5d %18.9f %18.9f %18.9f %18f %18f\n", 
          i,
          X[i], 
          Y[i],
          spline->Eval(X[i]),
          Y[i]-spline->Eval(X[i]),
          (Y[i]-spline->Eval(X[i]))/Y[i]
          );
    }
#endif
  return spline;
};
