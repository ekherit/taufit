#include<stdlib.h>
#include<math.h>
#include"vp.h"

int main(int argc,char *argv[])
 {
  double E,S;
  
  double ReH,ImH,ReL,ImL;
  double E_Re,E_Im;
  double Re,Im;
  double P;

  double D_Re,D_Im,D_ReIm;

  if(argc!=2)
   {
    printf("Usage: ./VP [Energy(GeV)]\n");
    return 1;
   }  

  E=atof(argv[1]); 
  S=E*E;
  if(E<0) S=-S;

  vp_init();  // Инициализируем структуры данных в библиотеке

  vp_lepton(&S,&ReL,&ImL);  // Считаем лептонную часть поляризации
  vp_hadron(&S,&ReH,&ImH);  // Считаем адронную часть поляризации
  vp_hadron_err(&S,&D_Re,&D_Im,&D_ReIm);
  E_Re=sqrt(D_Re);
  E_Im=sqrt(D_Im);

  Re=ReL+ReH;
  Im=ImL+ImH;
  P= (1-Re)*(1-Re) + Im*Im;


  printf("\n %15.5lf  %15.5le %15le %15.5le +/- %15.5le %15le +/- %15.5le ",E,ReL,ImL,ReH,E_Re,ImH,E_Im);  
 
  printf("1/P=%15.5lf        %le\n",1/P,D_ReIm); 

  printf("S=%lf  Re = %le Im=%le\n",S,Re,Im); 
  vp_done(); // Освобождаем структуры данных и память, занятую ими
 }

