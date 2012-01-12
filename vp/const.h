#ifndef _CONST_INC_
#define _CONST_INC_
//------------------------ Всякие разные константы

#define pi		3.1415926535897  // пи
#define alpha		(1./137.0359976)      // альфа
#define	alpi		(alpha/pi)     
#define	GeV2nb		0.389379292E6  // коэф. перевода сечения из GeV^2 в nb 
#define	nb2GeV		(1./GeV2nb)    

#define	me		0.000510998902  // масса электрона
#define	mmu		0.105658357   // масса мюона
#define	mpi		0.13957018   // масса пиона
#define	mtau		1.77699   // масса тау     
                           
#define	zeta3		1.20205690316 //  Дзета-функция Римана от 3  
                              
#define	treshold	(2*mpi)  // Порог рождения 2 pi
#define	treshold2	(4*mpi*mpi)

#define	deltaS		0.001
#define	Max		11.   // Максимальная энергия в ГэВ для работы программы
#define	Max2		(Max*Max)

#endif                                             
