all : draw_sigma taufit tauopt tausim

CXXFLAGS=`root-config --cflags` -I./ -I$(HOME)/work/ -std=c++14
LIBS= `root-config --libs` -lMinuit -lboost_program_options -lMinuit2 -lfmt 
#-L$(HOME)/local/lib -lfmt 

draw_sigma :  draw_sigma.o libvp.a
			g++ draw_sigma.o vp/libvp.a $(LIBS)   -o $(HOME)/work/bin/taudraw 

draw_sigma.o :  draw_sigma.cpp interpolate.h
			g++ $(CXXFLAGS) -c draw_sigma.cpp -o draw_sigma.o


taufit :  taufit.o libvp.a draw.o
			g++ taufit.o  draw.o vp/libvp.a $(LIBS)  -o $(HOME)/work/bin/taufit

taufit.o :  taufit.cpp fit/TauMassFitter.h draw.h
			g++ $(CXXFLAGS) -c taufit.cpp -o taufit.o

tauopt.o :  tauopt.cpp 
			g++ $(CXXFLAGS) -c tauopt.cpp -o tauopt.o

tauopt :  tauopt.o libvp.a
			g++ tauopt.o  vp/libvp.a $(LIBS)  -o $(HOME)/work/bin/tauopt

draw.o :  draw.cpp fit/TauMassFitter.h ScanPoint.h
			g++ $(CXXFLAGS) -c draw.cpp -o draw.o

tausim.o :  tausim.cpp 
			g++ $(CXXFLAGS) -c tausim.cpp -o tausim.o

tausim :  tausim.o libvp.a draw.o
			g++ tausim.o  draw.o vp/libvp.a $(LIBS)  -o $(HOME)/work/bin/tausim

libvp.a :
			make -C vp
			
clean :
				rm -f *.o rm *.so* *.a
