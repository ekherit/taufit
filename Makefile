all : draw_sigma taufit tauopt

CXXFLAGS=`root-config --cflags` -I./ -I$(HOME)/work/ -std=c++14 
LIBS= `root-config --libs` -lMinuit -lboost_program_options -lMinuit2 -lfmt

draw_sigma :  draw_sigma.o libvp.a
			g++ draw_sigma.o vp/libvp.a $(LIBS)   -o $(HOME)/work/bin/taudraw

draw_sigma.o :  draw_sigma.cpp interpolate.h
			g++ $(CXXFLAGS) -c draw_sigma.cpp -o draw_sigma.o

taufit :  taufit.o libvp.a
			g++ taufit.o  vp/libvp.a $(LIBS)  -o $(HOME)/work/bin/taufit

taufit.o :  taufit.cpp fit/TauMassFitter.h
			g++ $(CXXFLAGS) -c taufit.cpp -o taufit.o

tauopt.o :  tauopt.cpp
			g++ $(CXXFLAGS) -c tauopt.cpp -o tauopt.o

tauopt :  tauopt.o libvp.a
			g++ tauopt.o  vp/libvp.a $(LIBS)  -o $(HOME)/work/bin/tauopt

libvp.a :
			make -C vp
			
clean :
				rm -f *.o rm *.so* *.a
