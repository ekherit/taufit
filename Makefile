all : draw_sigma taufit

CXXFLAGS=`root-config --cflags` -I./ -I$(HOME)/work/ -std=c++0x
LIBS= `root-config --libs` -lMinuit -lboost_program_options -lMinuit2

draw_sigma :  draw_sigma.o libvp.a
			g++ draw_sigma.o vp/libvp.a $(LIBS)   -o $(HOME)/work/bin/taudraw

draw_sigma.o :  draw_sigma.cpp
			g++ $(CXXFLAGS) -c draw_sigma.cpp -o draw_sigma.o

taufit :  taufit.o libvp.a
			g++ taufit.o  vp/libvp.a $(LIBS)  -o $(HOME)/work/bin/taufit

taufit.o :  taufit.cpp fit/TauMassFitter.h
			g++ $(CXXFLAGS) -c taufit.cpp -o taufit.o

libvp.a :
			make -C vp
			
clean :
				rm -f *.o rm *.so* *.a
