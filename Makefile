all : draw_sigma

CXXFLAGS=`root-config --cflags` -I./ -I$(HOME)/work/
LIBS=`root-config --libs` -lMinuit -L./lib  -lvp

draw_sigma :  draw_sigma.o
			g++  $(LIBS) draw_sigma.o -o bin/draw_sigma

draw_sigma.o :  draw_sigma.cpp
			g++ $(CXXFLAGS) -c draw_sigma.cpp -o draw_sigma.o

fit :  fit.o
			g++  $(LIBS) fit.o -o bin/fit

fit.o :  fit.cpp
			g++ $(CXXFLAGS) -c fit.cpp -o fit.o

clean :
				rm *.o
