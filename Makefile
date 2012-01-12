all : draw_sigma

CXXFLAGS=`root-config --cflags` -I./
LIBS=`root-config --libs` -L./lib  -lvp

draw_sigma :  draw_sigma.o
			g++  $(LIBS) draw_sigma.o -o bin/draw_sigma

draw_sigma.o :  draw_sigma.cpp
			g++ $(CXXFLAGS) -c draw_sigma.cpp -o draw_sigma.o


clean :
				rm *.o
