.PHONY: all clean 

all: c f

c: $(patsubst %.c,%.o, $(wildcard ./*.c))
	gcc -o c $^ -lm -L. -lvp
	objcopy --strip-all c
f: $(patsubst %.F,%.o, $(wildcard ./*.F))
	f77 -o f $^ -L. -lvp
	objcopy --strip-all f
clean:
	rm -f *.o c f
