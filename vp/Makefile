Name   := libvp.so
CFLAGS := -fPIC

.PHONY: all clean cleanall

all: $(Name) libvp.a
$(Name): $(patsubst %.c, %.o, $(wildcard ./*.c))
	gcc -shared -Wl,-soname,$@.1 -o $@.1.0.1 $^ -lc -lm
	objcopy --strip-debug $@.1.0.1
	chmod -x $@.1.0.1
	/sbin/ldconfig -n .
	ln -sf $@.1 $@

libvp.a : $(wildcard *.o)
	ar rcs $@ *.o

%.o: %.c 
	gcc -c -MD $(CFLAGS) $<

include $(wildcard *.d)

cleanall: clean
	rm -f *.d
clean: 
	rm -f *.o $(Name)* libvp.a
		
