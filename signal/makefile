CC = /usr/bin/gcc

default: pysignal

pysignal: setup.py pysignal.pyx libsignal.a
	python setup.py build_ext --inplace && rm -f pysignal.c && rm -rf build

libsignal.a: signal.o
	ar rcs $@ $^
    
signal.o: signal.c signal.h
	$(CC) -c $<

clean:
	rm *.o *.a *.so
