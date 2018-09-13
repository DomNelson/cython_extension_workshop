CC=/usr/bin/gcc

default: pymultiply

pymultiply: setup.py pymultiply.pyx libmultiply.a
	python setup.py build_ext --inplace && rm -f pymultiply.c && rm -rf build

libmultiply.a: multiply.o
	ar rcs $@ $^
    
multiply.o: multiply.c multiply.h
	$(CC) -c $<

clean:
	rm *.o *.a *.so
