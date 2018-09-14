## Format here is
##
## output: input(s)
##     Command to turn inputs to outputs

## I need this to avoid using the anaconda compiler
CC=/usr/bin/gcc

default: pymultiply


## Use distutils and cython to build extension. Important to clean up after
## so that the extension is rebuilt every time the source is changed. If
## anaconda is complicating things, you can use the system python here, 
## something like `/usr/bin/python setup.py ...`

pymultiply: setup.py pymultiply.pyx libmultiply.a
	python setup.py build_ext --inplace && rm -f pymultiply.c && rm -rf build


## Add object file to library libmultiply. $@ is shorthand for the first target
## (output) and $^ is shorthand for all prerequisites (inputs), so this is
## equivalent to `ar rcs libmultiply.a multiply.o`

libmultiply.a: multiply.o
	ar rcs $@ $^
    

## Compile C source code to create object file. Flag -fPIC allows memory to
## be dynamically allocated, which is needed for (importable) shared objects
## on some systems. '$<' is shorthand for the first prerequisite (input),
## so this is equivalent to `/usr/bin/gcc -c -fPIC multiply.c`

multiply.o: multiply.c multiply.h
	$(CC) -c -fPIC $<


## Clean up all our build files

clean:
	rm *.o *.a *.so
