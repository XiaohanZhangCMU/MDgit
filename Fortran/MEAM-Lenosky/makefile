
F90 = ifort

F90FLAGS = -w -xW -g
#F90FLAGS = -w -xW -O2
#F90FLAGS = -w -xW -O2 -openmp

all: test

%.o:%.f90
	$(F90) -c $(F90FLAGS) $^

test: main.o lenosky_lib.o bazant_lib.o
	$(F90) -o $@ $^

clean:
	rm -f test *.o
