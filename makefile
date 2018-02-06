# /usr/home/caiwei/Codes/MD++/makefile
# Creation Date : Mon Oct 20 17:56:29 2003
# Last Modified : Tue Dec 13 23:30:08 2005

SYS = gpp
build = D
BIN = bin
MAKE = make
SPEC = 

help:
	cat README

run:
	bin/fs_$(SYS) Examples/example01-mo.script

all: fs sw eam lj


%:
	cd src; $(MAKE) $@ SYS=$(SYS) build=$(build) BIN=$(BIN) SPEC=$(SPEC); cd ..

clean:
	rm -f bin/*_*
	cd src; $(MAKE) clean; cd ..
	cd Fortran/MEAM-Marian; $(MAKE) clean;
	cd Fortran/MEAM-Baskes/meam/linux; $(MAKE) clean;
	cd Fortran/MEAM-Lenosky; $(MAKE) clean;
	cd Fortran/MEAM-Lammps; $(MAKE) clean;
	cd src/StencilToolkit; $(MAKE) clean;
