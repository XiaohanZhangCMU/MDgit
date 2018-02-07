#   makefile
#   $Revision: 1.1 $
#   $Date: 2000/11/20 20:58:21 $

#F77 = g77
#compiler_options = -Ad -O

F77 = ifort
#compiler_options = -O3 
#compiler_options = -O3 -xK -ipo

cc  = icc
#cc_compiler_options = -O3 
#cc_compiler_options = -O3 -xK -ipo

SPEC =

MEAMHOME = ..
UTILSHOME = ../../utils
INCS = -I$(UTILSHOME) -I$(MEAMHOME) $(SPEC)

LIB: utils.o cutils.o acfs.o dyn88a.o forces.o

dyn88:  main.o cutils.o acfs.o dyn88.o  utils.o forces.o dyn88a.o
	$(F77) $(compiler_options) -o dyn88 main.o dyn88.o forces.o utils.o cutils.o acfs.o $(INCS)

main.o:	../main.f
	$(F77) $(compiler_options) -c $(INCS)  ../main.f

template.o:	../template.f
	$(F77) $(compiler_options) -c $(INCS)  ../template.f

acfs.o:         $(UTILSHOME)/acfs.c
	$(cc) $(cc_compiler_options) -c $(INCS) $(UTILSHOME)/acfs.c 

dyn88.o:	$(MEAMHOME)/dyn88.f 
	$(F77) $(compiler_options) -c $(INCS) $(MEAMHOME)/dyn88.f

dyn88a.o:	$(MEAMHOME)/dyn88a.f 
	$(F77) $(compiler_options) -c $(INCS) $(MEAMHOME)/dyn88a.f

utils.o:	$(UTILSHOME)/utils.f 
	$(F77) $(compiler_options) -c  $(INCS) $(UTILSHOME)/utils.f

forces.o:	$(MEAMHOME)/forces.f 
	$(F77) $(compiler_options) -c $(INCS)  $(MEAMHOME)/forces.f

cutils.o:       $(UTILSHOME)/cutils.c
	$(cc) $(cc_compiler_options) -c  $(INCS)  $(UTILSHOME)/cutils.c

clean:
	rm -f *.o dyn88
