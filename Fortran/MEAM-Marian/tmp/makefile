PROG = disperx
OBJ =  disper.o input.o forces.o linvlc.o krameam.o rhomeam.o tred2.o tqli.o \
	meamfunctions.o screen.o dscreen.o dscrfor.o

${PROG} : ${OBJ}
	ifort -g -o ${PROG} ${OBJ} 

disper.o : disper.f
	ifort -c -O3 disper.f
input.o : input.f
	ifort -c -O3 input.f
forces.o : forces.f
	ifort -c -O3 forces.f
linvlc.o : linvlc.f
	ifort -c -O3 linvlc.f
krameam.o : krameam.f
	ifort -c -O3 krameam.f
rhomeam.o : rhomeam.f
	ifort -c -O3 rhomeam.f
tred2.o : tred2.f
	ifort -c -O3 tred2.f
tqli.o : tqli.f
	ifort -c -O3 tqli.f
dscrfor.o : dscrfor.f
	ifort -c -O3 dscrfor.f
dscreen.o : dscreen.f
	ifort -c -O3 dscreen.f
screen.o : screen.f
	ifort -c -O3 screen.f
meamfunctions.o : meamfunctions.f
	ifort -c -O3 meamfunctions.f

${OBJ} : moldy.h 
