FLAGS=-O3 -std=c++11 -I../

all:
	g++ -c -fPIC LocalBoundary1D.cpp LocalBoundary2D.cpp LocalBoundary3D.cpp
	mpic++ -c -fPIC MachineInfo.cpp RemoteBoundary1D.cpp RemoteBoundary2D.cpp RemoteBoundary3D.cpp
	mpic++ -shared -fPIC MachineInfo.o LocalBoundary1D.o LocalBoundary2D.o LocalBoundary3D.o RemoteBoundary1D.o RemoteBoundary2D.o RemoteBoundary3D.o -o libstk.so
	cp libstk.so ../
	mv libstk.so ../../bin

clean:
	rm -f *.o ../libstk.so ../../bin/libstk.so
