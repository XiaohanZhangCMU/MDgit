	subroutine input
*       ****************       *	

	include "moldy.h"

*****************************************************************************
C Read initial atomic relaxed configuration
*****************************************************************************

C	open (unit = 10, file = 'conf.dat', status = 'old')
	open (unit = 10, 
     >	file = '/home/caiwei/Codes/MD++/Fortran/MEAM-Marian/conf.dat', 
     >  status = 'old')
*** file 'conf.dat' contains relaxed atomic positions

	nm = 256

	do n = 1, nm
	read(10,*) jdummy, id(n), (xi(m,n), m = 1, 3)
*** relaxed atomic positions loaded in array xi
	enddo

	x = xi

	b0=4.0
	do k = 1, nm
        x0(k) = x(1,k)/b0
        y0(k) = x(2,k)/b0
        z0(k) = x(3,k)/b0
        enddo

	alatt=4.640
	Lx=b0*alatt
	Ly=b0*alatt
	Lz=b0*alatt

	nlcx = 4
	nlcy = 4
	nlcz = 4
	nlc = nlcx*nlcy*nlcz

	close (10)

C	meamfile = 'meamdata'
	meamfile = '/home/caiwei/Codes/MD++/Fortran/MEAM-Marian/meamdata.Pu'

	call meamdefine

	return
	end
