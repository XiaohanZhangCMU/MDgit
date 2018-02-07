	subroutine input
*       ****************       *	

	include "moldy.h"

*****************************************************************************
C Read initial atomic relaxed configuration
*****************************************************************************

	open (unit = 10, file = 'conf.dat', status = 'old')
*** file 'conf.dat' contains relaxed atomic positions

	do n = 1, nm
	read(10,*) jdummy, id(n), (xi(m,n), m = 1, 3)
*** relaxed atomic positions loaded in array xi
	enddo

	x = xi
	close (10)

	call meamdefine

	return
	end
