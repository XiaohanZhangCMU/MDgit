	subroutine forces

	include "moldy.h"
	
*** x comes in lattice parameter units

C	do k = 1, nm
C        x0(k) = x(1,k)/b0
C        y0(k) = x(2,k)/b0
C        z0(k) = x(3,k)/b0
C        enddo

*** x0 in unit-box units
C	write(*,*) 'forces: nm=',nm

	  call linvlc
          call screen
          call dscreen
          call rhomeam
          call krameam
          if (noscr.eq.0) call dscrfor

C	do kl = 1, nm
C	   f(1,kl) = fx(kl)
C	   f(2,kl) = fy(kl)
C	   f(3,kl) = fz(kl)
C	enddo

	return
	end
