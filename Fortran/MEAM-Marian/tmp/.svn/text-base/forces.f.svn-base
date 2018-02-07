	subroutine forces

	include "moldy.h"
	
*** x comes in lattice parameter units

	do k = 1, nm
        x0(k) = x(1,k)/b0
        y0(k) = x(2,k)/b0
        z0(k) = x(3,k)/b0
        enddo

*** x0 in unit-box units

	  call linvlc
          call screen
          call dscreen
          call rhomeam
          call krameam
          if (noscr.eq.0) call dscrfor

	do kl = 1, nm
	   f(1,kl) = fx(kl)
	   f(2,kl) = fy(kl)
	   f(3,kl) = fz(kl)
	enddo

	return
	end
