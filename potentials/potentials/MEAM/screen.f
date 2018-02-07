c#define DEBUG_VECLIST 
*
C This is based on the KRASW subroutine. It is developed to calculate RHO
C Questions/comments >>>> hanchen@llnl.gov
C finished date: Jly 31, 1995   *
*
* Subroutine copied from Mike Baskes' code
* MJC (8/18/00)
*
	subroutine  screen
c	******************
c
        include "moldy.h"
        parameter(scnres=1.e-6)
*
* Initialize data
*
        do i=1,nm
	   do j=1,nnbrs2
              scrnab(i,j) = 0.0
           end do
        end do

       b0lat = b0*alatt

c       start vectorized link cell method

      ix     = 1
      iy     = 1
      iz     = 1
      s      = 1.0

c
c**** Select particles in the current Link cell
c
c******************************************************
c  ibm does not support the following yet
      do 120 ic=1, nlc
      i=ltop(ic)
      if (i.eq.0) goto 599
      m = 0
99    m = m+1
      jaddr(m) = i
      xnbr(m)  = x0(i)
      ynbr(m)  = y0(i)
      znbr(m)  = z0(i)
c
      i        = linkmp(i)
      if (i.gt.0) goto 99
      mstart = m

c
c
c**** Now select all particles from neighboring cells
c
      do 4001 kc = 2,27
c
      sx = 0.0
      sy = 0.0
      sz = 0.0
      jx = ix + nix(kc)
      jy = iy + niy(kc)
      jz = iz + niz(kc)
*
* Remove periodic boundary conditions if necessary
*
cc      if(.not.pbcx) then
cc        if ((ix.eq.nlcx).and.(jx.gt.ix)) goto 4001
cc        if ((ix.eq.1).and.(jx.lt.ix)) goto 4001
cc      end if
cc      if(.not.pbcy) then
cc        if ((iy.eq.nlcy).and.(jy.gt.iy)) goto 4001
cc        if ((iy.eq.1).and.(jy.lt.iy)) goto 4001
cc      end if
cc      if(.not.pbcz) then
cc        if ((iz.eq.nlcz).and.(jz.gt.iz)) goto 4001
cc        if ((iz.eq.1).and.(jz.lt.iz)) goto 4001
cc      end if
      
      if ((ix.eq.nlcx).and.(jx.gt.ix)) then
      jx = 1
      sx = s
      elseif ((ix.eq.1).and.(jx.lt.ix)) then
      jx = nlcx
      sx = -s
      endif
      if ((iy.eq.nlcy).and.(jy.gt.iy)) then
      jy = 1
      sy = s
      elseif ((iy.eq.1).and.(jy.lt.iy)) then
      jy = nlcy
      sy = -s
      endif
      if ((iz.eq.nlcz).and.(jz.gt.iz)) then
      jz = 1
      sz = s
      elseif ((iz.eq.1).and.(jz.lt.iz)) then
      jz = nlcz
      sz = -s
      endif
c
c**** Index of neighbouring cell
c
      jc = jx + nlcx*( (jy-1) + nlcy*(jz-1) ) 
      j  = ltop(jc)
c
c**** Bypass this neighbouring cell if it is empty
c
      if (j.eq.0) goto 4009
199   m = m+1
      jaddr(m) = j
      xnbr(m) = x0(j) + sx 
      ynbr(m) = y0(j) + sy
      znbr(m) = z0(j) + sz
      j = linkmp(j)
      if (j.gt.0) goto 199
c
c**** Save number of particles in first half of link cells
c
4009  if (kc .eq. 14) then
      neihaf = m
      endif
4001  continue
c
c**** We have now found all the neighbouring particles of cell ic.
c
      max1= m
      if (max1.gt.nnbrs) goto 999
c
c	convert to real*8 coordenates
c
	do 89 jg=1,max1
	  x0v = b0lat*xnbr(jg)
	  y0v = b0lat*ynbr(jg)
	  z0v = b0lat*znbr(jg)
	  xnbr(jg)=x0v
	  ynbr(jg)=y0v
	  znbr(jg)=z0v
89	continue
c
c	start calculation
c   calculate pair contribution to total energy and forces
c
c	do on i-particles	**********

	do 6001 im=1,mstart
	iat=jaddr(im)
        idx=id(iat)
C update forces for core atoms only, but need skins as neighbors
	if(iat.gt.nm) goto 6001

	mm=0
c
c	do on j-particles

	do 6002 m=1,max1
	jat=jaddr(m)
        jdx=id(jat)
	if(m.eq.im) goto 6002

c
	dx12t=xnbr(im)-xnbr(m)
        if(abs(dx12t).gt.rcutmeam) goto 6002
	dy12t=ynbr(im)-ynbr(m)
        if(abs(dy12t).gt.rcutmeam) goto 6002
	dz12t=znbr(im)-znbr(m)
        if(abs(dz12t).gt.rcutmeam) goto 6002
c
	rijsq=dx12t*dx12t+dy12t*dy12t+dz12t*dz12t
        if(rijsq.gt.(rcutmeam*rcutmeam)) goto 6002

	rmagg = rijsq
*
* Create list for later calculations
*
        mm = mm + 1
        neighlabel(iat,mm) = jat
*
* Save distance for later calculations
*
cc        neighdist(iat,mm) = rmagg
*
* Find two-body screening -- if complete, forgo next loop
*
        screenij = rscrn(rmagg)

        if(noscr.eq.1.or.screenij.le.scnres) goto 900
ct        if(screenij.le.scnres) goto 900
*
* Three body contribution
*

          dxij = -dx12t
          dyij = -dy12t
          dzij = -dz12t
          rij = rmagg

          do 6004 n2 = 1,max1

             kat=jaddr(n2)
             kdx=id(kat)
             if(kat.eq.jat) goto 6004
             if(kat.eq.iat) goto 6004
c
             dxil=-(xnbr(im)-xnbr(n2))
             if(abs(dxil).gt.rcutmeam) goto 6004
             dyil=-(ynbr(im)-ynbr(n2))
             if(abs(dyil).gt.rcutmeam) goto 6004
             dzil=-(znbr(im)-znbr(n2))
             if(abs(dzil).gt.rcutmeam) goto 6004
c
             ril=dxil*dxil+dyil*dyil+dzil*dzil

             if(ril.gt.(rcutmeam*rcutmeam)) goto 6004

             rjl = rij + ril - 
     >             2.0*(dxij*dxil + dyij*dyil + dzij*dzil)

*
* Three body screening contribution
*
	     screenij = screenij*dscrn(rij,ril,rjl,cmin,cmax)

             if(screenij.le.scnres) goto 900

 6004    continue

 900     continue
*
* Set the screening between pairs of neighbors
*
	 scrnab(iat,mm)=screenij
 

         if(screenij.gt.scnres) rscrmax=max(rscrmax,sqrt(rij))
*
* Move to next atom pair
*

6002  continue

*
* Final number of neighbors for each atom
*

      numneigh(iat) = mm

      if(mm.gt.nnbrs2) then
         print*,'ERROR: nnbrs2 too small'
         call exit(1)
      end if

6001  continue

 599  continue
      
        ix= ix+1
        if ( ix .gt. nlcx ) then
        ix= 1
        iy= iy+1
        if ( iy .gt. nlcy ) then
        iy= 1
        iz=iz+1
        endif
        endif

 120  continue

       return
 999   print*,'ERROR IN NEIGHBORS SCREEN',nnbrs,max1
       end
