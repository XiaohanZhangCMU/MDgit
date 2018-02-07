c#define DEBUG_VECLIST 
*
C This is based on the KRASW subroutine. It is developed to calculate RHO
C Questions/comments >>>> hanchen@llnl.gov
C finished date: Jly 31, 1995   *
*
* Subroutine copied from Mike Baskes' code
* MJC (8/18/00)
*
	subroutine  dscreen
c	*******************
c
        include "moldy.h"
        parameter(scnres=1.d-6)

C replace b0lat by Lx, Ly, Lz
C        b0lat = b0*alatt

c
c     Increment for differentiation
c
        hmeam=1.0e-5
        hhmeam = hmeam * hmeam
        h2meam = 2.0 * hmeam
*
*  Initialize
*
        do i=1,nm
            do k=1,3
               do j=1,nnbrs2
                  dscrnab(k,i,j) = 0.0
               end do
            end do
        end do
c
c**** Select particles in the current Link cell
c
c******************************************************
      do 120 ic=1, nlc
      i=ltop(ic)
      if (i.eq.0) goto 599
      m = 0
99    m = m+1
      jaddr(m) = i
      xnbr(m)  = x0(i)*Lx
      ynbr(m)  = y0(i)*Ly
      znbr(m)  = z0(i)*Lz
c
      i        = linkmp(i)
      if (i.gt.0) goto 99
      mstart = m

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

        max1 = numneigh(iat)

	do 6002 m=1,max1
	jat=neighlabel(iat,m)
        jdx=id(jat)
	if(iat.eq.jat) goto 6002

c
	dx12t=xnbr(im)-x0(jat)*Lx
        if(dx12t.ge. Lx/2.) dx12t = dx12t - Lx
        if(dx12t.lt.-Lx/2.) dx12t = dx12t + Lx
	dy12t=ynbr(im)-y0(jat)*Ly
        if(dy12t.ge. Ly/2.) dy12t = dy12t - Ly
        if(dy12t.lt.-Ly/2.) dy12t = dy12t + Ly
	dz12t=znbr(im)-z0(jat)*Lz
        if(dz12t.ge. Lz/2.) dz12t = dz12t - Lz
        if(dz12t.lt.-Lz/2.) dz12t = dz12t + Lz
c
	rijsq=dx12t*dx12t+dy12t*dy12t+dz12t*dz12t

	rmagg = rijsq
        dxij=-dx12t
        dyij=-dy12t
        dzij=-dz12t
        rij = rmagg

        scree=scrnab(iat,m)
c
c    If not, set derivatives to zero and exit from the neighbor loop
c
        if (scree .gt. 1.0-scnres .or. scree .lt. scnres) then
              do 105 km=1,3
                dscrnab(km,iat,m)=0.0
 105          continue
              goto 6002
           endif
c
c  Perform numerical differenciation
c
           scrpx = rscrn(rij - h2meam * dxij + hhmeam)
           scrpy = rscrn(rij - h2meam * dyij + hhmeam)
           scrpz = rscrn(rij - h2meam * dzij + hhmeam)

           scrmx = rscrn(rij + h2meam * dxij + hhmeam)
           scrmy = rscrn(rij + h2meam * dyij + hhmeam)
           scrmz = rscrn(rij + h2meam * dzij + hhmeam)

           if(noscr.eq.1) goto 6005

           do 6004 n2 = 1,max1

             kat=neighlabel(iat,n2)

             kdx=id(kat)
             if(kat.eq.jat) goto 6004
             if(kat.eq.iat) goto 6004
c
             dxil=-(xnbr(im)-x0(kat)*Lx)
             if(dxil.ge. Lx/2.) dxil = dxil - Lx
             if(dxil.lt.-Lx/2.) dxil = dxil + Lx
             dyil=-(ynbr(im)-y0(kat)*Ly)
             if(dyil.ge. Ly/2.) dyil = dyil - Ly
             if(dyil.lt.-Ly/2.) dyil = dyil + Ly
             dzil=-(znbr(im)-z0(kat)*Lz)
             if(dzil.ge. Lz/2.) dzil = dzil - Lz
             if(dzil.lt.-Lz/2.) dzil = dzil + Lz
c
             ril=dxil*dxil+dyil*dyil+dzil*dzil

             rjl = rij + ril - 
     >             2.0*(dxij*dxil + dyij*dyil + dzij*dzil)
*
c    Incrementing ith atom position
c
             xpij = rij - h2meam * dxij + hhmeam
             xpil = ril - h2meam * dxil + hhmeam
             scrpx=scrpx*dscrn(xpij,xpil,rjl,cmin,cmax)

             ypij = rij - h2meam * dyij + hhmeam
             ypil = ril - h2meam * dyil + hhmeam
             scrpy=scrpy*dscrn(ypij,ypil,rjl,cmin,cmax)

             zpij = rij - h2meam * dzij + hhmeam
             zpil = ril - h2meam * dzil + hhmeam
             scrpz=scrpz*dscrn(zpij,zpil,rjl,cmin,cmax)
c
c    negative increment in xi's
c
             xpij = rij + h2meam * dxij + hhmeam
             xpil = ril + h2meam * dxil + hhmeam
             scrmx=scrmx*dscrn(xpij,xpil,rjl,cmin,cmax)

             ypij = rij + h2meam * dyij + hhmeam
             ypil = ril + h2meam * dyil + hhmeam
             scrmy=scrmy*dscrn(ypij,ypil,rjl,cmin,cmax)

             zpij = rij + h2meam * dzij + hhmeam
             zpil = ril + h2meam * dzil + hhmeam
             scrmz=scrmz*dscrn(zpij,zpil,rjl,cmin,cmax)

 6004    continue
 6005    continue
c
c    find and store derivative
c
          dscrnab(1,iat,m)=0.5*(scrpx-scrmx)/hmeam
          dscrnab(2,iat,m)=0.5*(scrpy-scrmy)/hmeam
          dscrnab(3,iat,m)=0.5*(scrpz-scrmz)/hmeam
*
* Move to next atom pair
*

6002  continue

6001  continue
 599  continue
 120  continue

       return
       end
