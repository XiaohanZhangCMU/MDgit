c#define DEBUG_VECLIST 
*
C This is based on the KRASW subroutine. It is developed to calculate RHO
C Questions/comments >>>> hanchen@llnl.gov
C finished date: Jly 31, 1995   *
*
* Subroutine copied from Mike Baskes' code
* MJC (8/18/00)
*
* Forces on screen atoms
*
	subroutine  dscrfor
c	*******************
c
        include "moldy.h"
        parameter(scnres=1.e-6)

        real*8 dil(3),dij(3)

        real*8 rcos(3)

       b0lat = b0*alatt

c
c     Increment for differentiation
c
        hmeam=1.0e-5
        hhmeam = hmeam * hmeam
        h2meam = 2.0 * hmeam

*
*  Initialize the number of visitors to be updated
*
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
      xnbr(m)  = x0(i)*b0lat
      ynbr(m)  = y0(i)*b0lat
      znbr(m)  = z0(i)*b0lat
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
        if(c8a(iat).le.0) goto 6001
*
* Define weighting factors for atomic densities
*
        t1 = tav(1,iat)/c8a(iat)
        t2 = tav(2,iat)/c8a(iat)
        t3 = tav(3,iat)/c8a(iat)
        ibar = ibarr
        ib = abs(ibar)
        asub=asubs
        esub=esubs
        rozeroi=rozros
        zmeam=zsmeam
ct        z0meam = zbar(ibar,zmeam,ncrys,t1,t2,t3)
        z0meam = zmeam
        s1 = 0.0
        s2 = 0.0
        s3 = 0.0
*
* Note the division of rho by z0(rozeroi here ......
*
        roz=rhotot(iat)/(z0meam*rozeroi)
        rho0=c8a(iat)
        if(roz.lt.0) print*,'ROZ=',roz,'i=',iat
        dfrho=dfrhoi(roz,asub,esub)/(z0meam*rozeroi)
*
* d rhobar/d rho0
*
        drbd0=rhotot(iat)/rho0+rho0*dang1(iat)
*
* d rhobar / d AA
*
        drbdAA=rho0*dang2(iat)

        tmp0 = (2.0/3.0) * cg8c(iat)

	mm=0
c
c	do on j-particles

        max1 = numneigh(iat)

	do 6002 m=1,max1
	jat=neighlabel(iat,m)
        jdx=id(jat)
	if(iat.eq.jat) goto 6002

        scree = scrnab(iat,m)
        if(scree.le.scnres.or.scree.ge.1.0-scnres) go to 6002
c
	dx12t=xnbr(im)-x0(jat)*b0lat
        if(dx12t.ge. b0lat/2.) dx12t = dx12t - b0lat
        if(dx12t.lt.-b0lat/2.) dx12t = dx12t + b0lat
	dy12t=ynbr(im)-y0(jat)*b0lat
        if(dy12t.ge. b0lat/2.) dy12t = dy12t - b0lat
        if(dy12t.lt.-b0lat/2.) dy12t = dy12t + b0lat
	dz12t=znbr(im)-z0(jat)*b0lat
        if(dz12t.ge. b0lat/2.) dz12t = dz12t - b0lat
        if(dz12t.lt.-b0lat/2.) dz12t = dz12t + b0lat
c
	rijsq=dx12t*dx12t+dy12t*dy12t+dz12t*dz12t

        rij = rijsq
	rmagg = sqrt(rijsq)

  
        dxij = -dx12t
        dyij = -dy12t
        dzij = -dz12t
*
* Save for later calculations
*
        dij(1) = dxij
        dij(2) = dyij
        dij(3) = dzij
          
        rcos(1) = dxij/rmagg
        rcos(2) = dyij/rmagg
        rcos(3) = dzij/rmagg
        re = res
        rozero=rozros
        beta0 = betas(1)
        beta1 = betas(2)
        beta2 = betas(3)
        beta3 = betas(4)
        rhs1=rhof(rmagg,beta0,re,rozero)
        rhs2=rhof(rmagg,beta1,re,rozero)
        rhs3=rhof(rmagg,beta2,re,rozero)
        rhs4=rhof(rmagg,beta3,re,rozero)
*
        if(ibar.le.0) then
           factor = 0.0
        else
           factor=(rho0/zmeam)**2
        end if

        s1 = 0.0
        s2 = 0.0
        s3 = 0.0

        ddd=(ts(2)-t1)*(rhsq(1,iat)-s1*factor)+
     >      (ts(3)-t2)*(rhsq(2,iat)-s2*factor)+
     >      (ts(4)-t3)*(rhsq(3,iat)-s3*factor)

        phiij = phiid(rmagg,idx)

        do 6004 n2 = 1,max1

             kat=neighlabel(iat,n2)
             kdx=id(kat)
             if(kat.eq.jat) goto 6004
             if(kat.eq.iat) goto 6004
c
             dxil=-(xnbr(im)-x0(kat)*b0lat)
             if(dxil.ge. b0lat/2.) dxil = dxil - b0lat
             if(dxil.lt.-b0lat/2.) dxil = dxil + b0lat
             dyil=-(ynbr(im)-y0(kat)*b0lat)
             if(dyil.ge. b0lat/2.) dyil = dyil - b0lat
             if(dyil.lt.-b0lat/2.) dyil = dyil + b0lat
             dzil=-(znbr(im)-z0(kat)*b0lat)
             if(dzil.ge. b0lat/2.) dzil = dzil - b0lat
             if(dzil.lt.-b0lat/2.) dzil = dzil + b0lat
*
* Save for later calculations
*
             dil(1) = dxil
	     dil(2) = dyil
	     dil(3) = dzil
c
             rilsq=dxil*dxil+dyil*dyil+dzil*dzil

             ril = rilsq

             rjl = rij + ril - 
     >             2.0*(dxij*dxil + dyij*dyil + dzij*dzil)
c
c    Find proper coordinates for compilation of slocal
c
              dxlm = 2.0 * dxil - dxij
              dylm = 2.0 * dyil - dyij
              dzlm = 2.0 * dzil - dzij
c
c     Find screening
c
              scr1=dscrn(rij,ril,rjl,cmin,cmax)
c
c
c     See if derivative of screening term is non-zero
c
              if (scr1.lt.scnres.or.scr1.gt.1.0-scnres) goto 6004

              do 3400 kk = 1 , 3
*
c    Incrementing ith atom position
c
                rpil = ril + h2meam * dil(kk) + hhmeam
                rpjl = rjl + h2meam * (dil(kk)-dij(kk)) + hhmeam
                scrpl=dscrn(rij,rpil,rpjl,cmin,cmax)
c
c    negative increment in xi's
c
                rpil = ril - h2meam * dil(kk) + hhmeam
                rpjl = rjl - h2meam * (dil(kk)-dij(kk)) + hhmeam
                scrml=dscrn(rij,rpil,rpjl,cmin,cmax)
c
c     Find derivative
c
                dscrnabk=0.5*((scrpl-scrml)*scree)
     >                    /(scr1*hmeam)

                rcosj = rcos(kk)
c
c     Initialize accumulation variables
c
                drho0s = rhs1
                drho1s = 0.0
                drho1sg = 0.0
                drho2s = -rhs3 * tmp0
                drho3s = 0.0
                do 250 km=1,3
                  rcosm = rcos(km)
                  tmp1 = 2.0 * a8b(km,iat)
                  tmp1g = 2.0 * ag(km,iat)
                  drho1s = drho1s + tmp1 * rhs2 * rcosm
                  drho1sg = drho1sg + tmp1g * rhs4 * rcosm
                  do 240 ll = 1, 3
                    rcosl = rcos(ll)
                    rml = rcosm * rcosl
                   tmp2 = 2.0 * b8c(km,ll,iat)
                    drho2s = drho2s + tmp2 * rhs3 * rml
                    do 230 n = 1, 3
                      rcosn = rcos(n)
                      rmln = rml * rcosn
                      tmp3 = 2.0 * d8d(km,ll,n,iat)
                      drho3s = drho3s + tmp3 * rhs4 * rmln
 230                continue
 240              continue
 250            continue

*
* THIS IS A PATCH UNTIL I FIGURE OUT WHERE THIS CHANGE
* IN SIGN IS COMING FROM. MJC 9/20/00
*
                drho1   = -drho1
                drho1g  = -drho1g
                drho1s  = -drho1s
                drho1sg = -drho1sg
                drho3   = -drho3
                drho3s  = -drho3s
c
                term2 = t1*drho1s + t2*drho2s + t3*(drho3s
     1               -legend*drho1sg)
     1                + ddd * drho0s / rho0
c
c     Sum forces on the central atom
c
                df = dscrnabk *
     .               (0.5 * phiij +
     .                dfrho * (drbd0 * drho0s + drbdAA * term2))


                if(kk.eq.1) then
                   fx(kat) = fx(kat) - df
                elseif(kk.eq.2) then
                   fy(kat) = fy(kat) - df
                elseif(kk.eq.3) then
                   fz(kat) = fz(kat) - df
                end if
         
3400          continue

6004      continue
*
* Move to next atom pair
*

6002   continue

6001   continue
599    continue
120    continue

       return
       end
