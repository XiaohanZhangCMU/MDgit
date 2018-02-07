c
	subroutine rhoMEAM
c	******************
c
        include "moldy.h"

        real*8 rcosv(3)

C replace b0lat by Lx, Ly, Lz
C        b0lat = b0*alatt

        petrip=0.
        rhocon= 1.e-10
c
	do i  = 1,nm
          atpe3b(i)=0.0
	  rhotot(i) = 0.0
          atpe(i) = 0.0
          embf(i) = 0.0
          c8a(i) = 0.0
          cg8c(i) = 0.0
          do j=1,3
             tav(j,i) = 0.0
	     ag(j,i) = 0.0
	     a8b(j,i) = 0.0
	     do k=1,3
                b8c(k,j,i) = 0.0
                do n=1,3
                   d8d(n,k,j,i) = 0.0
                end do
             end do
          end do
	end do

c	reset force vectors

	do 102 i=1,nm
	fx(i)= 0.
	fy(i)= 0.
102	fz(i)= 0.

	petrip = 0.0
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
c	if(iat.gt.nm) goto 6001

clc	mm=0
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

	rmagg =sqrt(rijsq)

	if (rmagg.eq.0.0) then
	print *, 'xi, yi, zi',iat, x0(iat)*b0, y0(iat)*b0, z0(iat)*b0
	print *, 'xj, yj, zj',jat, x0(jat)*b0, y0(jat)*b0, z0(jat)*b0
	endif


        scree = scrnab(iat,m)
*
*  Compute the pair potential phi(r), Eq. (4)
*

        phifs = phiid(rmagg,idx)

        atpe(iat) = atpe(iat) + 0.5*phifs*scree
	
C	write(*,*) 'atpe(',iat,')=',atpe(iat)
*
*  Compute the electron density: Eq. 8(a-d), 9(a)
*

        rej = res
        rzj = rozros
        betaj0 = betas(1)
        betaj1 = betas(2)
        betaj2 = betas(3)
        betaj3 = betas(4)
        rhj = rhof(rmagg,betaj0,rej,rzj)
        rogj1 = rhof(rmagg,betaj1,rej,rzj)
        rogj2 = rhof(rmagg,betaj2,rej,rzj)
        rogj3 = rhof(rmagg,betaj3,rej,rzj)

        rinve = 1./rmagg
*
*  WARNING!!! CHECK THE SIGN!!!
*
        rcosv(1) = -rinve*dx12t
        rcosv(2) = -rinve*dy12t
        rcosv(3) = -rinve*dz12t

        do itmp = 1, 3
           t1 = scree*rcosv(itmp)
*
* Eq. 8(b), sum over j, sum(xij*rho1j(R))
*
           a8b(itmp,iat) = a8b(itmp,iat) - rogj1*t1

           ag(itmp,iat) = ag(itmp,iat) - rogj3*t1

*
* Eq. 8(c), sum over j, sum(xij*xij*rho2j(r))
*
           do jtmp = 1,3

              t2=t1*rcosv(jtmp)
              b8c(itmp,jtmp,iat) = b8c(itmp,jtmp,iat)+rogj2*t2
*
* Eq. 8(d), sum over j, sum(xij*xij*xij*rho3j(r))
*
              do ktmp = 1,3

                 d8d(itmp,jtmp,ktmp,iat) = d8d(itmp,jtmp,ktmp,iat) 
     >                                   - rogj3*t2*rcosv(ktmp)
              end do
           end do
        end do

*
* Eq. 8(a), rho(0) = sum(rho0j(r))

        c8a(iat) = c8a(iat) + rhj*scree
*
* USED IN THE CALCULATION OF FORCES. PASS THIS INFORMATION TO KRAMEMA.F
* PASS ALSO c8a(iat) which is rho0
*
        tav(1,iat) = tav(1,iat) + ts(2)*rhj*scree
        tav(2,iat) = tav(2,iat) + ts(3)*rhj*scree
        tav(3,iat) = tav(3,iat) + ts(4)*rhj*scree
*
* Save for second term in Eq. 8(c)
*
        cg8c(iat) = cg8c(iat) + rogj2*scree

c
6002    continue
*
6001	continue
c

c	end of calculation

c	link cell index updat

599	continue


 120    continue

*
* NOW FIND F(rho) for EACH PARTICLE
*

	do iat = 1, nm

          rho0 = c8a(iat)

          if(rho0.gt.0) then
	     idx = id(iat)
             t1=tav(1,iat)/rho0
             t2=tav(2,iat)/rho0
             t3=tav(3,iat)/rho0
             ibar = ibarr
             asub = asubs
             esub = esubs
             rozeroi = rozros
             zmeam = zsmeam

ct             z0meam = zbar(ibar,zmeam,ncrys,t1,t2,t3)
             z0meam = zmeam
*
* Eqns  (8b), (8c), (8d)
* Sum over the coordinates to get rho1i, rho2i, rho3i  
*
             rho1 = 0.0
*
* Second term in Eq. 8(c)
*
             rho2 = -(cg8c(iat)**2)/3.0

	     rho3 = 0.0
*
* This is different from paper Phys. Rev. B46, (1992), 2727
* CHECK
*
             do 1399 itmp=1,3
                rho3 = rho3 - legend*ag(itmp,iat)**2             
 1399        continue
        
             do 1400 itmp=1,3

                rho1 = rho1 + a8b(itmp,iat)**2
 
                do 1400 jtmp=1,3

                   rho2 = rho2 + b8c(itmp,jtmp,iat)**2

                   do 1400 ktmp=1,3

                      rho3 = rho3 + d8d(itmp,jtmp,ktmp,iat)**2
 1400        continue

*
* Save values for later calculation of forces (I THINK!!)
*
             rhsq(1,iat) = rho1
             rhsq(2,iat) = rho2
             rhsq(3,iat) = rho3
*
* Eqn. 9(a). Sum over three terms to compute rhobari
*
             drh = t1*rho1+t2*rho2+t3*rho3
             rhobr = bar(rho0,drh,ibar,zmeam*rozeroi,
     >                   dang1(iat),dang2(iat)) 

             rhobar = rhobr/(z0meam*rozeroi)

             if(rhobr.lt.0) then
                print*,'RHO IS NEGATIVE',rhobr, 
     >               'ATOM =',iat
             end if

             rhotot(iat) = rhobr
             
* F(rho)

             embf(iat) = frhoi(rhobar,asub,esub)
*
* Energy for each atom
*
             atpe(iat) = atpe(iat) + embf(iat)                          
	     petrip = petrip + embf(iat)
	     atpe3b(iat) = embf(iat)

             if(rhotot(iat).lt.0) then
                print*,"WARNING: FOR ATOM ",iat,
     >                 "RHOBAR=",rhobar,embf(iat)
             end if
*
* if rho.lt.0
* 
	  else
             embf(iat) = 0.0
             rhotot(iat) = 1.e-10
	     rhsq(1,iat) = 0.0
	     rhsq(2,iat) = 0.0
	     rhsq(3,iat) = 0.0
           end if
*
* Finish loop over all atoms
*
	end do

	return
999     continue
c        if(mypid.eq.0) then
           write(101,1004)max1,nnbrs
           write(101,1005)mm,nnbrs1
           stop
c        endif
1004	format('stopping because nnbrs is too smal',
     #	/,'max1= ',i6,'  nnbrs= ',i6)
1005	format('     or because nnbrs1 is too smal',
     #	/,'mm= ',i6,'    nnbrs1= ',i6)
998   write(101,1006) iat,rhotot(iat),rhocon
1006  format(' node',i4,': stopping in rhovlc because rho is ',
     >       'too small',/
     >       ' ***** global ID =',i4,' rhotot(i) =',e17.8,' <',e17.8)
1001    write(101,*)"index for r out of bounds",index

	end

