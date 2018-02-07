c#define DEBUG_VECLIST 
*
C This is based on the KRASW subroutine. It is developed to calculate EAM
C Questions/comments >>>> hanchen@llnl.gov
C finished date: Jly 31, 1995   *
*
* This subroutine computes forces using the input from 
* rhoMEAM.f
* It has been taken for Baskes' code: pwsfor, rhosfor
* MJC (8/14/00)
*
	subroutine kraMEAM
c	******************
c
c
        include "moldy.h"
      
        parameter(scnres=1e-6)
*
        real*8 rcos(3), rr(3), drcos(3,3), dxt, dyt, dzt, rmagg

C replace b0lat by Lx, Ly, Lz
C       b0lat = b0*alatt

*
* Incrment for differentiation
*
         hmeam=1.0d-5
 
*
* Find derivative of the screening function
*
c	reset force vectors

	do 102 i=1,nm
	fx(i)= 0.
	fy(i)= 0.
	fz(i)= 0.
102     continue

	pepair = 0.0
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
*
* From Subroutine rhosfor
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
        if(roz.lt.0) print*,'ROZ=',roz,'I=',iat
        dfrho=dfrhoi(roz,asub,esub)/(z0meam*rozeroi)
*
* d rhobar/d rho0
*
        drbd0=rhotot(iat)/rho0+rho0*dang1(iat)
*
* d rhobar / d AA
*
        drbdAA=rho0*dang2(iat)
*
* Obtain the information about the neighbors of the given atoms
* 
	mm=0
c
c	do on j-particles

	dprcnt = 0.0

        max1 = numneigh(iat)

	do 6002 m=1,max1

	jat=neighlabel(iat,m)


        jdx=id(jat)
	if(iat.eq.jat) goto 6002
*
* I AM NOT SURE ABOUT THIS!!!
*
        scree=scrnab(iat,m)

        if(scree.le.scnres) goto 6002 

	dx12s=xnbr(im)-x0(jat)*Lx
        if(dx12s.ge. Lx/2.) dx12s = dx12s - Lx
        if(dx12s.lt.-Lx/2.) dx12s = dx12s + Lx
	dy12s=ynbr(im)-y0(jat)*Ly
        if(dy12s.ge. Ly/2.) dy12s = dy12s - Ly
        if(dy12s.lt.-Ly/2.) dy12s = dy12s + Ly
	dz12s=znbr(im)-z0(jat)*Lz
        if(dz12s.ge. Lz/2.) dz12s = dz12s - Lz
        if(dz12s.lt.-Lz/2.) dz12s = dz12s + Lz
c
* TRY THE SAME FOR J 	
* Define weighting factors for atomic densities
*
        t1j = tav(1,jat)/c8a(jat)
        t2j = tav(2,jat)/c8a(jat)
        t3j = tav(3,jat)/c8a(jat)
*
* Note the division of rho by z0(rozeroi here ......
*
        rozj=rhotot(jat)/(z0meam*rozeroi)
        rho0j=c8a(jat)
        if(rozj.lt.0) print*,'ROZJ=',rozj,'J=',jat
        dfrhoj=dfrhoi(rozj,asub,esub)/(z0meam*rozeroi)
*
* d rhobar/d rho0
*
        drbd0j=rhotot(jat)/rho0j+rho0j*dang1(jat)
*
* d rhobar / d AA
*
        drbdAAj=rho0j*dang2(jat)
*
* FINISHED TRY
	
	dxt = -dx12s 
	dyt = -dy12s
	dzt = -dz12s

	rijsq=dxt*dxt+dyt*dyt+dzt*dzt
	rmagg =sqrt(rijsq)

        rr(1) = dxt
        rr(2) = dyt
        rr(3) = dzt

        rcos(1) = dxt/rmagg
        rcos(2) = dyt/rmagg
        rcos(3) = dzt/rmagg

*
* Find directions cosines
*
        call dcmij(drcos,rr,rmagg)

*
* Parameters for calculating embedding potentials
*
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
        drhs1=drhof(beta0,re,rhs1)
        drhs2=drhof(beta1,re,rhs2)
        drhs3=drhof(beta2,re,rhs3)
        drhs4=drhof(beta3,re,rhs4)
*
        if(ibar.le.0) then
           factor = 0.0
        else
           factor=(rho0/zmeam)**2
        end if

        ddd=(ts(2)-t1)*(rhsq(1,iat)-s1*factor)+
     >      (ts(3)-t2)*(rhsq(2,iat)-s2*factor)+
     >      (ts(4)-t3)*(rhsq(3,iat)-s3*factor)
        tmp0 = (2.0/3.0)*cg8c(iat)

*
* TRY
*
        dddj=(ts(2)-t1j)*(rhsq(1,jat)-s1*factor)+
     >       (ts(3)-t2j)*(rhsq(2,jat)-s2*factor)+
     >       (ts(4)-t3j)*(rhsq(3,jat)-s3*factor)
        tmp0j = (2.0/3.0)*cg8c(jat)
*
* END TRY
*
        do 380 kk=1,3
           rcosj = rcos(kk)
*
* Initialize accumulation variables
*
           drho0 = drhs1*scree*rcosj
	   drho1 = 0.0
	   drho1g = 0.0
           drho1sg = 0.0
           drho2 = -tmp0*drhs3*scree*rcosj
           drho3 = 0.0
           drho0s = rhs1
           drho1s = 0.0
           drho2s = -rhs3*tmp0           
           drho3s = 0.0


           do 250 km = 1,3
              rcosm = rcos(km)
              tmp1 = 2.0*a8b(km,iat)
              tmp1g = 2.0*ag(km,iat)
              drho1 = drho1 + tmp1*(drhs2*rcosj*rcosm+
     >                rhs2*drcos(km,kk))*scree
              drho1g = drho1g + tmp1g*(drhs4*rcosj*rcosm+
     >                rhs4*drcos(km,kk))*scree
              drho1s = drho1s + tmp1*rhs2*rcosm
              drho1sg = drho1sg + tmp1g*rhs4*rcosm
              do 240 kl=1,3
                 rcosl = rcos(kl)
                 rml = rcosm*rcosl
                 drmllm=drcos(km,kk)*rcosl+
     >                  rcosm*drcos(kl,kk)
                 tmp2 = 2.0*b8c(km,kl,iat)
                 drho2=drho2 + tmp2*(drhs3*rcosj*rml +
     >                 rhs3*drmllm)*scree
                 drho2s = drho2s + tmp2*rhs3*rml
                 do 230 kn=1,3
                    rcosn = rcos(kn)
                    rmln = rml*rcosn
                    tmp3 = 2.0*d8d(km,kl,kn,iat)
                    drho3=drho3+tmp3*(drhs4*rcosj*rmln +
     >                    rhs4*(drmllm*rcosn +
     >                          rml*drcos(kn,kk)))*scree
                    drho3s=drho3s+tmp3*rhs4*rmln
 230             continue
 240          continue
 250       continue

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

           term1 = t1*drho1 + t2*drho2 +
     >             t3*(drho3-legend*drho1g) +
     >             ddd*drho0/rho0


           term2 = t1*drho1s + t2*drho2s +
     >             t3*(drho3s-legend*drho1sg) +
     >             ddd*drho0s/rho0

*
* Compile force on atom iat
*
           tmp4 = dscrnab(kk,iat,m)

         dfa = dfrho*(drho0-
     >                 drho0s*tmp4)*drbd0 
         dfb = dfrho*(term1-term2*tmp4)*drbdAA


            if(kk.eq.1) then
               fx(iat) = fx(iat) + dfa + dfb
            elseif(kk.eq.2) then 
              fy(iat) = fy(iat) + dfa + dfb
            elseif(kk.eq.3) then
              fz(iat) = fz(iat) + dfa + dfb
            else
                  print*,'WHAT THE ...?'
                  call exit(1)
            end if

 380    continue

*
* Change sign to cos ????
*
        rr(1) = -dxt
        rr(2) = -dyt
        rr(3) = -dzt

        rcos(1) = -dxt/rmagg
        rcos(2) = -dyt/rmagg
        rcos(3) = -dzt/rmagg

*
* Find directions cosines
*
        call dcmij(drcos,rr,rmagg)

        
        do 381 kk=1,3
           rcosj = rcos(kk)
*
* TRY
*    
           drho0j = drhs1*scree*rcosj
	   drho1j = 0.0
	   drho1gj = 0.0
           drho1sgj = 0.0
           drho2j = -tmp0j*drhs3*scree*rcosj
           drho3j = 0.0
           drho0sj = rhs1
           drho1sj = 0.0
           drho2sj = -rhs3*tmp0j           
           drho3sj = 0.0


           do 251 km = 1,3
              rcosm = rcos(km)
              tmp1j = 2.0*a8b(km,jat)
              tmp1gj = 2.0*ag(km,jat)
              drho1j = drho1j + tmp1j*(drhs2*rcosj*rcosm+
     >                rhs2*drcos(km,kk))*scree
              drho1gj = drho1gj + tmp1gj*(drhs4*rcosj*rcosm+
     >                rhs4*drcos(km,kk))*scree
              drho1sj = drho1sj + tmp1j*rhs2*rcosm
              drho1sgj = drho1sgj + tmp1gj*rhs4*rcosm
              do 241 kl=1,3
                 rcosl = rcos(kl)
                 rml = rcosm*rcosl
                 drmllm=drcos(km,kk)*rcosl+
     >                  rcosm*drcos(kl,kk)
                 tmp2j = 2.0*b8c(km,kl,jat)
                 drho2j=drho2j + tmp2j*(drhs3*rcosj*rml +
     >                 rhs3*drmllm)*scree
                 drho2sj = drho2sj + tmp2j*rhs3*rml
                 do 231 kn=1,3
                    rcosn = rcos(kn)
                    rmln = rml*rcosn
                    tmp3j = 2.0*d8d(km,kl,kn,jat)
                    drho3j=drho3j+tmp3j*(drhs4*rcosj*rmln +
     >                    rhs4*(drmllm*rcosn +
     >                          rml*drcos(kn,kk)))*scree
                    drho3sj=drho3sj+tmp3j*rhs4*rmln
 231             continue
 241          continue
 251       continue

*
* THIS IS A PATCH UNTIL I FIGURE OUT WHERE THIS CHANGE
* IN SIGN IS COMING FROM. MJC 9/20/00
*
           drho1j   = -drho1j
           drho1gj  = -drho1gj
           drho1sj  = -drho1sj
           drho1sgj = -drho1sgj
           drho3j   = -drho3j
           drho3sj  = -drho3sj

           term1 = t1j*drho1j + t2j*drho2j +
     >             t3j*(drho3j-legend*drho1gj) +
     >             dddj*drho0j/rho0j


           term2 = t1j*drho1sj + t2j*drho2sj +
     >             t3j*(drho3sj-legend*drho1sgj) +
     >             dddj*drho0sj/rho0j

*
* Compile force on atom iat
*
            tmp4 = -dscrnab(kk,iat,m)
*
            df2a = dfrhoj*(drho0j
     >                     -drho0sj*tmp4)*drbd0j 
            df2b = dfrhoj*(term1-term2*tmp4)*drbdAAj

               if(kk.eq.1) then
                  fx(iat) = fx(iat) - df2a - df2b
               elseif(kk.eq.2) then 
                  fy(iat) = fy(iat) - df2a - df2b
               elseif(kk.eq.3) then
                  fz(iat) = fz(iat) - df2a - df2b
            else
                  print*,'WHAT THE ...?'
                  call exit(1)
               end if

381     continue

*
* END TRY
*


*
* Calculate dphif between pairs
*
        phipmeam = phiid(rmagg+hmeam,idx)
        phimmeam = phiid(rmagg-hmeam,idx)
   
        dphif=(phipmeam-phimmeam)/(2.*hmeam)

        phiij=phiid(rmagg,idx)

*
* Compile force on i and jinx atoms
*
        
        dfxv=dphif*rcos(1)*scree+phiij*dscrnab(1,iat,m)
        dfyv=dphif*rcos(2)*scree+phiij*dscrnab(2,iat,m)
        dfzv=dphif*rcos(3)*scree+phiij*dscrnab(3,iat,m)
     
           fx(iat) = fx(iat) - dfxv
           fy(iat) = fy(iat) - dfyv
           fz(iat) = fz(iat) - dfzv

*
*
51	format(1x,2i4,1x,3(g15.8,1x))
c
6002	continue


6001	continue
*
c	end of calculation

599	continue

c 100    continue
c 110    continue
c      timenow = rtc()
c      diff_time = timenow - timestart
c      write (6,*) "diff_time =",
c     .             diff_time

 120    continue
c      l_timenow = rtc()
c      diff_time = l_timenow - l_timestart
c      write (6,*) "diff_time = ",
c     .             diff_time


	pepair = 0.0

	do iat = 1, nm
	  pepair = pepair + atpe(iat)
	end do

        atomepair=pepair
        pe = atomepair	

	atomepair = atomepair/nm
	atome = atomepair 

c	restore units

c 	  include inelastic energy losses if required
ct          if(mdsim.and.iel) call inelas
c
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
	end

