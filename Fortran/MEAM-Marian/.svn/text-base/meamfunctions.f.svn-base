* I will use only one type of atom!!!
* Defining variables for MEAM
*
       subroutine meamdefine

       include "moldy.h"

*      Read values from file

C        open(43,file='meamdata',status='old')
        open(43,file=meamfile,status='old')

        read(43,*)zsmeam
        read(43,*)alphas
        read(43,*)betas(1),betas(2),betas(3),betas(4)
        read(43,*)esubs,asubs
        read(43,*)ts(1),ts(2),ts(3),ts(4)
        read(43,*)rozros
        read(43,*)ibarr
        read(43,*)rcutmeam
        read(43,*)cmin,cmax
        read(43,*)repuls, attrac
        read(43,*)legend
        read(43,*)noscr

        read(43,*)ncrys
        read(43,*)equidist

        close(43)
       
        rcutmeam2 = rcutmeam*rcutmeam

        print*,zsmeam
        print*,alphas
        print*,betas(1),betas(2),betas(3),betas(4)
        print*,esubs,asubs
        print*,ts(1),ts(2),ts(3),ts(4)
        print*,rozros
        print*,ibarr
        print*,rcutmeam
        print*,cmin,cmax
        print*,repuls, attrac
        print*,legend

*      Nearest neighbors distance (res)

       print*,'FIND NEAREST NEIGHBORS DISTANCE'
C       ncrys = 1
       print*,'NCRYS =',ncrys 

C       equidist = 4.640
       print*,'equilibrium distance = ',equidist

c      fcc

       if(ncrys.eq.1.or.ncrys.eq.11.or.ncrys.eq.-1) then
          res = equidist/sqrt(2.)
          print*,'RES =',res, '   CRYSTAL = FCC'
c      bcc
       else if(ncrys.eq.2) then
          res = equidist*sqrt(3.)/2.
          print*,'RES =',res, '   CRYSTAL = BCC'
c      dia
       else if(ncrys.eq.0) then
          res = equidist*sqrt(3.)/4
          print*,'RES =',res, '   CRYSTAL = DIAMOND CUBIC'
c      hcp
       else if(ncrys.eq.3) then
          print*,"CAREFUL! WE DON'T HAVE A
     >            SUBROUTINE TO CONSTRUCT THE HCP"
          res = equidist
       else
          print*,'ERROR IN NND'
          call exit(1)
       end if

       return
       end

* Function to calculate the pair potential
* From Baskes' code, dyn87.f, forces.f
c*****************************************************************************
       double precision function phiid(r,it)
       include "moldy.h"
c
c     Identical types
c
*
*     Limit for ZBL
*
        xzbl = -3.
        charge = 94.0 
        re=res
        alpha=alphas
        xmjc=alpha*(r/re-1.d0)
        if(xmjc.le.xzbl) then
           phiid=zbl(r,charge,charge)
        return
        elseif(xmjc.ge.-1.) then
           phizbl=0.
           frac=1.
        else
         frac=xcut((x+1.d0)/(xzbl+1.d0))  !bug!! Wei, 12/15/2006 x->xmjc
         phizbl=zbl(r,charge,charge)
        endif
        beta=betas(1)
        rozero=rozros
        esub=esubs
        repul=repuls
        attra=attrac
        asub=asubs
        ibar=ibarr
        ztmp=zsmeam
        t1=ts(2)
        t2=ts(3)
        t3=ts(4)
        z0meam=zbar(ibar,ztmp,ncrys,t1,t2,t3)
c        z0meam=ztmp
C        write(*,*) 'zsmeam=',zsmeam,'z0meam=',z0meam
        rh0=rhof(r,beta,re,1.d0)
        roz=rh0
c        write(*,*) 'rh0=',rh0
*
* COMMENT THE FOLLOWING FOR NOW SINCE WE WILL START WITH FCC
*
CT         if (lattces(it) .eq. 'dia') then
c
c     Diamond lattice
c
C           beta3=betas(4,it)
C           rh3=rhof(r,beta3,re,1.d0)
C           arg=(32.0/9.0)*ts(4,it)*rh3**2
C           rh0=bar(z*rh0,arg,ibar,z,dum1,dum2)
C           roz=rh0/z0

        if(ncrys.eq.0) then
C          diamond cubic           
           beta3=betas(4)
           rh3=rhof(r,beta3,re,1.d0)
c           write(*,*) 'rh3=',rh3
           arg=(32.0/9.0)*ts(4)*rh3**2
           rh0=bar(z0meam*rh0,arg,ibar,z0meam,dum1,dum2)
           roz=rh0/z0meam
c           write(*,*) 'rho0=',rho0,'rh3=',rh3,'beta3=',
c     C                beta3,'r=',r,'re=',re
CT         elseif (lattces(it) .eq. 'hcp') then
c
c     HCP lattice
c
CT           beta3=betas(4,it)
CT           rh3=rhof(r,beta3,re,1.d0)
CT           arg=(1.0/3.0)*ts(4,it)*rh3**2
CT           rh0=bar(z*rh0,arg,ibar,z,dum1,dum2)
CT           roz=rh0/z0
           elseif (ncrys.eq.3) then
C          hcp
              beta3=betas(4)
              rh3=rhof(r,beta3,re,1.d0)
              arg=(1.0/3.0)*ts(4)*rh3**2
              rh0=bar(z*rh0,arg,ibar,z,dum1,dum2)
              roz=rh0/z0meam
           endif
CT         elseif(lattces(it) .eq. 'dim') then
c
c     dimer
c
CT           beta1=betas(2,it)
CT           beta2=betas(3,it)
CT           beta3=betas(4,it)
CT           rh1=rhof(r,beta1,re,1.d0)
CT           rh2=rhof(r,beta2,re,1.d0)
CT           rh3=rhof(r,beta3,re,1.d0)
CT           t1=ts(2,it)
CT           t2=ts(3,it)
CT           t3=ts(4,it)
CT           arg=t1*rh1**2+2./3.*t2*rh2**2+t3*(1-legend)*rh3**2
CT           rh0=bar(z*rh0,arg,ibar,z,dum1,dum2)
CT           roz=rh0/z0
CT         endif
*
*
c
c     Eqn. (4)
c

 100    phiid=(2.0/ztmp)*(erose(r,re,alpha,esub,repul,attra)
     1    -frhoi(roz,asub,esub))
c        write(*,*) 'phiid=',phiid
        if(r.lt.attra) then
           phiid=phiid+repul*(r-attra)**2
        endif
        phiid=frac*phiid+(1.d0-frac)*phizbl
       return
       end

*
*  Function for ZBL
*
        double precision function zbl(r,z1val,z2val)
        implicit double precision(a-h,o-z)
        dimension c(4),d(4)
        data c/0.028171,0.28022,0.50986,0.18175/
        data d/0.20162,0.40290,0.94229,3.1998/
        data azero/0.4685/
        data cc/14.3997/
        nz1 = nint(z1val)
        nz2 = nint(z2val)
        a=azero/(nz1**0.23+nz2**0.23)
        zbl=0.0
        x=r/a
        do 1 i=1,4
1       zbl=zbl+c(i)*dexp(-d(i)*x)
        zbl=zbl*nz1*nz2/r*cc
        return
        end
c*****************************************************************************
       double precision function xcut(x)
       implicit double precision(a-h,o-z)
         data n,m/2,4/
         xcut=(1.d0-x**m)**n
         return
       end
c*****************************************************************************
       double precision function zbar(ibar,z,lat,t1,t2,t3)
       implicit double precision(a-h,o-z)
         if(ibar.le.0) then
           zbar=z
         else
*
*     FOR THE CASE WE ARE DEALING WITH NOW IBAR = 0
*     LEAVE THIS FOR NOW  
CT           call s(lat,s1,s2,s3)
CT           zbar=bar(z,s1*t1+s2*t2+s3*t3,ibar,z,dum1,dum2)
         print*,'ERROR ibar is not 0'
         call exit(1)
         endif

         return
       end
c*****************************************************************************
       double precision function rhof(r,abc,re,rozero)
       implicit double precision(a-h,o-z)
c
c
c     Eqns. (16a) & (16b)
c
         rhof=rozero*dexp(-abc*(r/re-1.d0))
         return
       end
c*****************************************************************************
       double precision function erose(r,re,alpha,esub,repuls,attrac)
       implicit double precision(a-h,o-z)

        xzbl = -3.0
c
c     Eqns. (14a)
c
        x=alpha*(r/re-1.d0)
        if(x.ge.0) then
        an3=attrac
        elseif (x.lt.0) then
        an3=repuls
        endif
        if(x.ge.-1.d0) then
         erose=-esub*(1.+x+an3*x**3/(r/re))*dexp(-x)
        elseif(x.le.xzbl) then
        erose=0.0
        else
        f=xcut((x+1)/(xzbl+1))
         erose=-f*esub*(1.+x+an3*x**3/(r/re))*dexp(-x)
        endif
        return
        end
c*****************************************************************************
       double precision function frhoi(rhotp,asub,esub)
       implicit double precision(a-h,o-z)
       data const/-100./
c
c     Eqn. (15)
c
        if(rhotp.gt.0.0) then
            frhoi=asub*esub*rhotp*dlog(rhotp)
        else
            frhoi=asub*esub*rhotp*const
        endif
        return
       end
c************************************************************************
c  this routine calculates rhobar
c
       double precision function bar(rho0,A,ibar,z,dang1,dang2)
       implicit double precision(a-h,o-z)
c   rhobar = rho0 * ang
c   derivative of angular term wrt rho0 and A the angular density
         bar=1e-20
         if(rho0.le.0.0) return
         ib=abs(ibar)
         if (ib.eq.0 .or. ib.eq.4) then
c   ang = sqrt(1+A/rho0^2)
           ang = (1+A/rho0**2)
           if(ang.lt.0.0) then
          dang1=0.
          dang2=0.
          return
          else
          ang=sqrt(ang)
              dang1 =  - A/(rho0**3*ang)
              dang2 = 1/(2*rho0**2*ang)
           endif
         elseif (ib.eq.1) then
c  ang = exp(0.5 A/rho0^2)
           ang = dexp(0.5*A/rho0**2)
           dang1 =  -ang*A/rho0**3
           dang2 = ang/(2*rho0**2)
         elseif (ib.eq.2) then
c  ang = exp(0.5 A/z^2)
           ang = dexp(0.5*A/z**2)
           dang1 = 0
           dang2 = ang/(2.0*z**2)
         elseif (ib.eq.3) then
c   ang = 2/(1+exp(-A/rho0^2)
           ex=dexp(-A/rho0**2)
           ang = 2.0/(1.d0+ex)
           dang1 = -ang**2*ex*A/rho0**3
           dang2 = ang**2/2.0*ex/rho0**2
         endif
         bar=rho0*ang
         return
       end
c*****************************************************************************
      double precision function dfrhoi(rho,asub,esub)
      implicit double precision(a-h,o-z)
        data const/-100./
        if(rho.gt.0.0) then
          dfrhoi=asub*esub*(log(rho)+1.)
        else
        dfrhoi=const*asub*esub
        endif
        return
      end
c*****************************************************************************

      subroutine dcmij(drcos,rr,rs)
      implicit double precision(a-h,o-z)
c

c     Find direction cosines

c

        dimension drcos(3,3),rr(3)

        r1 = 1.d0 / rs

        r13 = r1 * r1 * r1

        do 120 i=1,3

          tmp1 = r13 * rr(i)

          do 110 j=1,3

            drcos(i,j) = -tmp1 * rr(j)

 110      continue

          drcos(i,i) = drcos(i,i) + r1

 120    continue

        return

      end
c*****************************************************************************
      double precision function drhof(abc,re,rho)
      implicit double precision(a-h,o-z)
        drhof=-abc*rho/re
        return
      end
c*****************************************************************************
      double precision function rscrn(rij)
c   calculate r cut-off function
        include "moldy.h"
        parameter(frcut = 0.9)
        fac=(rij-frcut*rcutmeam2)/((1.d0-frcut)*rcutmeam2)
        if(fac .le. 0.0) then
          rscrn=1.d0
        elseif(fac .ge. 1.d0) then
          rscrn=0.0
        else
          rscrn=xcut(fac)
        endif
        return
      end
c***************************************************************************
      double precision function dscrn(rij,ril,rjl,cmin,cmax)
      implicit double precision(a-h,o-z)
      parameter(sconst = 1.3)
        cmn=cmin
        cmx=cmax
        dscrn=1.d0
        if ((ril.gt.sconst*rij).or.(rjl.gt.sconst*rij)) goto 200
c
        argb=1.d0-((ril-rjl)/rij)**2
        if(argb .eq. 0.0) go to 200
        argt=4.0*(ril/rij)-((ril-rjl+rij)/rij)**2
        arg=argt/argb
        if(arg .ge. cmx .or. arg .lt. 0.0) then
           goto 200
        elseif (arg .le. cmn) then
           dscrn=0.0
        else
          fac=(cmx-arg)/(cmx-cmn)
          dscrn=xcut(fac)
        endif
 200    continue
        return
      end
c*****************************************************************************
