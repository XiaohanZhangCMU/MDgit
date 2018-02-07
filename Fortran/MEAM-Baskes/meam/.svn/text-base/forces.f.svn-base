c   forces.f
c$Author: baskes $
c$Date: 2004/09/29 19:26:02 $
c$Source: /home/baskes/programs/MEAM/v88/RCS/forces.f,v $
c$Revision: 1.23 $
c$Log: forces.f,v $
cRevision 1.23  2004/09/29 19:26:02  baskes
cfixed B1 + NN
c
cRevision 1.22  2004/02/24 16:10:55  baskes
cadded ialloy
c
cRevision 1.20.1.4  2004/02/24 04:17:32  baskes
cadded ialloy=1 change in t_bar
c
cRevision 1.20.1.3  2003/12/09 21:18:21  baskes
cmerged with dyn96
c
cRevision 1.20.1.2  2003/12/06 04:20:09  baskes
cupdating
c
cRevision 1.20.1.1  2003/10/18 17:20:24  baskes
cfixed zbl
c
cRevision 1.20  2002/04/26 20:55:47  baskes
cchanged averaging
c
cRevision 1.19  2002/04/07 01:21:51  baskes
cfixed screening again
c
cRevision 1.18  2002/03/12 17:17:04  baskes
cfixed screening
c
cRevision 1.17  2002/01/04 23:16:06  baskes
cmodified averaging
c
cRevision 1.16  2002/01/03 23:24:16  baskes
cadded sheardot
c
cRevision 1.15  2001/11/16 00:24:04  baskes
cfixed sign error in dang2 for small rhobar
c
cRevision 1.14  2001/09/17 16:23:14  baskes
cfixed 2NN
cadded ipscree
c
cRevision 1.13  2001/09/14 03:11:28  baskes
cadded beta quartz
c
cRevision 1.12  2001/09/09 14:28:57  baskes
cadded SiO2
c
cRevision 1.11  2001/05/15 21:33:03  baskes
cremoved default repuls
cadded 2NN MEAM
cadded table
c
cRevision 1.10  2001/03/16 22:02:25  baskes
cremoved old attra code
c
cRevision 1.9  2001/01/17 01:02:58  baskes
cfull slocal
c
cRevision 1.8  2001/01/08 23:30:50  baskes
cfixed real seconds
c
cRevision 1.7  2001/01/08 23:00:45  baskes
cfixed header
c
c
c  this routine compute the forces acting on each of the particles
c  as well as the energy of each particle and the contribution of each
c  particle to the stress tensor. Modified 08/05/93.  This version 
c  converges for the ni/si(001) interface. 
c  modified 2/1/94 to include averaged ti
c  modified 6/94 to use polynomial screening cutoff
c  modified 7/30/94 to include cutoff in r at rcut
c  modified 3/29/95 to include arbitrary rhobar
c  modified 8/2/95 to fix dscrn
c  modified 3/1/96 to add noscr
c  modified 3/26/97 to fix ph11 to be 2/z1
c  modified 06/  /97 to streamline distance finding in screening
c                    routines (TPS)
c  modified 07/22/97 to make forces on screening atoms in dscrfor
c                    consistent with their effects on interatomic
c                    potentials in rhosfor and to correct problems
c                    with compilations of slocal in dscrfor (TPS)
c  modified 07/24/97 to restore original ordering of subroutines and 
c                    make two branches: one with nint function calls 
c                    in interatomic distance finding routines, one
c                    without. I find that the code runs faster on 
c                    my machine without calls to nint(). (TPS)
c
c                    This version does not have calls to nint().
c
c	8/8/97 added correction to t3 term (legend)
c	9/1/97 included better negative sqrt for ibar=4
c	3/26/98 added attrac to erose
c	6/16/98 added C11b structure
      subroutine force(time)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      include 'wrkspc.inc'
      common /err/  ibarer,ibare
        common /scre/ noscr
        common /angular/  dang1(natmax),dang2(natmax)
	common /free/ ifree
        dimension rcos(3)
	logical full
	data ifree/0/
	data full/.true./
	data pi/3.1415926/
c
c
c  count the number of times force is called
c  also time the call to force
c
cprint *,'entered force ',ialloy
c        write(*,*) 'force...  natoms = ',natoms, '  rcut=',rcut
c        write(*,*) 'cmin=',cmin(1,1,1)

      ibare=0
      nforce = nforce + 1
      timin = seconds()
c   calculate shear
        if(sheardot.ne.0.0) then
	shear1=(time-t0)*sheardot
        ctbeta=mod(shear1*perlen(3),perlen(1))/perlen(3)
	shear=ctbeta
	beta=pi/2.+atan(shear)
        beta=180./pi*beta
	endif
c
c       clear the arrays used for computing running totals
c
      do 200 i = 1,natoms
        e(i) = 0.
        c(i)=0.
        cg(i)=0.
        do 200 m=1,3
        tav(m,i)=0.
        tsqav(m,i)=0.
        f(m,i) = 0.
          a(m,i)=0.
          ag(m,i)=0.
          do 200 l=1,3
            b(m,l,i)=0.
        slocal(m,l,i) = 0.
            do 200 n=1,3
              d(n,m,l,i)=0.
200     continue
        do 210 k1 = 1,3
          stresst(k1,1) = 0.0
          stresst(k1,2) = 0.0
          stresst(k1,3) = 0.0
210     continue
c
c
c  determine if a new neighbor list is required
c
        call chkdis
c
c   set check variable for updating screening functions
c
c        newlsttst=mod(nforce,nscree)
c        if (nforce.eq.1) newlsttst=1
c
c       calculate total density at each atom
c
      neitot = 0
      neiind(0) = neitot
      do 1000 i=1,natoms
c
c  obtain the information about the neighbors of atom i
c
        call gneigh(i,full)
        nneighs(i) = nneigh
        do 1700 j=1,nneips
          neitot = neitot + 1
          jneighs(neitot) = jneigh(j)
          dneighs(1,neitot) = dneigh(1,j)
          dneighs(2,neitot) = dneigh(2,j)
          dneighs(3,neitot) = dneigh(3,j)
1700    continue
        neiind(i) = neitot  
C        write(*,*) 'neiind(',i,')=',neiind(i),' nneips=',nneips
1000  continue
c
c  Calculate the screening function between all atom pairs if required
c  by a new neighbor list, or by newlsttst
c
c      if ((newlsttst.ne.1).and.(newlst.eq.0)) goto 2100
c
cprint *,'calculating screening'
      call screen
      call dscreen
2100  continue 
c
c       add up pairwise energy and calc the densities
c
      do 3000 i=1,natoms
c
c  obtain the information about the neighbors of the given atom
c

c         if(i.eq.11) then
c            return
c         endif

        nneips = neiind(i) - neiind(i-1)
        nneigh = nneighs(i)
        do 3100 j=1,nneips
          jj = neiind(i-1) + j
          jneigh(j) = jneighs(jj)
          dneigh(1,j) = dneighs(1,jj)
          dneigh(2,j) = dneighs(2,jj)
          dneigh(3,j) = dneighs(3,jj)
          rneigh(j) = sqrt(dneigh(1,j)**2 + dneigh(2,j)**2 +
     1                 dneigh(3,j)**2)
3100    continue
c
c     compute the contribution to rho(i) from particle j
c
c     compute the contribution to rho(j) from particle i
c
        ity = itype(i)
        do 1300 j=1,nneips
          jty=itype(jneigh(j))
          scree=scrnab(i,j)
	if (i.eq.ipscree .or. ipscree.gt.natoms ) then
	print *,'i,j,screen= ',i,jneigh(j),scree
	endif
          if (scree .lt. scnres) goto 1300
c
c     calc phif between pairs and store in e
c
cprint *,'calling phi ',i,j
          phifs=phif(rneigh(j),ity,jty)

c          if((i.le.10)) then
c             write(*,*) 'i=',i,'phifs=',phifs,' scree=',scree
c             return
c          endif

          e(i)=e(i)+0.5*phifs*scree

c      do i=1,natoms
c          if(i.eq.1) then
c             write(*,*) 'e(',i,')=',e(i),phifs,scree,rneigh(j),ity,jty
c          endif
c      enddo


c      e(jneigh(j))=e(jneigh(j))+0.5*phifs*scree
          fp(i)=e(i)
c      fp(jneigh(j))=e(jneigh(j))
c
c     compute each partial electron density at atom i due to each
c     neighbor atom j and visa versa
        rej=res(jty,jty)
        rzj=rozros(jty)
        betaj0=betas(1,jty)
        betaj1=betas(2,jty)
        betaj2=betas(3,jty)
        betaj3=betas(4,jty)
	irho=irhos(jty)
	if(ialloy.eq.0) then
	fac1=1.
	fac2=1.
	fac3=1.
	else
	fac1=ts(2,jty)
	fac2=ts(3,jty)
	fac3=ts(4,jty)
	endif
	nrho=0
        rhj=rhof(rneigh(j),betaj0,rej,rzj)
	nrho=2
        rogj1=rhof(rneigh(j),betaj1,rej,rzj*fac1)
	nrho=4
        rogj2=rhof(rneigh(j),betaj2,rej,rzj*fac2)
	nrho=6
        rogj3=rhof(rneigh(j),betaj3,rej,rzj*fac3)
c
c     Eqns. (8a),(8b),(8c),(8d)
c
c   the minus signs on odd  term due to diff in direc. cosines
c   what is meant here is that for x and xxx a minus sign is needed
c   while for xx a plus sign is ok
c
        r1 = 1.0 / rneigh(j)
        do 25 m=1,3
          rcos(m) = r1 * dneigh(m,j)
 25     continue
        do 50 m=1,3
          t1=scree*rcos(m)
          a(m,i)=a(m,i)-rogj1*t1
          ag(m,i)=ag(m,i)-rogj3*t1
          do 40 l=1,3
            t2=t1*rcos(l)
            b(m,l,i)=b(m,l,i)+rogj2*t2
            do 30 n=1,3
              d(m,l,n,i)=d(m,l,n,i)-rogj3*t2*rcos(n)
 30         continue
 40       continue
 50     continue
c
c
c
        c(i)=c(i)+rhj*scree
        tav(1,i)=tav(1,i)+ts(2,jty)*rhj*scree
        tav(2,i)=tav(2,i)+ts(3,jty)*rhj*scree
        tav(3,i)=tav(3,i)+ts(4,jty)*rhj*scree
	if(ialloy.ne.0) then
        tsqav(1,i)=tsqav(1,i)+ts(2,jty)**2*rhj*scree
        tsqav(2,i)=tsqav(2,i)+ts(3,jty)**2*rhj*scree
        tsqav(3,i)=tsqav(3,i)+ts(4,jty)**2*rhj*scree
	endif
        cg(i)=cg(i)+rogj2*scree
1300    continue
3000  continue   
cprint *,'after 3000 ',tav(2,108),tsqav(2,108)
c
c     now find f(rho) for each particle
c
      do 2000 i=1,natoms
        rho0=c(i)
        if(rho0.gt.0) then
          ity=itype(i)
c
c     define weighting factors etc. for atomic densities
c
c     t1=ts(2,ity)
c     t2=ts(3,ity)
c     t3=ts(4,ity)
c
c          write(*,*) 'c(',i,')=',c(i)

        if(ialloy.eq.0) then
          t1=tav(1,i)/c(i)
          t2=tav(2,i)/c(i)
          t3=tav(3,i)/c(i)
        else
          t1=0
          t2=0
          t3=0
          if(tsqav(1,i).gt.0) t1=tav(1,i)/tsqav(1,i)
          if(tsqav(2,i).gt.0) t2=tav(2,i)/tsqav(2,i)
          if(tsqav(3,i).gt.0) t3=tav(3,i)/tsqav(3,i)
        endif
          ibar=ibarr(ity)
          asub=asubs(ity)
          esub=esubs(ity,ity)
          rozeroi=rozros(ity)
          z=zs(ity)
c
c      z0=z0s(ity)
c
cprint *,'before bar ',t1
          z0=zbar(ibar,z,lattces(ity),t1,t2,t3,nn,z0s(ity))
c         z0=zbar(ibar,z,lattces(ity),t1,t2,t3)
c
c     Eqns (8a),(8b),(8c) & (8d)
c
          rho1=0.
          rho2=-(cg(i)**2)/3.0
          rho3=0.
	do 1399 m=1,3
1399	rho3=rho3-legend*ag(m,i)**2
          do 1400 m=1,3
            rho1=rho1+a(m,i)**2
            do 1400 l=1,3
              rho2=rho2+b(m,l,i)**2
              do 1400 n=1,3
                rho3=rho3+d(n,m,l,i)**2
1400      continue
c
c     ~ Eqn. (9a)
c
          rhsq(1,i)=rho1
          rhsq(2,i)=rho2
          rhsq(3,i)=rho3
          drh=(t1*rho1+t2*rho2+t3*rho3)
cif(i.eq.1) then
cprint *,'before bar ',t1,t2,t3,rho1,rho2,rho3
cprint *,'before bar ',i,t2,rho0,drh
cendif
          rhobr=bar(rho0,drh,ibar,z*rozeroi,dang1(i),dang2(i))

c debug Wei Cai, 11/28/2006 (for Ni, rhobar=rho0 but for Si they are not the same)
c        write(*,*) i, rho0, rhobr, rho1, rho2, rho3, drh


c
c     note the division of rho by z0*rozeroi here & below......
c
          rhobar=rhobr/(z0*rozeroi)
          rho(i)=rhobr
c          write(*,*) 't1=',t1,' rho1=',rho1
c          write(*,*) 't2=',t2,' rho2=',rho2
c          write(*,*) 't3=',t3,' rho3=',rho3
c          write(*,*) 'rho(',i,')=',rho(i),' drh=',drh

c	if(i.eq.1) then
cprint *,' after 2000 ',rho(1),plocal(1),rhobar,rozeroi,rhobr
cprint *,' after 2000 ',rho(1)
c	endif
          plocal(i)=frhoi(rhobar,asub,esub)
          e(i) = e(i) + plocal(i)
	if (rho(i).lt.0.) then
	print *,'**********warning**********'
	print *, 'for atom ',i,' rhobar= ',rhobar,plocal(i)
	endif
        else
          plocal(i)=0
          rho(i)=1.e-20
          rhsq(1,i)=0
          rhsq(2,i)=0
          rhsq(3,i)=0
        endif
2000  continue
cprint *,' after 2000 ',rhobr
c
c     compute forces:
c
c
c     forces due to pair-wise potential
c
      call pwsfor
c
c     forces due to embedding energies
c    
      call rhosfor
c
c     forces due to dscreen terms, i.e. 3 body
c
      if(noscr.eq.0) call dscrfor
      if(ifree.ne.0) return
c
c  add the constraints to the forces
c
cprint *,' before fixfor ',f(1,1)
      call fixfor(time)
cprint *,' after fixfor ',f(1,1)
      call modfor
c
c  compute the stress tensor (commented out by Wei Cai, 12/4/2006, only have potential energy part)
c                            (kinetic energy part in MD++)
c
c      do 3050 i=1,natoms
c        stresst(1,1) = stresst(1,1)
c     1             + amass(itype(i))*rv(4,i)*rv(4,i)
c        stresst(1,2) = stresst(1,2)
c     1             + amass(itype(i))*rv(4,i)*rv(5,i)
c        stresst(1,3) = stresst(1,3)
c     1             + amass(itype(i))*rv(4,i)*rv(6,i)
c        stresst(2,1) = stresst(2,1)
c     1             + amass(itype(i))*rv(4,i)*rv(5,i)
c        stresst(2,2) = stresst(2,2)
c     1             + amass(itype(i))*rv(5,i)*rv(5,i)
c        stresst(2,3) = stresst(2,3)
c     1             + amass(itype(i))*rv(5,i)*rv(6,i)
c        stresst(3,1) = stresst(3,1)
c     1             + amass(itype(i))*rv(4,i)*rv(6,i)
c        stresst(3,2) = stresst(3,2)
c     1             + amass(itype(i))*rv(5,i)*rv(6,i)
c        stresst(3,3) = stresst(3,3)
c     1             + amass(itype(i))*rv(6,i)*rv(6,i)
c
c  the 1/2 is factored in during loop 4000
c
3050  continue
c
c  sum up the local contributions to the stress tensor to get the total stress
c  tensor and compute the other half of the symmetric matrix
c  for MEAM tensor is not necessarily symmetric
c
      sum1=0.0
      sum2=0.0
      sum3=0.0
      do 4000 i = 1,natoms
        sum1=sum1+f(1,i)
        sum2=sum2+f(2,i)
        sum3=sum3+f(3,i)
        slocal(1,1,i) = 0.5*slocal(1,1,i)
        slocal(1,2,i) = 0.5*slocal(1,2,i)
        slocal(1,3,i) = 0.5*slocal(1,3,i)
        slocal(2,1,i) = 0.5*slocal(2,1,i)
        slocal(2,2,i) = 0.5*slocal(2,2,i)
        slocal(2,3,i) = 0.5*slocal(2,3,i)
        slocal(3,1,i) = 0.5*slocal(3,1,i)
        slocal(3,2,i) = 0.5*slocal(3,2,i)
        slocal(3,3,i) = 0.5*slocal(3,3,i)
      avslocal(1,1,i) = avslocal(1,1,i)+slocal(1,1,i)
      avslocal(1,2,i) = avslocal(1,2,i)+slocal(1,2,i)
      avslocal(1,3,i) = avslocal(1,3,i)+slocal(1,3,i)
      avslocal(2,1,i) = avslocal(2,1,i)+slocal(2,1,i)
      avslocal(2,2,i) = avslocal(2,2,i)+slocal(2,2,i)
      avslocal(2,3,i) = avslocal(2,3,i)+slocal(2,3,i)
      avslocal(3,1,i) = avslocal(3,1,i)+slocal(3,1,i)
      avslocal(3,2,i) = avslocal(3,2,i)+slocal(3,2,i)
      avslocal(3,3,i) = avslocal(3,3,i)+slocal(3,3,i)
        stresst(1,1) = stresst(1,1) + slocal(1,1,i)
        stresst(1,2) = stresst(1,2) + slocal(1,2,i)
        stresst(1,3) = stresst(1,3) + slocal(1,3,i)
        stresst(2,1) = stresst(2,1) + slocal(2,1,i)
        stresst(2,2) = stresst(2,2) + slocal(2,2,i)
        stresst(2,3) = stresst(2,3) + slocal(2,3,i)
        stresst(3,1) = stresst(3,1) + slocal(3,1,i)
        stresst(3,2) = stresst(3,2) + slocal(3,2,i)
        stresst(3,3) = stresst(3,3) + slocal(3,3,i)
 4000 continue
      call fixtmp(time)
c
c  add the pv value to the end of the e array when ibdtyp = 2
c
      if (ibdtyp.eq.2)
     1         e(natoms+1) = dpress*perlen(1)*perlen(2)*perlen(3)
     2                     + 0.5*( dstress(1)*perlen(1)**2
     3                           + dstress(2)*perlen(2)**2
     4                           + dstress(3)*perlen(3)**2)   
c
c  determine the timing for this call
c
9001  format(5g15.5)
      frctim = seconds() - timin
      frctmx = max(frctmx,frctim)
      frctmn = min(frctmn,frctim)
cprint *,' leaving force ',e(1)

      return
      end
c************************************************************************
c  this routine calculates rhobar
c
      function bar(rho0,A,ibar,z,dang1,dang2)
        implicit real*8 (a-h,o-z)
	data x0,n/-0.99,99/
c   n=-x_0/(1+x_0)
c   rhobar = rho0 * ang
c   derivative of angular term wrt rho0 and A the angular density
        common /err/  ibarer,ibare
        data ibarer/0/
        bar=1e-20
        if(rho0.le.0.0) return
        ib=abs(ibar)
        if (ib.eq.0 .or. ib.eq.4) then
c   ang = sqrt(1+A/rho0^2)
          ang = (1+A/rho0**2)
          if(ang.lt.1+x0) then
	  x=ang-1
	  ang=sqrt((1+x0)*(x0/x)**n)
	  dang1=ang*n/rho0
	  dang2=-ang*n/(2*A)
             ibare=ibare+1
             ibarer=ibarer+1
cprint *,'in bar ',rho0,A,x0
	  else
	  ang=sqrt(ang)
             dang1 =  - A/(rho0**3*ang)
             dang2 = 1/(2*rho0**2*ang)
          endif
        elseif (ib.eq.1) then
c  ang = exp(0.5 A/rho0^2)
          ang = exp(0.5*A/rho0**2)
          dang1 =  -ang*A/rho0**3
          dang2 = ang/(2*rho0**2)
        elseif (ib.eq.2) then
c  ang = exp(0.5 A/z^2)
          ang = exp(0.5*A/z**2)
          dang1 = 0
          dang2 = ang/(2.0*z**2)
        elseif (ib.eq.3) then
c   ang = 2/(1+exp(-A/rho0^2)
          ex=exp(-A/rho0**2)
          ang = 2.0/(1.0+ex)
          dang1 = -ang**2*ex*A/rho0**3
          dang2 = ang**2/2.0*ex/rho0**2
        endif
        bar=rho0*ang
        return
      end
c**********************************************************************
      subroutine screen
c
c     Find the screening function for all interactions
c
c     screen 07/12/93
c
c       Calculate screening between atom i and atom j resulting from all
c       other neighbors (l) of i in the vicinity of i and j
c
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      include 'wrkspc.inc'
        common /scre/ noscr
c
c     Loop through all atoms i
c
c        icount = 0  ! number of times to call rscrn,dscrn
      do 100 i=1,natoms
        ity=itype(i)
c
c     Obtain the information about the neighbors of the given atom
c
        nlow = neiind(i-1)
        nhigh = neiind(i)
c          if(i.eq.1) then
c             write(*,*) 'atom(',i,') nn=',nhigh-nlow
c          endif
c
c     Count through all neighbors, j
c
        do 110 j=nlow+1,nhigh
          jindx = jneighs(j)
          jty = itype(jindx)
c
c     Find distance between i and j
c
          dxij = -dneighs(1,j)
          dyij = -dneighs(2,j)
          dzij = -dneighs(3,j)
          rij = dxij*dxij + dyij*dyij + dzij*dzij
c
c     Find two-body screening -- if complete, forgo next loop
c
          screenij = rscrn(rij)
          if (noscr .eq. 1 .or. screenij .le. scnres) goto 900
c
c     Loop through mutual neighbors
c
          do 200 l=nlow+1,nhigh
            if (l.eq.j) goto 200
            
c            icount = icount + 1
            lindx = jneighs(l)
            lty=itype(lindx)
c
c     Find the separations between i and l, and j and l
c     Include the periodic effects
c
            dxil = -dneighs(1,l)
            dyil = -dneighs(2,l)
            dzil = -dneighs(3,l)
            ril = dxil*dxil + dyil*dyil + dzil*dzil
            rjl = rij + ril 
     1            - 2.0 * (dxij*dxil + dyij*dyil + dzij*dzil)
c
            screenij = screenij * dscrn(rij,ril,rjl,ity,jty,lty)
c            icount = icount + 1
c     
c     This exit condition is slightly different from the original, but
c     I think that it gives the intended result. (TPS, 6/11/97)
c
            if (screenij .le. scnres) goto 900
200       continue
c
c     Set the screening between the pair
c
 900      scrnab(i,j-nlow)=screenij
          if(screenij.gt.scnres) rscrmax=max(rscrmax,sqrt(rij))
c
c     Move on the the next pair
c
 110    continue
 100  continue

c      write(*,*) 'screen: icount=',icount,' scnres=',scnres
      return
      end
c********************************************************************
      subroutine dscreen
c
c     Find the derivative of the screening function for all interactions
c
      include 'param.inc'
        include 'dyn.inc'
        include 'meam.inc'
        include 'wrkspc.inc'
        common /scre/ noscr
        dimension scrp(3),scrm(3)
c
c     Increment for differentiation
c
c        icount = 0
        h=1.0e-5
        hh = h * h
        h2 = 2.0 * h
        do 100 i=1,natoms
          ity = itype(i)
          nlow = neiind(i-1)
          nhigh = neiind(i)
c
c    Work with j referencing directly into *neighs(*) arrays
c
          do 110 j=nlow+1,nhigh
c
c    See if derivative of screening term is different from 0
c
            scree=scrnab(i,j-nlow)
c
c    If not, set derivatives to zero and exit from the neighbor loop
c
            if (scree .gt. 1.0-scnres .or. scree .lt. scnres) then
              do 105 m=1,3
                dscrnab(m,i,j-nlow)=0.0 
 105          continue
              goto 110
            endif
c
c    Two-body screening:
c
            jty = itype(jneighs(j))
            dxij = -dneighs(1,j)
            dyij = -dneighs(2,j)
            dzij = -dneighs(3,j)
            rij = dxij*dxij + dyij*dyij + dzij*dzij
c
c    Perform numerical differentiation
c
            do 1020 m=1,3
              dij = -dneighs(m, j)
c
c    positive increment in xi's
c
              scrp(m)=rscrn (rij - h2 * dij + hh)
c
c    negative increment in xi's
c
              scrm(m)=rscrn (rij + h2 * dij + hh)
c              icount = icount + 2
c            if(i.eq.1) then
c               write(*,*) 'i=',i,'jindx=',jneighs(j)
c               write(*,*) 'r2ij_p(',m,')=',rij-h2*dij+hh
c               write(*,*) 'r2ij_m(',m,')=',rij+h2*dij+hh
c               write(*,*) 'scrp(',m,')=',scrp(m)
c               write(*,*) 'scrm(',m,')=',scrm(m)
c            endif
 1020       continue


            if(noscr.eq.1) go to 121
            do 120 l=nlow+1,nhigh
              if (l.eq.j) goto 120
c
c    Three body screening:
c
              lty = itype(jneighs(l))
              dxil = -dneighs(1,l)
              dyil = -dneighs(2,l)
              dzil = -dneighs(3,l)
              ril = dxil*dxil + dyil*dyil + dzil*dzil
              rjl = rij + ril -
     .              2.0 * (dxij*dxil + dyij*dyil + dzij*dzil)
c
c    Incrementing ith atom position
c
              do 102 m=1,3
c
c    Perform numerical differentiation
c
                dil = -dneighs(m, l)
                dij = -dneighs(m, j)
c
c    positive increment in xi's
c
                rpij = rij - h2 * dij + hh
                rpil = ril - h2 * dil + hh
                scrp(m)=scrp(m)*dscrn(rpij,rpil,rjl,ity,jty,lty)
c
c    negative increment in xi's
c
                rpij = rij + h2 * dij + hh
                rpil = ril + h2 * dil + hh
                scrm(m)=scrm(m)*dscrn(rpij,rpil,rjl,ity,jty,lty)
c                icount = icount + 2
 102          continue
 120        continue
 121        continue
c
c    find and store derivative
c
            do 107 m=1,3
              dscrnab(m,i,j-nlow)=0.5*(scrp(m)-scrm(m))/h
 107        continue

c            if((i.eq.1)) then
c               write(*,*) 'i=',i,'jindx=',jneighs(j)
c               write(*,*) 'dscrnab = ',dscrnab(1,i,j-nlow)
c               write(*,*) 'dscrnab = ',dscrnab(2,i,j-nlow)
c               write(*,*) 'dscrnab = ',dscrnab(3,i,j-nlow)
c            endif
            
 110      continue
 100    continue
c        write (*,*) 'dscreen: icount=',icount
        return
      end
c********************************************************************
      function rscrn(rij)
c   calculate r cut-off function
      include 'param.inc'
        include 'dyn.inc'
        include 'meam.inc'
	common /fr/ frcut,xncut,xmcut
        fac=(rij-frcut*rcutsq)/((1.0-frcut)*rcutsq)
        if(fac .le. 0.0) then
          rscrn=1.0
        elseif(fac .ge. 1.0) then
          rscrn=0.0
        else
          rscrn=xcut(fac)
        endif
        return
      end
c***************************************************************************
      function dscrn (rij,ril,rjl,ityp,jtyp,ltyp)
      include 'param.inc'
        include 'dyn.inc'
        include 'meam.inc'
        cmn=cmin(ityp,ltyp,jtyp)
        cmx=cmax(ityp,ltyp,jtyp)
        dscrn=1.0
        if ((ril.gt.sconst*rij).or.(rjl.gt.sconst*rij).or.cmx.eq.0.0)
     1	 goto 200
c
c    I don't understand this functionality, and I wish this code
c    was better documented. (TPS, 6/97)
c
c   If l atom is outside of ellipse there is no screening
c   dotx = 2 r_ij dot r_x
           dotil=ril+rij-rjl
        if(dotil.le.2.e-3) go to 200
           dotjl=rij+rjl-ril
        if(dotjl.le.2.e-3) go to 200
c       argb=1.0-((ril-rjl)/rij)**2
c       if(abs(argb) .le. 1.0e-6) go to 200
c       argt=4.0*(ril/rij)-((ril-rjl+rij)/rij)**2
        argb=rij**2-(ril-rjl)**2
        argt=4.0*(ril*rij)-(ril-rjl+rij)**2
        arg=argt/argb
        if(arg .ge. cmx ) then
           dscrn=1.0
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
      function xcut(x)
        implicit real*8 (a-h,o-z)
	common /fr/ frcut,xncut,xmcut
c      data n,m/2,4/
c      data n,m/2,2/
        xcut=(1.0-x**xmcut)**xncut
        return
      end
c*****************************************************************************
      function rhof(r,abc,re,rozero)
        implicit real*8 (a-h,o-z)
        common /rhon/ irho,nrho
c
c
c     Eqns. (16a) & (16b)
c
	if(irho.eq.0) then
        rhof=rozero*exp(-abc*(r/re-1.0))
	else
        rhof=rozero*(r/re)**nrho*exp(-abc*(r/re-1.0))
	endif
        return
      end
c*****************************************************************************
      function drhof(r,abc,re,rho)
        implicit real*8 (a-h,o-z)
        common /rhon/ irho,nrho
	if(irho.eq.0) then
        drhof=-abc*rho/re
	else
        drhof=(nrho/r-abc/re)*rho
	endif
        return
      end
c*****************************************************************************
      function erose(r,re,alpha,esub,repuls,attrac)
        implicit real*8 (a-h,o-z)
ccommon /czbl/ xzbl
c
c     Eqns. (14a)
c
        x=alpha*(r/re-1.0)
	if(x.ge.0) then
	an3=attrac
	elseif (x.lt.0) then
	an3=repuls
	endif
cif(x.ge.-1.0) then
        erose=-esub*(1.+x+an3*x**3/(r/re))*exp(-x)
celseif(x.le.xzbl) then
cerose=0.0
celse
cf=xcut((x+1)/(xzbl+1))
c       erose=-f*esub*(1.+x+an3*x**3/(r/re))*exp(-x)
cendif
c        write(*,*) 'x=',x,' erose=',erose
	return
	end
c*****************************************************************************
      function frhoi(rho,asub,esub)
        implicit real*8 (a-h,o-z)
        common /err/  ibarer,ibare
	data const/-100./
c
c     Eqn. (15)
c
        if(rho.gt.0.0) then
           frhoi=asub*esub*rho*log(rho)
	else
        frhoi=asub*esub*rho*const
	endif	
        return
      end
c*****************************************************************************
      function dfrhoi(rho,asub,esub)
        implicit real*8 (a-h,o-z)
        common /err/  ibarer,ibare
	data const/-100./
        if(rho.gt.0.0) then
          dfrhoi=asub*esub*(log(rho)+1.)
	else
        dfrhoi=const*asub*esub
        endif
	return
      end
c*****************************************************************************
      subroutine calphiid(r,it,phiid_var)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      
      phiid_var = phiid(r,it)
      return
      end

      function phiid(r,it)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
	common /czbl/ xzbl,xzbl0
c
c     Identical types
c
cprint *,'in phi ',r
        re=res(it,it)
        alpha=alphas(it,it)
          beta0=betas(1,it)
          beta1=betas(2,it)
          beta2=betas(3,it)
          beta3=betas(4,it)
        rozero=rozros(it)
        esub=esubs(it,it)
        repul=repuls(it,it)
        attra=attrac(it,it)
        asub=asubs(it)
        ibar=ibarr(it)
        irho=irhos(it)
        z=zs(it)
        t1=ts(2,it)
        t2=ts(3,it)
        t3=ts(4,it)
	z0=zbar(ibar,z,lattces(it),t1,t2,t3,nn,z0s(it))
        
c        write(*,*) 'z=',z,' z0=',z
c       z0=zbar(ibar,z,lattces(it),t1,t2,t3)
        nrho=0
        rh0=rhof(r,beta0,re,1.0d0)
        nrho=2
        rh1=rhof(r,beta1,re,1.0d0)
        nrho=4
        rh2=rhof(r,beta2,re,1.0d0)
        nrho=6
        rh3=rhof(r,beta3,re,1.0d0)
        x=alpha*(r/re-1.0)

c        write(*,*) 'r=',r,' beta0=',beta0,' re=',re,' rh0=',rh0
c        write(*,*) 'phiid: ',re,x,alpha,r,re,alats(1),res(1,1),it
c        write(*,*) 'rh0=',rh0,' rh1=',rh1,' rh2=',rh2,' rh3=',rh3
	if(x.le.xzbl0) then
c
c     dimer
c
          arg=t1*rh1**2+2./3.*t2*rh2**2+t3*(1-legend)*rh3**2
	zdimer=1
cprint *,'at ZBL ',x,xzbl0,arg
          roz=bar(zdimer*rh0,arg,ibar,zdimer,dum1,dum2)/z0

	fdimer=frhoi(roz,asub,esub)
	phizbl=zbl(r,ielement(it),ielement(it))
	   if(x.lt.xzbl) then
	   frac=0.0
	   else
           frac=xcut((x-xzbl0)/(xzbl-xzbl0))
	   endif
	else
	frac=1.
	phizbl=0.
	fdimer=0
	endif
	if(nn) then
	roz=z*(rh0-b2nn(it)*rhof(r*a2nn(it),beta0,re,1.0d0))
     1       /z0s(it)
	else
        roz=rh0
	endif


        if (lattces(it) .eq. 'dia') then
c
c     Diamond lattice
c
          arg=(32.0/9.0)*t3*rh3**2
          roz=bar(z*rh0,arg,ibar,z,dum1,dum2)/z0
          
c          write(*,*) 'arg=',arg,' rh0=',rh0,' ibar=',ibar,' z=',z
        elseif (lattces(it) .eq. 'hcp') then
c
c     HCP lattice
c
          arg=(1.0/3.0)*t3*rh3**2
c         rh0=bar(z*rh0,arg,ibar,z,dum1,dum2)
c         roz=bar(z*rh0,arg,ibar,z,dum1,dum2)/z0
c       if(nn) then
croz=zs(it)*(rh0-b2nn(it)*rhof(r*a2nn(it),beta1,re,1.0d0))
c    1       /z0s(it)
c       else
c       roz=bar(z*rh0,arg,ibar,z,dum1,dum2)
c       endif
          roz=bar(z*rh0,arg,ibar,z,dum1,dum2)/z0
        elseif(lattces(it) .eq. 'dim') then
c
c     dimer
c
          arg=t1*rh1**2+2./3.*t2*rh2**2+t3*(1-legend)*rh3**2
          roz=bar(z*rh0,arg,ibar,z,dum1,dum2)/z0
        endif
c
c     Eqn. (4)
c

crose=erose(r,re,alpha,esub,repul,attra)
cembed=frhoi(roz,asub,esub)
cprint *,'in phiid ',r,rh0,roz,rose,embed
100    continue
	fff=frhoi(roz,asub,esub)
	phiid=(2.0/z)*(erose(r,re,alpha,esub,repul,attra)-fff)
c        write(*,*) 'roz=',roz,' asub=',asub,' esub=',esub
c        write(*,*) 'erose=',erose(r,re,alpha,esub,repul,attra),
c     1       ' frhoi=',fff,' z=',z
c       write(*,*) 'r=',r,' frac=',frac
c       write(*,*) ' phiid=',phiid,' phizbl=',phizbl
c       write(*,*) '**',(2.0/z)*(erose(r,re,alpha,esub,repul,attra)-fff)

       phiid=(2.0/z)*(erose(r,re,alpha,esub,repul,attra)-fff)
c       write(*,*) '***',phiid

c       fdimer=0.d0 !corresponding to LAMMPS implementation
       phiid=frac*phiid+(1-frac)*(phizbl-2*fdimer)
c        write(*,*) 'phiid=',phiid,' ******* fdimer=',fdimer

c        write(*,*) 'phiid(',r,it,')=',phiid
c        write(*,*) 'frac=',frac,'phizbl=',phizbl,'fdimer=',fdimer
        
      return 
      end

      function phif(r,it,jt)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
	common /czbl/ xzbl,xzbl0
      common /alloy/ all(nelmax,nelmax)
        common /scre/ noscr
      character *4 all
	if((it.eq.jt.and.itable.lt.0).or.
     1        (it.ne.jt.and.itable.lt.-1))  then
	 p=r*ntable/rcut
	 k=p
	  if(k.gt.ntable-1) then
	  phif=0.
	  return
	  endif
	 p=p-k
	 phif=(1-p)*phitab(k,it,jt)+p*phitab(k+1,it,jt)
	 return
	endif

      if (it .eq. jt) then
c
c   Identical types
c
	phif=phiid(r,it)

	return
c
c   Mixed types
c
      else
        re12=res(it,jt)
        alpha12=alphas(it,jt)
        i1=min(it,jt)
        i2=max(it,jt)
        ibar1=ibarr(i1)
        ibar2=ibarr(i2)
        beta1=betas(1,i1)
        beta2=betas(1,i2)
        re1=res(i1,i1)
        re2=res(i2,i2)
        rozro1=rozros(i1)
        rozro2=rozros(i2)
        alpha1=alphas(i1,i1)
        alpha2=alphas(i2,i2)
        esub1=esubs(i1,i1)
        esub2=esubs(i2,i2)
        esub12=esubs(i1,i2)
	repul1=repuls(i1,i1)
	repul2=repuls(i2,i2)
	repul12=repuls(i1,i2)
	attra1=attrac(i1,i1)
	attra2=attrac(i2,i2)
	attra12=attrac(i1,i2)
        asub1=asubs(i1)
        asub2=asubs(i2)
        z1=zs(i1)
        z2=zs(i2)
c
	x=alpha12*(r/re12-1.0)
        if(x.le.xzbl0) then
c
c     Cases for all opposite neighbors
c
          nrho=0
          irho=irhos(i1)
          rh1=rhof(r,beta1,re1,rozro1)
          irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
          z01=zbar(ibar1,z1,lattces(i1),ts(2,i2),ts(3,i2),ts(4,i2),
     1        nn,z0s(i1))
          z02=zbar(ibar2,z2,lattces(i2),ts(2,i1),ts(3,i1),ts(4,i1),
     1        nn,z0s(i2))
c
c     dimer
c
          zij=1.
          irho=irhos(i1)
          nrho=2
          rh11=rhof(r,betas(2,i1),re1,rozro1)
          nrho=4
          rh21=rhof(r,betas(3,i1),re1,rozro1)
          nrho=6
          rh31=rhof(r,betas(4,i1),re1,rozro1)
          t1=ts(2,i1)
          t2=ts(3,i1)
          t3=ts(4,i1)
          arg=t1*rh11**2+2./3.*t2*rh21**2+t3*(1-legend)*rh31**2
          rh1=bar(rh1,arg,ibar2,zij*rozro1,dum1,dum2)
          irho=irhos(i2)
          nrho=2
          rh11=rhof(r,betas(2,i2),re2,rozro2)
          nrho=4
          rh21=rhof(r,betas(3,i2),re2,rozro2)
          nrho=6
          rh31=rhof(r,betas(4,i2),re2,rozro2)
          t1=ts(2,i2)
          t2=ts(3,i2)
          t3=ts(4,i2)
          arg=t1*rh11**2+2./3.*t2*rh21**2+t3*(1-legend)*rh31**2
          rh2=bar(rh2,arg,ibar1,zij*rozro2,dum1,dum2)
          rh1z=rh1/(z02*rozro2)
          rh2z=rh2/(z01*rozro1)
 	  fadimer=frhoi(rh2z,asub1,esub1)
	  fbdimer=frhoi(rh1z,asub2,esub2)
	  phizbl=zbl(r,ielement(it),ielement(jt))
	   if(x.lt.xzbl) then
	    frac=0.0
	   else
            frac=xcut((x-xzbl0)/(xzbl-xzbl0))
	   endif
        else
          phizbl=0.0
	  fadimer=0.0
	  fbdimer=0.0
	  frac=1.0
        endif
c
c find the reference structure
        if(all(it,jt).eq.'L12') then
          zij=12.0
	  nrho=0
	  irho=irhos(i1)
          rh1=rhof(r,beta1,re1,rozro1)
	  irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
	  z02=zbar(ibar2,z2,lattces(i2),ts(2,i1),ts(3,i1),ts(4,i1),
     1       nn,z0s(i2))
c         z02=zbar(ibar2,z2,lattces(i2),ts(2,i1),ts(3,i1),ts(4,i1))
          rh1z=zij*rh1/(z02*rozro2)
c
c   rhobarme is not just the
c   sum of zeroth components of the electron densities.  The
c   second component enters for the L12 structure
c
          rh20=(4.*rh2+8.*rh1)
          beta21=betas(3,i1)
          beta22=betas(3,i2)
          nrho=4
          irho=irhos(i1)
          rh21=rhof(r,beta21,re1,rozro1)
          irho=irhos(i2)
          rh22=rhof(r,beta22,re2,rozro2)
          t1av=(4*rh2*ts(2,i2)+8*rh1*ts(2,i1))/rh20
          t2av=(4*rh2*ts(3,i2)+8*rh1*ts(3,i1))/rh20
          t3av=(4*rh2*ts(4,i2)+8*rh1*ts(4,i1))/rh20
	  z01=zbar(ibar1,z1,lattces(i1),t1av,t2av,t3av,nn,z0s(i1))
c         z01=zbar(ibar1,z1,lattces(i1),t1av,t2av,t3av)
	  if(ialloy.eq.0) then
           rho22=(8./3.)*(rh21-rh22)**2
	  else
	   denom=8*rh1*ts(3,i1)**2+4*rh2*ts(3,i2)**2
	   rho22=(8./3.)*(rh21*ts(3,i1)-rh22*ts(3,i2))**2
	   if(denom.gt.0.) rho22=rho22/denom*rh20
	  endif
          rhobarme=bar(rh20,t2av*rho22,ibar1,z01*rozro1,dum1,dum2)
          rh2z=(rhobarme)/(z01*rozro1)
	  ph11=phiid(r,i1)
c
c       Eqn. (23)
c 
          fa=frhoi(rh2z,asub1,esub1)
          fb=frhoi(rh1z,asub2,esub2)
          phif=erose(r,re12,alpha12,esub12,repul12,attra12)/3.0-
     1           fa/4.0 -fb/12.0 -ph11
        elseif(all(it,jt).eq.'L10') then
          zij=12.0
          nrho=0
          irho=irhos(i1)
          rh1=rhof(r,beta1,re1,rozro1)
          irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
c
c   rhobarme is not just the
c   sum of zeroth components of the electron densities.  The
c   second component enters for the L10 structure
c
          beta21=betas(3,i1)
          beta22=betas(3,i2)
          nrho=4
          irho=irhos(i1)
          rh21=rhof(r,beta21,re1,rozro1)
          irho=irhos(i2)
          rh22=rhof(r,beta22,re2,rozro2)
          rho22=(8./3.)*(rh21-rh22)**2
c   at the 2 atom
          rh2z=(4.*rh2+8.*rh1)
          t1av=(4*rh2*ts(2,i2)+8*rh1*ts(2,i1))/rh2z
          t2av=(4*rh2*ts(3,i2)+8*rh1*ts(3,i1))/rh2z
          t3av=(4*rh2*ts(4,i2)+8*rh1*ts(4,i1))/rh2z
          z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av,nn,z0s(i2))
          rhobarme=bar(rh2z,t2av*rho22,ibar2,z02*rozro2,dum1,dum2)
          rh2z=(rhobarme)/(z02*rozro2)
          ph22=phiid(r,i2)
c   at the 1 atom
          rh1z=(4.*rh1+8.*rh2)
          t1av=(4*rh1*ts(2,i1)+8*rh2*ts(2,i2))/rh1z
          t2av=(4*rh1*ts(3,i1)+8*rh2*ts(3,i2))/rh1z
          t3av=(4*rh1*ts(4,i1)+8*rh2*ts(4,i2))/rh1z
          z01=zbar(ibar1,z1,lattces(i1),t1av,t2av,t3av,nn,z0s(i1))
          rhobarme=bar(rh1z,t2av*rho22,ibar1,z01*rozro1,dum1,dum2)
          rh1z=(rhobarme)/(z01*rozro1)
          ph11=phiid(r,i1)
c
c       Eqn. (23)
c
          fa=frhoi(rh1z,asub1,esub1)
          fb=frhoi(rh2z,asub2,esub2)
          phif=erose(r,re12,alpha12,esub12,repul12,attra12)/4.0-
     1           fa/8.0 -fb/8.0 -ph11/4.0 -ph22/4.0
c	B32 structure
	elseif(all(it,jt).eq.'B32') then
          zij=8.0
          nrho=0
          irho=irhos(i1)
          rh1=rhof(r,beta1,re1,rozro1)
          irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
c
c   rhobarme is not just the
c   sum of zeroth components of the electron densities.  The
c   third component enters for the B32 structure
c
          rh20=(4.*rh2+4.*rh1)
          beta31=betas(4,i1)
          beta32=betas(4,i2)
          nrho=6
          irho=irhos(i1)
          rh31=rhof(r,beta31,re1,rozro1)
          irho=irhos(i2)
          rh32=rhof(r,beta32,re2,rozro2)
          t1av=(4*rh2*ts(2,i2)+4*rh1*ts(2,i1))/rh20
          t2av=(4*rh2*ts(3,i2)+4*rh1*ts(3,i1))/rh20
          t3av=(4*rh2*ts(4,i2)+4*rh1*ts(4,i1))/rh20
          z01=zbar(ibar1,z1,lattces(i1),t1av,t2av,t3av,nn,z0s(i1))
          z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av,nn,z0s(i2))
          if(ialloy.eq.0) then
           rho33=(32./9.)*(rh31-rh32)**2
          else
           denom=4*rh1*ts(3,i1)**2+4*rh2*ts(3,i2)**2
           rho33=(32./9.)*(rh31*ts(3,i1)-rh32*ts(3,i2))**2
           if(denom.gt.0.) rho33=rho33/denom*rh20
          endif
          rhobarme=bar(rh20,t3av*rho33,ibar2,z02*rozro2,dum1,dum2)
          rh1z=(rhobarme)/(z2*rozro2)
          rhobarme=bar(rh20,t3av*rho33,ibar1,z01*rozro1,dum1,dum2)
          rh2z=(rhobarme)/(z1*rozro1)
          ph11=phiid(r,i1)
          ph22=phiid(r,i2)
	  temp1=erose(r,re12,alpha12,esub12,repul12,attra12)
	  temp2=frhoi(rh2z,asub1,esub1)
	  temp3=frhoi(rh1z,asub2,esub2)
c
          fa=frhoi(rh2z,asub1,esub1)
          fb=frhoi(rh1z,asub2,esub2)
          phif=erose(r,re12,alpha12,esub12,repul12,attra12)/2.0-
     1           fa/4.0- fb/4.0 -ph11/2.0 -ph22/2.0
cstop
c
c     Diamond lattice with vacancies SiO_2
c
        elseif(all(it,jt).eq.'SiO2' )  then
            zij=4.
c  At an O
          roo=sqrt(8./3.)*r
          nrho=0
          irho=irhos(i1)
          rh0=rhof(r,betas(1,i1),re1,rozro1)
          irho=irhos(i2)
          rhoo0=rhof(roo,betas(1,i2),re2,rozro2)
          nrho=2
          irho=irhos(i1)
          rh1=rhof(r,betas(2,i1),re1,rozro1)
          irho=irhos(i2)
          rhoo1=rhof(roo,betas(2,i2),re2,rozro2)
          nrho=4
          irho=irhos(i1)
          rh2=rhof(r,betas(3,i1),re1,rozro1)
          irho=irhos(i2)
          rhoo2=rhof(roo,betas(3,i2),re2,rozro2)
          nrho=6
          irho=irhos(i1)
          rh3=rhof(r,betas(4,i1),re1,rozro1)
          irho=irhos(i2)
          rhoo3=rhof(roo,betas(4,i2),re2,rozro2)
          rh20=(2.*rh0+6.*rhoo0)
          t1av=(2*rh0*ts(2,i1)+6*rhoo0*ts(2,i2))/rh20
          t2av=(2*rh0*ts(3,i1)+6*rhoo0*ts(3,i2))/rh20
          t3av=(2*rh0*ts(4,i1)+6*rhoo0*ts(4,i2))/rh20
          z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av,nn,z0s(i2))
c         z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av)
          rho21=(2/sqrt(3.)*rh1+4/sqrt(2.)*rhoo1)**2
          rho22=3*(2./3.*rh2+2*rhoo2)**2+2*(2./3.*rh2)**2
     1          -(2*rh2+6*rhoo2)**2/3
          rho23=(2./3./sqrt(3.)*rh3-sqrt(2.)*rhoo3)**2
     1          +6*(2./3./sqrt(3.)*rh3-sqrt(2.)/2.*rhoo3)**2
     1          +6*(2./3./sqrt(3.)*rh3)**2
     1          -legend*(2/sqrt(3.)*rh3+4/sqrt(2.)*rhoo3)**2
          rhobarme=bar(rh20,t1av*rho21+t2av*rho22+t3av*rho23
     1                 ,ibar2,z02*rozro2,dum1,dum2)
          rhox=rhobarme/(z02*rozro2)
cprint *,'Si',r,roo,rh0,rh1,rh2,rh3
cprint *,'O',rhoo0,rhoo1,rhoo2,rhoo3
cprint *,' in phif Si',irho,arg,arg1,arg2,arg3,z01,rhsi,rozro1
c  At a Si
          irho=irhos(i2)
          nrho=0
          rhsi=rhof(r,betas(1,i2),re2,rozro2)
          rhsi=rhsi*zij
          nrho=6
          rh3=rhof(r,betas(4,i2),re2,rozro2)
          arg=(32.0/9.0)*ts(4,i2)*rh3**2
          z01=zbar(ibar1,z1,lattces(i1),ts(2,i1),ts(3,i1),ts(4,i1),
     1        nn,z0s(i1))
c         z01=zbar(ibar1,z1,lattces(i1),ts(2,i1),ts(3,i1),ts(4,i1))
            rhsi=bar(rhsi,arg,ibar1,z01*rozro1,dum1,dum2)
          rhsi=rhsi/(z01*rozro1)
          fa=frhoi(rhsi,asub1,esub1)
          fb=frhoi(rhox,asub2,esub2)
          phif=(3*erose(r,re12,alpha12,esub12,repul12,attra12)
     1          -fa-2*fb )/zij
cprint *,z01,z02,rhsi,rhox
cprint *, rhobarme,rh20,rozro2,rozro1
cprint *,phif
cstop
        elseif(all(it,jt).eq.'bqz' )  then
          zij=4.
c  At an O
          nrho=0
          irho=irhos(i1)
          rh0=rhof(r,betas(1,i1),re1,rozro1)
          nrho=2
          rh1=rhof(r,betas(2,i1),re1,rozro1)
          nrho=4
          rh2=rhof(r,betas(3,i1),re1,rozro1)
          nrho=6
          rh3=rhof(r,betas(4,i1),re1,rozro1)
          rh20=2.*rh0
          t1av=ts(2,i1)
          t2av=ts(3,i1)
          t3av=ts(4,i1)
          z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av,nn,z0s(i2))
c         z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av)
          rho21=0.21*rh1**2
          rho22=2.27*rh2**2
          rho23=(0.567-0.21*legend)*rh3**2
          rhobarme=bar(rh20,t1av*rho21+t2av*rho22+t3av*rho23
     1                 ,ibar2,z02*rozro2,dum1,dum2)
          rhox=rhobarme/(z02*rozro2)
c  At a Si
          irho=irhos(i2)
          nrho=0
          rhsi=rhof(r,betas(1,i2),re2,rozro2)
          rhsi=rhsi*zij
          nrho=6
          rh3=rhof(r,betas(4,i2),re2,rozro2)
          arg=(32.0/9.0)*ts(4,i2)*rh3**2
          z01=zbar(ibar1,z1,lattces(i1),ts(2,i1),ts(3,i1),ts(4,i1),
     1        nn,z0s(i1))
c         z01=zbar(ibar1,z1,lattces(i1),ts(2,i1),ts(3,i1),ts(4,i1))
          rhsi=bar(rhsi,arg,ibar1,z01*rozro1,dum1,dum2)
cprint *,'bqz',r,rhsi,rhobarme,rh20,rho21,rho22,rho23,rh0,i1,i2
          rhsi=rhsi/(z01*rozro1)
          fa=frhoi(rhsi,asub1,esub1)
          fb=frhoi(rhox,asub2,esub2)
          phif=(3*erose(r,re12,alpha12,esub12,repul12,attra12)
     1          -fa-2*fb )/zij
        elseif(all(it,jt).eq.'C11') then
          zij=10.0
          nrho=0
          irho=irhos(i1)
          rh1=rhof(r,beta1,re1,rozro1)
          irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
c
c   rhobarme is not just the 
c   sum of zeroth components of the electron densities.  All 
c   components enters for the C11 structure
c
	  nrho=2
          irho=irhos(i1)
          rh11=rhof(r,betas(2,i1),re1,rozro1)
          irho=irhos(i2)
          rh12=rhof(r,betas(2,i2),re2,rozro2)
	  nrho=4
          irho=irhos(i1)
          rh21=rhof(r,betas(3,i1),re1,rozro1)
          irho=irhos(i2)
          rh22=rhof(r,betas(3,i2),re2,rozro2)
	  nrho=6
          irho=irhos(i1)
          rh31=rhof(r,betas(4,i1),re1,rozro1)
          irho=irhos(i2)
          rh32=rhof(r,betas(4,i2),re2,rozro2)
c  at a Mo
	  rh10=10.*rh2
          z01=zbar(ibar1,z1,lattces(i1),ts(2,i2),ts(3,i2),ts(4,i2),
     1        nn,z0s(i1))
c         z01=zbar(ibar1,z1,lattces(i1),ts(2,i2),ts(3,i2),ts(4,i2))
	  rho12=(2./3.)*rh22**2
          rhobarme=bar(rh10,ts(3,i2)*rho12,ibar1,z01*rozro1,dum1,dum2)
          rh2z=rhobarme/(z01*rozro1)
c  at a Si
          rh20=(5.*rh2+5.*rh1)
          t1av=(5*rh2*ts(2,i2)+5*rh1*ts(2,i1))/rh20
          t2av=(5*rh2*ts(3,i2)+5*rh1*ts(3,i1))/rh20
	  t3av=(5*rh2*ts(4,i2)+5*rh1*ts(4,i1))/rh20
          z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av,nn,z0s(i2))
c         z02=zbar(ibar2,z2,lattces(i2),t1av,t2av,t3av)
          rho21=(rh11-rh12)**2
          rho22=(1./6.)*(rh21+rh22)**2
          rho23=(29./8.-legend)*(rh31-rh32)**2
          rhobarme=bar(rh20,t1av*rho21+t2av*rho22+t3av*rho23
     1                 ,ibar2,z02*rozro2,dum1,dum2)
          rh1z=rhobarme/(z02*rozro2)
	  ph22=phiid(r,i2)
c
c	Eqn. (23)
c
          fa=frhoi(rh2z,asub1,esub1)
          fb=frhoi(rh1z,asub2,esub2)
          phif=(3./10.)*(erose(r,re12,alpha12,esub12,repul12,attra12)
     1          -fa/3.0 -2.*fb/3.0 -(5./3.)*ph22)
        elseif( all(it,jt).eq.'MH' ) then
c
c
c
          zij=6
          rmm=sqrt(2.)*r
          nrho=0
          irho=irhos(i1)
          rhm=rhof(rmm,beta1,re1,rozro1)
          rh1=rhof(r,beta1,re1,rozro1)
          irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
          rh20=12*rhm+zij*rh2
c   t's at the M
          t1av=(zij*rh2*ts(2,i2)+12*rhm*ts(2,i1))/rh20
          t2av=(zij*rh2*ts(3,i2)+12*rhm*ts(3,i1))/rh20
          t3av=(zij*rh2*ts(4,i2)+12*rhm*ts(4,i1))/rh20
          z01=zbar(ibar1,z1,lattces(i1),t1av,t2av,t3av,nn,z0s(i1))
c         z01=zbar(ibar1,z1,lattces(i1),t1av,t2av,t3av)
          zh=zbar(ibar2,z2,lattces(i2),ts(2,i1),ts(3,i1),ts(4,i1),
     1       nn,z0s(i2))
c         zh=zbar(ibar2,z2,lattces(i2),ts(2,i1),ts(3,i1),ts(4,i1))
          rh1z=zij*rh1/(zh*rozro2)
          rh2z=(12*rhm+zij*rh2)/(z01*rozro1)
	  ph11=phiid(r,i1)
          fa=frhoi(rh2z,asub1,esub1)
          fb=frhoi(rh1z,asub2,esub2)
          phif=(2*erose(r,re12,alpha12,esub12,repul12,attra12)
     1          -fa-fb)/zij - ph11
        else
c
c     Cases for all opposite neighbors
c
	  nrho=0
          irho=irhos(i1)
          rh1=rhof(r,beta1,re1,rozro1)
          irho=irhos(i2)
          rh2=rhof(r,beta2,re2,rozro2)
          z01=zbar(ibar1,z1,lattces(i1),ts(2,i2),ts(3,i2),ts(4,i2),
     1        nn,z0s(i1))
          z02=zbar(ibar2,z2,lattces(i2),ts(2,i1),ts(3,i1),ts(4,i1),
     1        nn,z0s(i2))
	  b2=0.0
	  ph11=0.0
c
c     Diamond lattice
c
          if(all(it,jt).eq.'dia' )  then
            zij=4.
            rh1=rh1*zij
            rh2=rh2*zij
            beta3=betas(4,i1)
	    nrho=6
            irho=irhos(i1)
            rh3=rhof(r,beta3,re1,rozro1)
            arg=(32.0/9.0)*ts(4,i1)*rh3**2
            rh1=bar(rh1,arg,ibar2,zij*rozro1,dum1,dum2)
            beta3=betas(4,i2)
            irho=irhos(i2)
            rh3=rhof(r,beta3,re2,rozro2)
            arg=(32.0/9.0)*ts(4,i2)*rh3**2
            rh2=bar(rh2,arg,ibar1,zij*rozro2,dum1,dum2)
          elseif(all(it,jt).eq.'dim' )  then
c
c     dimer
c
            zij=1.
            beta1=betas(2,i1)
            beta2=betas(3,i1)
            beta3=betas(4,i1)
            irho=irhos(i1)
	    nrho=2
            rh11=rhof(r,beta1,re1,rozro1)
	    nrho=4
            rh21=rhof(r,beta2,re1,rozro1)
	    nrho=6
            rh31=rhof(r,beta3,re1,rozro1)
            t1=ts(2,i1)
            t2=ts(3,i1)
            t3=ts(4,i1)
            arg=t1*rh11**2+2./3.*t2*rh21**2+t3*(1-legend)*rh31**2
            rh1=bar(rh1,arg,ibar2,zij*rozro1,dum1,dum2)
            beta1=betas(2,i2)
            beta2=betas(3,i2)
            beta3=betas(4,i2)
            irho=irhos(i2)
	    nrho=2
            rh11=rhof(r,beta1,re2,rozro2)
	    nrho=4
            rh21=rhof(r,beta2,re2,rozro2)
	    nrho=6
            rh31=rhof(r,beta3,re2,rozro2)
            t1=ts(2,i2)
            t2=ts(3,i2)
            t3=ts(4,i2)
            arg=t1*rh11**2+2./3.*t2*rh21**2+t3*(1-legend)*rh31**2
            rh2=bar(rh2,arg,ibar1,zij*rozro2,dum1,dum2)
          elseif(all(it,jt).eq.'B1' )  then
c
c     B1 lattice
c
            zij=6.
            rh1=rh1*zij
            rh2=rh2*zij
            a2=sqrt(2.)
		if(nn) then
		p=a2*r*ntable/rcut
   		k=p
		  if(k.gt.ntable-1) then
		  ph11=0.
		  else
		  p=p-k
		  ph11=(1-p)*phitab(k,it,it)+p*phitab(k+1,it,it)
		  endif
		else
		ph11=phiid(a2*r,i1)
		endif
        	if(noscr.eq.0) then
		 cmx=cmax(i1,i2,i1)
		 cmn=cmin(i1,i2,i1)
          	 arg=1.
          	 if(cmx.gt.0.0) then
			fac=(cmx-arg)/(cmx-cmn)
		 else
			fac=0.
		 endif
          	 if(fac.ge.1) then
           	  xfac=0.
          	 elseif (fac.le.0.) then
           	  xfac=1.
          	 else
           	  xfac=xcut(fac)
          	 endif
        	else
          	 xfac=1.
        	endif
	      b2=-xfac**2*12./6.
              irho=irhos(i2)
	      rh2=rh2-zij*b2*rhof(r*a2,beta1,re1,rozro1)
          elseif(all(it,jt).eq.'B2' )  then
c
c     B2 lattice
c
            zij=8.
            rh1=rh1*zij
            rh2=rh2*zij
          endif
c   end of only opposite neighbors
          rh1z=rh1/(z02*rozro2)
          rh2z=rh2/(z01*rozro1)
	  fa=frhoi(rh2z,asub1,esub1)
	  fb=frhoi(rh1z,asub2,esub2)
          phif=(2*erose(r,re12,alpha12,esub12,repul12,attra12)
     1          -fa-fb)/zij +b2*ph11/2.0
        endif
c  end of reference structure
      endif
c  end of Mixed types
	phif=frac*phif+(1-frac)*(phizbl-fadimer-fbdimer)
      return
      end





c*****************************************************************************

      subroutine pwsfor
c
c     Find PairWiSeFORces
c
      include 'param.inc'
        include 'dyn.inc'
        include 'meam.inc'
        include 'wrkspc.inc'
        dimension rcos(3)
c
c     Increment for differentiation
c
        h=1.0e-5
        h=1.0e-4
c
c    find derivative of the screening function
c
        do 3000 i=1,natoms
          nlow = neiind(i-1)
          nhigh = neiind(i)
          ity = itype(i)
          do 3100 j=nlow+1,nhigh
            scree=scrnab(i,j-nlow)
c
            if (scree .le. scnres) go to 3100
c
            jindx = jneighs(j)
            jty=itype(jindx)
c
c     calculate distance between i and j
c
            dxij = -dneighs(1,j)
            dyij = -dneighs(2,j)
            dzij = -dneighs(3,j)
            rs=sqrt(dxij**2+dyij**2+dzij**2)
            rcos(1)=dxij/rs
            rcos(2)=dyij/rs
            rcos(3)=dzij/rs
c
c     calc dphif between pairs
c
            dphif=(phif(rs+h,ity,jty)-phif(rs-h,ity,jty))/(2.0*h)
            phiij=phif(rs,ity,jty)
3200        continue
c
c     Compile force on i and jindx atoms
c     Compile slocal on i and jindx atoms
c
c
c     note the minus sign on j term as the direc. cosine is opposite
c
            do 3300 l=1,3
              df=dphif*rcos(l)*scree-phiij*dscrnab(l,i,j-nlow)
              f(l,i)=f(l,i)+df
              slocal(l,1,i) = slocal(l,1,i)-df*dxij
              slocal(l,2,i) = slocal(l,2,i)-df*dyij
              slocal(l,3,i) = slocal(l,3,i)-df*dzij
3300        continue
3100      continue
3000    continue
      end



c*****************************************************************************

      subroutine rhosfor
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      include 'wrkspc.inc'
      dimension drcos(3,3),rcos(3),rr(3)
        common /angular/  dang1(natmax),dang2(natmax)
c
c     Compute forces due to density with d(screen)/dr terms.
c
      do 550 i=1,natoms
        if(c(i).le.0) go to 550
        ity=itype(i)
c
c     Define weighting factors etc. for atomic densities
c
	if(ialloy.eq.0) then
          t1=tav(1,i)/c(i)
          t2=tav(2,i)/c(i)
          t3=tav(3,i)/c(i)
	else
	  t1=0
	  t2=0
	  t3=0
	  if(tsqav(1,i).gt.0) t1=tav(1,i)/tsqav(1,i)
	  if(tsqav(2,i).gt.0) t2=tav(2,i)/tsqav(2,i)
	  if(tsqav(3,i).gt.0) t3=tav(3,i)/tsqav(3,i)
	endif
        ibar=ibarr(ity)
        ib=abs(ibar)
        call s(lattces(ity),s1,s2,s3)
        asub=asubs(ity)
        esub=esubs(ity,ity)
        rozeroi=rozros(ity)
        z=zs(ity)
        z0=zbar(ibar,z,lattces(ity),t1,t2,t3,nn,z0s(ity))
c       z0=zbar(ibar,z,lattces(ity),t1,t2,t3)
c
c     Note the division of rho by z0*rozeroi here ......
c
        roz=rho(i)/(z0*rozeroi)
        rho0=c(i)
        dfrho=dfrhoi(roz,asub,esub)/(z0*rozeroi)
c
c d rhobar / d rho0
c
          drbd0=rho(i)/rho0+rho0*dang1(i)
c
c d rhobar / d AA
c
          drbdAA=rho0*dang2(i)
c
c     Obtain the information about the neighbors of the given atom
c
        nlow = neiind(i-1)
        nhigh = neiind(i)
        do 500 j=nlow+1,nhigh
          scree=scrnab(i,j-nlow)
          if (scree .le. scnres) go to 500
          jindx = jneighs(j)
          jty=itype(jindx)
c
c     Find which neighbor of j is i. (Will seg fault on error.)
c
          llow = neiind(jindx-1)
          lngh = llow
3200      continue
            lngh = lngh + 1
          if (jneighs(lngh) .ne. i) goto 3200
c
c     Find interatomic distances, normalized directional vectors, etc.
c
          dxij = -dneighs(1,j)
          dyij = -dneighs(2,j)
          dzij = -dneighs(3,j)
          rs = sqrt(dxij*dxij + dyij*dyij + dzij*dzij)
          rr(1)=dxij
          rr(2)=dyij
          rr(3)=dzij
          rcos(1)=dxij/rs
          rcos(2)=dyij/rs
          rcos(3)=dzij/rs
          call dcmij(drcos,rr,rs)

c          if(i.eq.1) then
c             write(*,*) 'i=',i,' jindx=',jindx
c             write(*,*) 'rr=(',rr(1),',',rr(2),',',rr(3),')'
c             write(*,*) 'rcos=(',rcos(1),',',rcos(2),',',rcos(3),')'
c           write(*,*) 'drcos=(',drcos(1,1),',',drcos(2,1),',',drcos(3,1)
c           write(*,*) drcos(1,2),',',drcos(2,2),',',drcos(3,2)
c           write(*,*) drcos(1,3),',',drcos(2,3),',',drcos(3,3),')'
c          endif

c
c     Parameters for calculating embedding potentials
c
          re=res(jty,jty)
          rozero=rozros(jty)
          beta0=betas(1,jty)
          beta1=betas(2,jty)
          beta2=betas(3,jty)
          beta3=betas(4,jty)
        irho=irhos(jty)
        if(ialloy.eq.0) then
        fac1=1.
        fac2=1.
        fac3=1.
	else
        fac1=ts(2,jty)
        fac2=ts(3,jty)
        fac3=ts(4,jty)
        endif
	nrho=0
          rhs1=rhof(rs,beta0,re,rozero)
          drhs1=drhof(rs,beta0,re,rhs1)
	nrho=2
          rhs2=rhof(rs,beta1,re,rozero*fac1)
          drhs2=drhof(rs,beta1,re,rhs2)
	nrho=4
          rhs3=rhof(rs,beta2,re,rozero*fac2)
          drhs3=drhof(rs,beta2,re,rhs3)
	nrho=6
          rhs4=rhof(rs,beta3,re,rozero*fac3)
          drhs4=drhof(rs,beta3,re,rhs4)
c
          if(ibar.ge.0) then
            factor=0.0
          else
c
c   For zbar not = to z
c
            factor=(rho0/z)**2
          endif
	if(ialloy.eq.0) then
          ddd=((ts(2,jty)-t1)*(rhsq(1,i)-s1*factor)+
     1        (ts(3,jty)-t2)*(rhsq(2,i)-s2*factor)+
     1        (ts(4,jty)-t3)*(rhsq(3,i)-s3*factor))/rho0
	else
	  ddd=0
	  if(tsqav(1,i).gt.0) then
	    ddd= ddd
     1    +(ts(2,jty)-t1*ts(2,jty)**2)*(rhsq(1,i)-s1*factor)/tsqav(1,i)
	  endif
	  if(tsqav(2,i).gt.0) then
	    ddd= ddd
     1    +(ts(3,jty)-t2*ts(3,jty)**2)*(rhsq(2,i)-s2*factor)/tsqav(2,i)
	  endif
	  if(tsqav(3,i).gt.0) then
	    ddd= ddd
     1    +(ts(4,jty)-t3*ts(4,jty)**2)*(rhsq(3,i)-s3*factor)/tsqav(3,i)
	  endif
	endif
          tmp0 = (2.0 / 3.0) * cg(i)
          do 380 kk=1,3
            rcosj = rcos(kk)
c
c     Initialize accumulation variables
c
            drho0 = drhs1 * scree * rcosj
            drho1 = 0.0
            drho1g = 0.0
            drho1sg = 0.0
            drho2 = -tmp0 * drhs3 * scree * rcosj
            drho3 = 0.0
            drho0s = rhs1
            drho1s = 0.0
            drho2s = -rhs3 * tmp0
            drho3s = 0.0
            do 250 m=1,3
              rcosm = rcos(m)
              tmp1 = 2.0 * a(m,i)
              tmp1g = 2.0 * ag(m,i)
              drho1 = drho1 + tmp1 * (drhs2 * rcosj * rcosm +
     1                rhs2 * drcos(m,kk)) * scree
c            write(*,*) 'drho1: tmp1=',tmp1,' drhs2=',drhs2,' rhs2=',rhs2
              drho1g = drho1g + tmp1g * (drhs4 * rcosj * rcosm +
     1                rhs4 * drcos(m,kk)) * scree
              drho1s = drho1s + tmp1 * rhs2 * rcosm
              drho1sg = drho1sg + tmp1g * rhs4 * rcosm
              do 240 l = 1, 3
                rcosl = rcos(l)
                rml = rcosm * rcosl
                drmllm = drcos(m,kk) * rcosl + rcosm * drcos(l,kk)
                tmp2 = 2.0 * b(m,l,i)
                drho2 = drho2 + tmp2 * (drhs3 * rcosj * rml +
     1                  rhs3 * drmllm) * scree
                drho2s = drho2s + tmp2 * rhs3 * rml
                do 230 n = 1, 3
                  rcosn = rcos(n)
                  rmln = rml * rcosn
                  tmp3 = 2.0 * d(m,l,n,i)
                  drho3 = drho3 + tmp3 * (drhs4 * rcosj * rmln +
     1                    rhs4 * (drmllm * rcosn +
     1                    rml * drcos(n,kk))) * scree
                  drho3s = drho3s + tmp3 * rhs4 * rmln
 230            continue
 240          continue
 250        continue
c
            term1 = t1 * drho1 + t2 * drho2  +t3*(drho3-legend*drho1g)
            term2 = t1 * drho1s + t2 * drho2s+t3*(drho3s-legend*drho1sg)
            term1 = term1  + ddd * drho0
            term2 = term2  + ddd * drho0s
            
c           write(*,*) 't1=',t1,' drho1=',drho1,' t2=',t2,' drho2=',drho2
c           write(*,*) 't3=',t3,' drho3=',drho3,' t2=',t2,' drho2=',drho2
c           write(*,*) 'legend=',legend,' drho1g=',drho1g,' ddd=',ddd
c           write(*,*) 'rho0=',rho0,' drho1s=',drho1s,' drho2s=',drho2s
c           write(*,*) 'rho3s',rho3s,' drho1sg=',drho1sg

c
c
c     Compile force and slocal on i atom
c
            tmp4 = dscrnab(kk,i,j - nlow)
            df = dfrho * ((drho0 - drho0s * tmp4) * drbd0 +
     1           (term1 - term2 * tmp4) * drbdAA)
            f(kk,i) = f(kk,i) + df
            slocal(kk,1,i) = slocal(kk,1,i)-df*dxij
            slocal(kk,2,i) = slocal(kk,2,i)-df*dyij
            slocal(kk,3,i) = slocal(kk,3,i)-df*dzij

c            if(i.eq.1) then
c               write(*,*) 'i=',i,' jindx=',jindx,' kk=',kk,' df=',df
c             write(*,*) 'dfrho=',dfrho,' drho0=',drho0,' drho0s=',drho0s
c             write(*,*) 'tmp4=',tmp4,' drbd0=',drbd0,' term1=',term1
c             write(*,*) 'term2=',term2,' tmp4=',tmp4,' drbdAA=',drbAA             
c               return
c            endif
c
c     Compile force and slocal on jindx atom
c
            tmp4 = dscrnab(kk,jindx,lngh - llow)
            df = dfrho * ((drho0 + drho0s * tmp4) * drbd0 +
     1           (term1 + term2 * tmp4) * drbdAA)
            f(kk,jindx) = f(kk,jindx) - df
            
c           if((i.eq.1)) then
c              write(*,*) 'f(',jindx,',',kk,')=',f(kk,jindx),' df=',df
c           endif

            slocal(kk,1,jindx) = slocal(kk,1,jindx)-df*dxij
            slocal(kk,2,jindx) = slocal(kk,2,jindx)-df*dyij
            slocal(kk,3,jindx) = slocal(kk,3,jindx)-df*dzij
 380      continue
 500    continue
	if(rho(i).lt.0.0) then
	print *,' atom ',i,' forces ',f(1,i),f(2,i),f(3,i)
	endif
 550  continue
      return
      end

c**********************************************************************

      subroutine dscrfor
c
c     Calculate forces on screeing atoms
c
      include 'param.inc'
        include 'dyn.inc'
        include 'meam.inc'
        include 'wrkspc.inc'
        common /angular/  dang1(natmax),dang2(natmax)
        dimension rcos(3)

c        icount = 0
c
c     Increment for differentiation
c
        h=1.0e-5
        hh = h * h
        h2 = 2.0 * h
c
c     Compute forces due to density with d(screen(ij))/dr(k) terms.
c
        do 550 i=1,natoms
          if(c(i).le.0) go to 549
          ity=itype(i)
c
c     Define weighting factors etc. for atomic densities
c
        if(ialloy.eq.0) then
          t1=tav(1,i)/c(i)
          t2=tav(2,i)/c(i)
          t3=tav(3,i)/c(i)
        else
          t1=0
          t2=0
          t3=0
          if(tsqav(1,i).gt.0) t1=tav(1,i)/tsqav(1,i)
          if(tsqav(2,i).gt.0) t2=tav(2,i)/tsqav(2,i)
          if(tsqav(3,i).gt.0) t3=tav(3,i)/tsqav(3,i)
        endif
          ibar=ibarr(ity)
          ib=abs(ibar)
          call s(lattces(ity),s1,s2,s3)
          asub=asubs(ity)
          esub=esubs(ity,ity)
          rozeroi=rozros(ity)
          z=zs(ity)
          z0=zbar(ibar,z,lattces(ity),t1,t2,t3,nn,z0s(ity))
c         z0=zbar(ibar,z,lattces(ity),t1,t2,t3)
c
c     Note the division of rho by z0*rozeroi here .....
c
          roz=rho(i)/(z0*rozeroi)
          rho0=c(i)
          dfrho=dfrhoi(roz,asub,esub)/(z0*rozeroi)
c
c   rhobar=rho0*ang(rho0,AA)  AA = sum of rho1-3
c   rhobar'=rho0'*ang+ rho0*dang
c   dang1 = derivative of angular term wrt rho0
c   dang2 = derivative of angular term wrt AA
c
          drbd0=rho(i)/rho0+rho0*dang1(i)
          drbdAA=rho0*dang2(i)
c
          tmp0 = (2.0 / 3.0) * cg(i)
c
c     Obtain the information about the neighbors of the given atom
c
          nlow = neiind(i-1)
          nhigh = neiind(i)
          do 500 j=nlow+1,nhigh
             jindx=jneighs(j)
            scree=scrnab(i,j-nlow)
            if (scree .le. scnres .or. scree .ge. 1.0-scnres) go to 500
            jty=itype(jneighs(j))
            dxij = -dneighs(1,j)
            dyij = -dneighs(2,j)
            dzij = -dneighs(3,j)
            rij = dxij*dxij + dyij*dyij + dzij*dzij
            rs = sqrt (rij)
            rcos(1) = dxij / rs
            rcos(2) = dyij / rs
            rcos(3) = dzij / rs
            re=res(jty,jty)
            rozero=rozros(jty)
            beta0=betas(1,jty)
            beta1=betas(2,jty)
            beta2=betas(3,jty)
            beta3=betas(4,jty)
        irho=irhos(jty)
        if(ialloy.eq.0) then
        fac1=1.
        fac2=1.
        fac3=1.
        else
        fac1=ts(2,jty)
        fac2=ts(3,jty)
        fac3=ts(4,jty)
        endif
        nrho=0
          rhs1=rhof(rs,beta0,re,rozero)
        nrho=2
          rhs2=rhof(rs,beta1,re,rozero*fac1)
        nrho=4
          rhs3=rhof(rs,beta2,re,rozero*fac2)
        nrho=6
          rhs4=rhof(rs,beta3,re,rozero*fac3)
            if(ibar.ge.0) then
              factor=0.0
            else
c
c   For zbar not = to z
c
              factor=(rho0/z)**2
            endif
        if(ialloy.eq.0) then
          ddd=((ts(2,jty)-t1)*(rhsq(1,i)-s1*factor)+
     1        (ts(3,jty)-t2)*(rhsq(2,i)-s2*factor)+
     1        (ts(4,jty)-t3)*(rhsq(3,i)-s3*factor))/rho0
        else
          ddd=0
          if(tsqav(1,i).gt.0) then
            ddd= ddd
     1    +(ts(2,jty)-t1*ts(2,jty)**2)*(rhsq(1,i)-s1*factor)/tsqav(1,i)
          endif
          if(tsqav(2,i).gt.0) then
            ddd= ddd
     1    +(ts(3,jty)-t2*ts(3,jty)**2)*(rhsq(2,i)-s2*factor)/tsqav(2,i)
          endif
          if(tsqav(3,i).gt.0) then
            ddd= ddd
     1    +(ts(4,jty)-t3*ts(4,jty)**2)*(rhsq(3,i)-s3*factor)/tsqav(3,i)
          endif
        endif
c
c    Pair potential for i and j atoms
c
            phiij=phif(rs,ity,jty)
c            icount = icount + 1

            do 395 l=nlow+1,nhigh

              if (l.eq.j) goto 395
              lindx=jneighs(l)
              lty=itype(lindx)  
              dxil = -dneighs(1,l)
              dyil = -dneighs(2,l)
              dzil = -dneighs(3,l)
              ril = dxil*dxil + dyil*dyil + dzil*dzil
              rjl = rij + ril -
     .              2.0 * (dxij*dxil + dyij*dyil + dzij*dzil)
c
c    Find proper coordinates for compilation of slocal
c
              dxlm = 2.0 * dxil - dxij
              dylm = 2.0 * dyil - dyij
              dzlm = 2.0 * dzil - dzij
c
c     Find screening
c
              scr1=dscrn(rij,ril,rjl,ity,jty,lty)
c
c
c     See if derivative of screening term is non-zero
c
              if (scr1.lt. scnres .or. scr1 .gt. 1.0-scnres) go to 395
c
c     If so, resort to numeric differentiation
c
              do 3400 kk=1,3
                dil = -dneighs(kk, l)
                dij = -dneighs(kk, j)
c
c     Positive increment in xl's
c  
                rpil = ril + h2 * dil + hh
                rpjl = rjl + h2 * (dil - dij) + hh
                scrp1=dscrn(rij,rpil,rpjl,ity,jty,lty)
c
c     Negative increment in xl's
c      
                rmil = ril - h2 * dil + hh
                rmjl = rjl - h2 * (dil - dij) + hh
                scrm1=dscrn(rij,rmil,rmjl,ity,jty,lty)

c                icount = icount + 2
c
c     Find derivative
c
                dscrnabk=0.5*((scrp1-scrm1)*scree)/(scr1*h)
                rcosj = rcos(kk)
                
c
c     Initialize accumulation variables
c
                drho0s = rhs1
                drho1s = 0.0
                drho1sg = 0.0
                drho2s = -rhs3 * tmp0
                drho3s = 0.0
                do 250 m=1,3
                  rcosm = rcos(m)
                  tmp1 = 2.0 * a(m,i)
                  tmp1g = 2.0 * ag(m,i)
                  drho1s = drho1s + tmp1 * rhs2 * rcosm
                  drho1sg = drho1sg + tmp1g * rhs4 * rcosm
                  do 240 ll = 1, 3
                    rcosl = rcos(ll)
                    rml = rcosm * rcosl
                    tmp2 = 2.0 * b(m,ll,i)
                    drho2s = drho2s + tmp2 * rhs3 * rml
                    do 230 n = 1, 3
                      rcosn = rcos(n)
                      rmln = rml * rcosn
                      tmp3 = 2.0 * d(m,ll,n,i)
                      drho3s = drho3s + tmp3 * rhs4 * rmln
 230                continue
 240              continue
 250            continue
c
            term2 = t1 * drho1s + t2 * drho2s+t3*(drho3s-legend*drho1sg)
            term2 = term2  + ddd * drho0s
c
c     Sum forces on the cental atom
c
                df = dscrnabk *
     .               (0.5 * phiij +
     .                dfrho * (drbd0 * drho0s + drbdAA * term2))
                f(kk,lindx) = f(kk,lindx) - df
                slocal(kk,1,lindx) = slocal(kk,1,lindx)-df*dxlm
                slocal(kk,2,lindx) = slocal(kk,2,lindx)-df*dylm
                slocal(kk,3,lindx) = slocal(kk,3,lindx)-df*dzlm

c            if(lindx.eq.41) then
c             write(*,*) 'i=',i,' jindx=',jindx,' lindx=',lindx,' kk=',kk
c               write(*,*) 'df=',df
c               write(*,*) 'f(',kk,',',lindx,')=',f(kk,lindx)
c               write(*,*) 'dscrnabk=',dscrnabk,' phiij=',phiij
c               write(*,*) 'dfrho=',dfrho,' drh0s=',drho0s
c             write(*,*) 'term2=',term2,' drbdAA=',drbdAA,' drbd0=',drbd0
c              write(*,*) 'scr1=',scr1,' scrp1=',scrp1,' scrm1=',scrm1
c              write(*,*) 'rij=',rij,' rpil=',rpil,' rpjl=',rpjl
c              write(*,*) 'rij=',rij,' rmil=',rmil,' rmjl=',rmjl
c              write(*,*) 'ril=',ril,' rjl=',rjl
c             write(*,*) 'rhs1=',rhs1,' rhs2=',rhs2,' rhs3=',rhs3
c     1                 ,' rhs4=',rhs4
c               return
c            endif

3400          continue
 395        continue
 500      continue
 549      continue
 550    continue

c        write(*,*) 'dscrfor: icount=',icount
        return
      end
c*****************************************************************************
      subroutine dcmij(drcos,rr,rs)
        implicit real*8 (a-h,o-z)
c
c     Find direction cosines
c
        dimension drcos(3,3),rr(3)
        r1 = 1.0 / rs
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
      function zbar(ibar,z,lat,t1,t2,t3,nn,z0s)
        implicit real*8 (a-h,o-z)
        character *4 lat
	logical nn
        if(ibar.ge.0) then
          zbar=z
        elseif (nn) then
	  zbar=z0s
	else
          call s(lat,s1,s2,s3)
          zbar=bar(z,s1*t1+s2*t2+s3*t3,ibar,z,dum1,dum2)
        endif
        return
      end
c*****************************************************************************
      subroutine s(lat,s1,s2,s3)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
        character *4 lat
        s1=0.0
        s2=0.0
        s3=0.0
        if(lat .eq. 'dia') then
          s3=32.0/9.0
        elseif (lat .eq. 'hcp') then
          s3=1.0/3.0
        elseif(lat .eq. 'dim') then
          s1=1.0
          s2=2.0/3.0
          s3=1.0
	s3=s3-legend*s1
        endif
        return
      end                               



