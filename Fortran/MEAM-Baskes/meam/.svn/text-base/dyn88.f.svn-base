c   dyn88.f
c$Author: baskes $
c$Date: 2004/09/29 19:26:24 $
c$Source: /home/baskes/programs/MEAM/v88/RCS/dyn88.f,v $
c$Revision: 1.30 $
c$Log: dyn88.f,v $
cRevision 1.30  2004/09/29 19:26:24  baskes
cfixed B1 + NN
c
cRevision 1.29  2004/02/24 16:10:32  baskes
cadded ialloy
c
cRevision 1.28  2003/12/09 21:17:58  baskes
cmerged with dyn96
c
cRevision 1.27  2003/12/09 17:04:25  baskes
cupdating dyn96 compatability
c
cRevision 1.26  2003/12/06 04:20:22  baskes
cupdating
c
cRevision 1.25  2003/01/18 18:02:13  baskes
cfixed initialization of 2NN MEAM with iscr=1
c
cRevision 1.24  2002/11/04 18:28:39  baskes
cchanged print
c
cRevision 1.23  2002/04/26 20:55:35  baskes
cadded vcm
c
cRevision 1.22  2002/04/07 01:24:36  baskes
cmoved restart stuff to utils
cadded constant velocity of periodic vector
c
cRevision 1.21  2002/02/22 19:48:29  baskes
cfixed header in output
c
cRevision 1.20  2002/01/04 23:15:14  baskes
cmodified averaging
c
cRevision 1.19  2001/12/04 22:40:36  baskes
cCM fixed
c
cRevision 1.18  2001/12/04 18:52:25  baskes
cadd icm zero CM velocity
c
cRevision 1.17  2001/09/17 16:25:04  baskes
cfixed NN
cadded ipscree
c
cRevision 1.16  2001/09/09 14:28:43  baskes
cadded SiO2
c
cRevision 1.15  2001/06/12 15:35:01  baskes
cadded acf
c
cRevision 1.14  2001/05/16 21:04:56  baskes
cfix arg for hcp
covercome z inconsistency
c
cRevision 1.13  2001/05/15 21:33:45  baskes
cremoved default repuls
caddeed 2NN MEAM
cadded table
c
cRevision 1.12  2001/02/27 18:29:37  baskes
cfixed lattces
c
cRevision 1.11  2001/01/17 01:02:30  baskes
cfull slocal
c
cRevision 1.10  2001/01/08 23:31:18  baskes
cfixed real seconds
c
cRevision 1.9  2001/01/08 22:59:47  baskes
cfixed header
c
c**********************************************************************
c
c               dynamic code version 8.7 (unix) msd/smf  8/30/99
c      (temperature control: tmpcard & tmprgn added to input namelist)
c	meafile added to read meacard 3/20/97
c	fixed hcp initialization 4/15/97
c	8/8/97 added legendre correction
c	10/4/97 added repulsive term to Rose function
c	3/26/98 added attrac to erose
c	11/13/98 added iconf=4 (energy output)
c       monoclinic structure added 8/9/99
c	CM velocity fixed  1/21/00
c	change in fix for negative rho^2 3/11/00
c
c*******************************************************************
c               conventions
c       9000+ reserved for format statement #'s
c*******************************************************************
c
	block data
	include 'param.inc'
	include 'dyn.inc'
	include 'meam.inc'
c	character*80 meamf,forcef,printf,rstrtf,conff
c	character*24 timestr
	data eqtim/0.0/
	data itable/0/,nn/.false./
	data rscrmax /0.0/,ipscree/0/
	data dt,tol,iaccur,inte,accbnd/0.001,0.001,0,1,1000.,400./
	data ipinter/0/,ipatoms/3/,ipitera/-1/,irstrt/0/,iconf/0/,
     1    iconst/0/
c	data contin/.false./,lastconf/.true./
	data xminp,yminp,zminp/3*-1.e6/,xmaxp,ymaxp,zmaxp/3*1.e6/
        data nminp,nmaxp/1,10000/
	data typep/0/,sortp/.false./
        data legend /0.0/
	data ialloy/0/
	end

      subroutine dyn88
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      logical contin,lastconf
c      character*80 printf,rstrtf,conff,meamf,forcef,statf
c      character*24 timestr
	common/free/ ifree
	common/print1/iustat
      namelist /prntcard/ ipinter,ipatoms,ipitera,printf,rstrtf,conff,
     1    iconst,xminp,yminp,zminp,xmaxp,ymaxp,zmaxp,nminp,nmaxp,
     2    typep,sortp,iconf,meamf,forcef,ipscree,iout,statf
      namelist /continue/ contin,lastconf
      data rstrtf,conff,forcef/3*'none'/,printf/'dynprint'/
c
c  determine the initial timing values
c
      call initsec()
c
c  read the parameters controlling output
c
      read(5,prntcard)
c
c  determine which output files are desired and open them
c
      open(unit=6,file=printf,form='FORMATTED',status='UNKNOWN')
c
	if(ifree.eq.0) then
      write(6,9100)
9100  format(' *******  dynamo version 8.7 (unix) *******'/,
     1   '$Revision: 1.30 $',/'$Date: 2004/09/29 19:26:24 $')
	else
      if (iout.eq.1) then
         iustat=6
      else
         iustat=23
         open(unit=23,file=statf,form='FORMATTED',status='UNKNOWN')
      endif
      write(6,9107)
9107  format(' *******  dynamo version 9.6 (unix) *******'/,
     1   '$Revision: 1.30 $',/'$Date: 2004/09/29 19:26:24 $')
	endif
c
      write(6,9005)printf
9005  format(/,/,1x,'opened print file named ',a)
c
c  print out the  real time for the beginning of the job
c
      call mytime(timestr)
      write(6,9101) timestr
9101  format(/,/,1x,'begin run:',a24)
c
c	open restart and configuration files
c
      if (.not.(rstrtf(1:4).eq.'NONE'.or.rstrtf(1:4).eq.'none')) then
         irstrt = 1
         open(unit=21,file=rstrtf,form='FORMATTED',status='UNKNOWN')
         write(6,9006)rstrtf
9006     format(1x,'opened restart file named ',a)
      end if
      if (.not.(conff(1:4).eq.'NONE'.or.conff(1:4).eq.'none')) then
         if (iconf.eq.0) iconf = 1
         open(unit=20,file=conff,form='UNFORMATTED',status='UNKNOWN')
         write(6,9007)conff,iconf
9007     format(1x,'opened configuration file named ',a,/
     $             ' configuration file type:',i2)
       else
         iconf = 0
       end if
      if (.not.(forcef.eq.'NONE'.or.forcef.eq.'none')) then
         iforce = 1
         open(unit=30,file=forcef,form='FORMATTED',status='UNKNOWN')
         write(6,9008) forcef
9008     format(1x,'opened force file named ',a)
        else
        iforce=0
      end if
c
c  call dynamo to perform the computation
c
1000  continue
      call dynamo(contin,lastconf)
      contin = .false.
      read(5,continue,end=9000,err=9000)
      if (.not.contin) goto 9000
c
c  mark the output file
c
      write(6,9201)
9201  format(1x,/,1x,'******  continuation run  *******',/,1x)
c
c  close and reopen the configuration file and restart file
c
      if (irstrt.eq.1) then
         close(21)
         open(unit=21,file=rstrtf,form='FORMATTED',status='UNKNOWN')
         write(6,9006)rstrtf
      end if
      if (iconf.ne.0) then
         close(20)
         open(unit=20,file=conff,form='FORMATTED',status='UNKNOWN')
         write(6,9007)conff
      end if
      goto 1000
9000  continue
c
c  compute final timing information
c
      t2tot = seconds()
      t2cpu = 0.
      t2io = 0.
      t2sys = 0.
c
      write(6,9010) t2tot,t2cpu,t2io,t2sys
9010  format(/,/,1x,'timing:',7x,'total',9x,'cpu',10x,'io',6x,'system',
     1  /,8x,4f12.3)
      call mytime(timestr)
      write(6,9102) timestr
9102  format(1x,'end of run:',a24)
c
      if (inte.ne.-1.or.mnterm.eq.0) then
         stop ''
       else if (mnterm.eq.1) then
         stop '*nfmax*'
       else
         stop '*timeup*'
       endif
      end
c************************************************************************
c
c
c       this is the main routine
c
      subroutine dynamo(contin,lastconf)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
        common /err/  ibarer,ibare
	common/head/ header
	common/free/ifree
      common /average1/  potbar,potrms,vcm(3),potrefbar,potrefrms,nfree
      logical contin,lastconf
      character*80 header
      data header /'        '/
      data nmeth/2/,dradn/-1.0/
	data zero/0.0001/
	data ncalcacf/1/
      namelist /headcard/ header
      namelist /neicard/ nmeth,dradn
      namelist /avecard/ eqtim,nbuff,nval,nlimitacf,ncalcacf
c
c
c  clear the force statistics variables
c
	sqrtz=sqrt(zero)
      nforce = 0
      ngtlst = 0
      nneimx = 0
      mxlstu = 0
      mxnnei = 0
      frctmx = -1.e5
      frctmn = 1.e5
      gnetmx = -1.e5
      gnetmn = 1.e5
c
c  read and print the descriptive header for this run
c
      if (.not.contin) read(5,headcard)
      write(6,9110) header
9110  format(' ****************************************** ',//,
     11x,a,//,
     1' ******************************************  ')
c
c describe the output control in effect
c
      if(.not.contin)write(6,9120)ipinter,ipatoms,ipitera,iconst
9120  format(/,/,'  print control:  ',/,
     1'    print interaction functions: ',i4,/,
     2'    print atom positions: ',i4,/,
     3'    print iteration info: ',i4,/,
     4'    print full constraint list: ',i4)
c
c  print the units convention
c
      if(.not.contin)write(6,9130)
 9130 format(/,/,' explanation of units',/,' positions in angstroms',/,
     1' time in picoseconds',/,' velocities in angstroms/picosecond',/,
     2' energies in ev',/,' temperature in degrees kelvin',/,
     3' pressure in bars ',/,
     4' mass in ev-psec^2/angstroms^2',/,/)
c
c       set up interactions
c
      if (.not.contin) then
         call inter
c         call modint
       end if
	if(ipinter.eq.1) call prphi
c
c  determine the initial positions and velocities
c
      call initlat(contin,lastconf)
c
c  read in the type of boundaries to be used
c
	call getbnd
c
c  determine the type of neighbor finding method to be used
c
      if (.not.contin) read(5,neicard)
      write(6,9611) nmeth
9611  format(1x,/,'****** using neighbor method: ',i2)
c
c  if nmeth=2, set up the storage required by the neighbor list,
c  initialize the values of the neighbor list index and compute the
c  cut off radii for the neighbor list using the default
c  expression for dradn if it is not specified
c
c  nmeth=2
      if (abs(nmeth).eq.2) then
         nnstrt = 0
         nnindx(0) = nnstrt
         nnindx(1) = nnstrt
         newlst = 1
         if (dradn.lt.0.0) dradn = 0.1*sqrt(rcutsq)
c         if (dradn.lt.0.0) dradn = 0.35*sqrt(float(natoms)/500.0)
         rctsqn = (sqrt(rcutsq) + dradn)**2
         write(6,9612) rctsqn,dradn
9612  format(1x,/,'    rctsqn:',g13.5,'  dradn:',g13.5)
       end if
c
c  nmeth = 3
c
      if (nmeth.eq.3) then
         newlst = 1
         natsav = natoms
         if (dradn.lt.0.0) dradn = 0.35
         rctsqn = (sqrt(rcutsq) + dradn)**2
         write(6,9612) rctsqn,dradn
       endif
c
c       calculate energy of initial lattice without defects
c
      call calce(ekin1,eperf,etot1,temper,press,1,t0)
c
c       read in defects
c
      write(6,9140)
 9140 format(/,/' ******  reading in defects ')
      call readef(ndef)
c
c  if defects have been created, force the determination of a
c  new neighbor list
c
      if (ndef.ne.0) newlst = 1
      write(6,9150)ndef
 9150 format('  finished reading in ',i5,' defects ')
c
c  delete vacancies
c
      call delvac
c
c  check that all types used have been defined
c
      call chktyp
c
c  make sure 0.lt.natoms.le.natmax
c
      if(natoms.le.0)then
         write(6,*) ' there are no atoms'
         stop
      endif
      if(natoms.gt.natmax)then
         write(6,9160)natoms,natmax
 9160    format(' natoms=',i10,' is greater than natmax=',i10)
         stop
      endif
c
c  read in the fixcard constraints
c
      call setfix(contin,lastconf,nfixcard)
c
c
c  read in and set up the temperature control
c
      call settmp(contin)
c
c  initiaze the averaging variables
c
      if (.not.contin) read(5,avecard)
      if (.not.contin) write(6,avecard)
c       calculate initial state
c
      if(nfixcard.ne.0.or.ndef.ne.0)then
          call calce(ekin1,pot1,etot1,temper,press,1,t0)
      else
          pot1 = eperf
      endif
      volint = perlen(1)*perlen(2)*perlen(3)/natoms
      write(6,9205)ekin1,pot1,etot1,temper,press,volint
 9205 format(' ******   initial state: ',/,
     1'  kinetic, potential and total energies ',3e17.9,/,
     2'  temperature:',g13.6,'  pressure:',g13.6,'  volume:',g13.6)
	if(ifree.ne.0) call freevel
c
c       print out types, positions, velocities, electron densities,
c       forces and energies on particles
c
      write(6,9170)natoms
 9170 format(/,/1x,i5,' particles ')
      if(abs(ipatoms).eq.1.or.abs(ipatoms).eq.3) call patoms
c
c  write the job parameters to the output file 20 for use by the
c  analysis programs
c
      if (iconf.ne.0) then
        write(20) header
        write(20) natoms,iconf
        write(20) (perub(i),i=1,3),(perlb(i),i=1,3)
       end if
c
c  call newton to integrate the equations of motion
c               or to minimize the potential
c
	nskipd=0
	nfree=0
	call average(-1,t)
	call average(0,t)
	if(ifree.eq.0) then
      call newton(contin)
	else
      call frnewton(contin)
	endif
c
c  if a minimization, put the final positions in a configuration file
c  if there is one
c
      if (iconf.ne.0.and.inte.eq.-1) then
c
      if (iconf.eq.1) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,3),itype(j),j=1,natoms)
       end if
      if (iconf.eq.2) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
       end if
      if (abs(iconf).eq.3) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,3),itype(j),j=1,natoms)
	if(iconf.gt.0) then
         write(20) (((slocal(i,j,k),i=1,j),j=1,3),k=1,natoms)
	else
         write(20) (((slocal(i,j,k),i=1,3),j=1,3),k=1,natoms)
	endif
       end if
      if (iconf.eq.4) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,6),itype(j),e(j),j=1,natoms)
       end if
c
      endif
c  write the force file
        if(iforce.eq.1) then
        write(30,9777) ((f(i,j),i=1,3),j=1,natoms)
9777    format(3e20.8)
        endif
c
c  if this run is just to calculate energy (inte=2), skip down to
c  timing information
c
c      if(inte.eq.2)goto 5000
c
c  print out the final position, velocities, and forces as well as
c  the final energies, temerature, and pressure
c
	if(ifree.ne.0) volintot=volint*natoms
	call finalpr(ekin1,pot1,etot1,eperf,volint)
c
c  compute and output the average values.
c
        if(inte.eq.1)  call average(3,t)
c
      if (nfree.gt.0) call freepot
c
c  create the restart file (21)
c
c       use format compatible with vax/creator capability
      if (irstrt.eq.1) then
	rewind(21)
      write(21,9501) header
9501  format(a)
      write(21,9502) natoms,ntypes,tend,beta
9502  format(2i10,e15.8,f10.5)
      write(21,9503) (perub(i),i=1,3),(perlb(i),i=1,3)
9503  format(3e25.16)
      write(21,9504) (amass(i),ielement(i),i=1,ntypes)
9504  format(e25.16,i10)
      write(21,9505) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
9505  format(3e25.16/3e25.16/i10)
      write(21,9503) (bndvel(i),i=1,3)
      end if
c
c  output the statistics on call to force
c
c  branch to here if only calculated energy once (inte=2)
5000  continue
      write(6,9401)
9401  format(/,/,1x,'****** statistics on FORCE and GNEIGH')
      write(6,9402)nforce
9402  format(1x,'total number of calls to FORCE:',i5)
      if (abs(nmeth).eq.2) write(6,9403) ngtlst
9403  format(1x,'total number of neighbor list updates:',i5)
	if(natoms.gt.2) then
      mxtmp = 2*int(1.+float(mxlstu)/float(natoms-2))
      mxnnei = max0(mxnnei,mxtmp)
      write(6,9404) nneimx,mxnnei
9404  format(1x,'maximum number of actual neighbors found:',i5,/,
     1       1x,'maximum number of possible neighbors found:',i5)
	endif
      write(6,9405)frctmx,frctmn
9405  format(1x,'maximum and minimum times for calls to FORCE:',/,
     1       5x,2g15.5)
      write(6,9406)gnetmx,gnetmn
9406  format(1x,'maximum and minimum times for calls to GNEIGH:',/,
     1       5x,2g15.5)
	write(6,9407) rscrmax
9407	format(1x,'maximum distance between partially',
     1	' unscreened atoms: ',f10.5)
	write(6,9408) ibarer,ibare
9408	format(1x,'number of negative square root calls',i5/
     1   'number in last call to force',i5)
      return
      end
c
c
c***************************************************************
c	prints out phi
c
	subroutine prphi
      include 'param.inc'
      include 'dyn.inc'
      dimension p(3,3)
	ntype=min(ntypes,3)
	write(6,9000) ((it,jt,jt=it,ntype),it=1,ntype)
9000	format('  pair potentials'/5x,'r',4x,6(4x,2i1,4x))
	do 9809 i=1,100
        r=.05*i
	do 1 it=1,ntype
	do 1 jt=it,ntype
1	p(it,jt)=phif(r,it,jt)
        write(6,9808) r,((p(it,jt),jt=it,ntype),it=1,ntype)
9808    format(7f10.4)
9809    continue         
	return
	end
c***************************************************************
c
c  patoms prints out the atom positions etc.
c
      subroutine patoms
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      include 'wrkspc.inc'
      dimension rmult(natmax),dum1(natmax),idum1(natmax),idum2(natmax)
        naver=0
        ave=0.
c
c  sort the output list if requested
c
      do 1000 i = 1,natoms
      idum1(i) = i
1000  idum2(i) = i
      if (.not.sortp) goto 2000
      permax = max(perlen(1),perlen(2),perlen(3))
      if(permax.ge.1000.)then
         xbot = 1.e6
         ybot = 1.e6
         zbot = 1.e6
         do 1005 i = 1,natoms
         xbot = min(xbot,rv(1,i))
         ybot = min(ybot,rv(2,i))
         zbot = min(zbot,rv(3,i))
1005     continue
         call m01aaf(perlen,1,3,idum1,idum2,ifail)
         rmult(1) = 10**(4*idum1(1)-3)
         rmult(2) = 10**(4*idum1(2)-3)
         rmult(3) = 10**(4*idum1(3)-3)
         do 1010 i=1,natoms
         dum1(i) = rmult(1)*(rv(1,i)-xbot)
     1           + rmult(2)*(rv(2,i)-ybot)
     2           + rmult(3)*(rv(3,i)-zbot)
1010     continue
      else
         do 1020 i=1,natoms
         dum1(i) = rv(1,i)**2 + rv(2,i)**2 + rv(3,i)**2
1020     continue
      endif
      call m01aaf(dum1,1,natoms,idum1,idum2,ifail)
      do 1030 i = 1,natoms
1030  idum2(idum1(i)) = i
      if(ifail.ne.0)write(6,9020)ifail
 9020 format('    trouble in sorting routine   ',i5)
2000  continue
c
c
c  now print out the atoms in the specified region and/or type
c
      write(6,9180)
 9180 format(1x,'atomic data:',
     1       'number, type, positions, energy')
      if (ipatoms.gt.0) write(6,9181)
9181  format(13x,'velocities, kinetic energy (if not minimization)',/,
     2   13x,'forces',/,
     3   13x,'el. density, f(rho), phi(r)',/,
     4   13x,'rho0, rho1^2, rho2^2, rho3^2'/)
	nmaxp=min(nmaxp,natoms)
      do 2100 i = nminp,nmaxp
      j = idum2(i)
      if (rv(1,j).lt.xminp.or.rv(1,j).gt.xmaxp) goto 2100
      if (rv(2,j).lt.yminp.or.rv(2,j).gt.ymaxp) goto 2100
      if (rv(3,j).lt.zminp.or.rv(3,j).gt.zmaxp) goto 2100
      if (typep.ne.0.and.itype(j).ne.typep) goto 2100
      write(6,9190) j,itype(j),rv(1,j),rv(2,j),rv(3,j),e(j)
9190  format(1x,i5,1x,i2,4g15.7)
        naver=naver+1
        ave=ave+e(j)
      if (ipatoms.lt.0) goto 2100
      if (inte.ne.-1) then
         akin = 0.5*amass(itype(j))*(rv(4,j)**2+rv(5,j)**2+rv(6,j)**2)
         write(6,9191) rv(4,j),rv(5,j),rv(6,j),akin
9191     format(9x,4g15.5)
       endif
      write(6,9192) (f(kk,j),kk=1,3)
9192  format(9x,3g15.5)
      ity = itype(j)
      write(6,9193) rho(j),plocal(j),fp(j)
9193  format(9x,3g15.5)
	write(6,9191) c(j),rhsq(1,j),rhsq(2,j),rhsq(3,j)
2100  continue
c
c  print out the periodic lengths
c
      write(6,9200) (perlen(i),i=1,3)
 9200 format(1x,'perlen:',3g15.7)
        if(naver.gt.0) ave=ave/naver
        write(6,9201) ave
9201    format(' average energy ',e20.12)
      return
      end
c***************************************************************
c
c  this routine sets the parameters defining the interactions
c
      subroutine inter
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
       common /fr/ frcut,xncut,xmcut
	common /alloy/ all(nelmax,nelmax)
	common /scre/ noscr
        common /czbl/ xzbl,xzbl0
c      character*80 header,meamf,meafile
      character*4 name,all,lattce
      character*7 kodes(nelmax)
	data xncut,xmcut/2.,4./
      data conmas/1.0365e-4/
      data rz0/1.0/,alpha0/1.0/
        data xzbl,xzbl0/-3.,-1./
	data meamf/'meamf'/
	data meafile/'none'/
	data noscr/0/
	data hartre/27.2/,bohr/0.529/
	data frcut/0.9/
      namelist /meacard/ ntypes,lattces,enames,esubs,alats,alphas,asubs,
     1                   betas,ts,deltas,zs,rozros,kderiv,kodes,rcut,
     2                   rz0,alpha0,all,res,cmin,cmax,itable,nn,ialloy,
     3                   noscr,meafile,legend,repuls,attrac,xzbl,xzbl0,
     4			frcut,xncut,xmcut
c
c     subroutine to read in meam data and compute some parameters
c
c
	write(*,*) 'call inter...'
      kmeamf=13
      open (unit=kmeamf,file=meamf,status='old',form='formatted')
c
 105  format(a)
 110  format(a)
	do 777 n=1,nelmax
	rozros(n)=0.
	do 777 m=1,nelmax
	alphas(m,n)=0.
	repuls(m,n)=0.
	attrac(m,n)=0.
	res(m,n)=0.
	do 777 l=1,nelmax
	cmin(l,m,n)=2.0
	cmax(l,m,n)=2.8
777	continue
	write(*,*) ' inter... (2)'
        read  (5,meacard)
	write(*,*) ' inter... (3)'
	if(meafile.ne.'none') then
	kmeaf=55
        open (unit=kmeaf,file=meafile,status='old',form='formatted')
	read(kmeaf,meacard)
	endif
      rcutsq=rcut**2
	write(*,*) ' inter... (4)'
      write (6,100) legend
 100  format(//,1x,' ******  MEAM data',/' legendre ', f10.5)
      write (6,101) xzbl,xzbl0
 101  format(//,1x,'  ZBL limits ', 2f10.5)
      write (6,115) ntypes,kderiv,rcut,frcut,xncut,xmcut
 115  format(/,1x,'ntypes=',i3,' kderiv=',i3,/' rcut=',f8.4,
     1    ' frcut=',f8.4,' ncut=',f8.4,' mcut=',f8.4)
	write (6,116) noscr
116	format(' noscr=  ',i1)
	write (6,117) ialloy
117	format(' ialloy=  ',i1)
c
c     take data from namelist or from library if kodes(i)='library'.
c     note that kderiv,deltas and rcut must be defined in namelist data.
c
      do 460 i=1,ntypes
         if (kodes(i) .ne. 'library') go to 460
      rewind kmeamf
      read  (kmeamf,105) header
      write (6,110) header
      read  (kmeamf,105) header
      write (6,110) header
      read  (kmeamf,105) header
      write (6,110) header
      read  (kmeamf,105) header
      write (6,110) header
c
c       this is the only place where the mass is in amu
c       here we convert to eV-psec**2/angstrom**2
c       this is the unit used throughout the program
c       restart assumes this mass unit
c
 410  read  (kmeamf,*) name,lattces(i),zs(i),ielement(i),amamu
      read  (kmeamf,*)   alphas(i,i),(betas(l,i),l=1,4),alats(i),
     1                   esubs(i,i),asubs(i)
      read  (kmeamf,*) (ts(k,i),k=1,4),tempr0,ibarr(i)
         if (name .eq. enames(i)) go to 430
         if (name .ne. 'zz') go to 410
      write (6,420) i,enames(i)
 420  format(1x,'error in meamdat, could not find, i=',i4,2x,a4)
      stop 'error in meamdat, could not find element in library'    
 430  write (6,440)    enames(i),lattces(i),zs(i),ielement(i),amamu
 440  format(1x,a,2x,a,2x,f5.2,i5,f10.4)
      write (6,450) alphas(i,i),(betas(l,i),l=1,4),alats(i),
     1                   esubs(i,i),asubs(i)
 450  format(1x,8f8.4)
	if(rozros(i).eq.0) rozros(i)=tempr0
      write (6,451) (ts(k,i),k=1,4),rozros(i),ibarr(i)
 451  format(1x,5f8.4,i5)
	if (ibarr(i).ge.10) then
	irhos(i)=1
	ibarr(i)=ibarr(i)-10
	elseif(ibarr(i).le.-10) then
	irhos(i)=1
	ibarr(i)=ibarr(i)+10
	else
	irhos(i)=0
	endif
	if(irhos(i).gt.0) print *,'powers of r in density for atom ',i
      amass(i) = conmas*amamu
c      write (6,441) i,amass(i)
c 441  format(1x,'i=',i3,' amass(i)=',e25.16)
 460  continue
      close(unit=kmeamf)
      write (6,120)
 120  format(/,5x,'i',3x,'lattice',6x,'name',5x,'kodes')
      do 130 i=1,ntypes
      write (6,125) i,lattces(i),enames(i),kodes(i)
 125  format(1x,i5,7x,a,8x,a,3x,a)
 130  continue
c
      write (6,150) 
 150  format(/,2x,'i',5x,'esubs(i,j)')
      do 160 i=1,ntypes
      write (6,155) i,(esubs(i,j),j=1,ntypes)
 155  format(1x,i2,5x,4f10.5)
 160  continue
c
      write (6,165) 
 165  format(/,2x,'i',5x,'alphas(i,j)')
      do 175 i=1,ntypes
	if(ibarr(i).eq.0.or. abs(ibarr(i)).eq.4) then
	all(i,i)='sq root'
	elseif(abs(ibarr(i)).eq.1) then
	all(i,i)='ex1'
	elseif(abs(ibarr(i)).eq.2) then
	all(i,i)='ex2'
	elseif(abs(ibarr(i)).eq.3) then
	all(i,i)='ex3'
	else
	write(6,*) 'unknown ibar', ibarr(i), ' for element ',i
	stop
	endif
      write (6,170) i,(alphas(i,j),j=1,ntypes)
 170  format(1x,i2,5x,4f10.5)
 175  continue
c
      write (6,180) 
 180  format(/,2x,'i',4x,'betas(i,j)')
      do 190 i=1,4
      write (6,185) i,(betas(i,j),j=1,ntypes)
 185  format(1x,i2,5x,4f10.5)
 190  continue
c
      write (6,195) 
 195  format(/,2x,'i',8x,'ts(i,j)')
      do 205 i=1,4
      write (6,200) i,(ts(i,j),j=1,ntypes)
 200  format(1x,i2,5x,4f10.5)
 205  continue
c
      write (6,210) 
 210  format(/,2x,'i',4x,'deltas(i,j)')
      do 220 i=1,ntypes
      write (6,215) i,(deltas(i,j),j=1,ntypes)
 215  format(1x,i2,5x,4f10.5)
 220  continue
c
      write (6,297) 
 297  format(/,2x,'i',7x,'repuls(i,j)')
      do 313 i=1,ntypes
      write (6,300) i,(repuls(i,j),j=1,ntypes)
 313  continue
      write (6,298) 
 298  format(/,2x,'i',7x,'attrac(i,j)')
      do 314 i=1,ntypes
      write (6,300) i,(attrac(i,j),j=1,ntypes)
 314  continue
c
	if(nn) itable=1
c
      write (6,135)
 135  format(/,3x,'i',5x,'alats',5x,'asubs',5x,'zs',7x,'z0s',8x,
     1     'rozros')
      do 250 i=1,ntypes
      lattce=lattces(i)
        cmx=cmax(i,i,i)
        cmn=cmin(i,i,i)
c     fcc
         if (lattce .eq. 'fcc') then
	if(zs(i).ne.12) then
	zs(i)=12
	print *,' neighbor inconsistency - neighbors set to ',zs(i)
	endif
	omegas(i)=alats(i)**3/4
	res(i,i)=alats(i)/sqrt(2.)
	a2nn(i)=sqrt(2.)
	if(noscr.eq.0) then
	  arg=1.
          fac=(cmx-arg)/(cmx-cmn)
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
        b2nn(i)=-xfac**4*6./12.
	if (ibarr(i).lt.0) then
	rho0=zs(1)*(1-b2nn(i)*exp(-betas(1,i)*(a2nn(i)-1)))
        arg=0.
	z0s(i)=bar(rho0,arg,ibarr(i),zs(i)*rozros(i),temp1,temp2)
	endif
c
c     bcc
         elseif (lattce .eq. 'bcc') then
	if(zs(i).ne.8) then
	zs(i)=8
	print *,' neighbor inconsistency - neighbors set to ',zs(i)
	endif
	omegas(i)=alats(i)**3/2
	res(i,i)=alats(i)*sqrt(3.)/2
	a2nn(i)=2/sqrt(3.)
        if(noscr.eq.0) then
          arg=2.
          fac=(cmx-arg)/(cmx-cmn)
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
        b2nn(i)=-xfac**4*6./8.
        if (ibarr(i).lt.0) then
	rho0=zs(1)*(1-b2nn(i)*exp(-betas(1,i)*(a2nn(i)-1)))
        arg=0.
	z0s(i)=bar(rho0,arg,ibarr(i),zs(i)*rozros(i),temp1,temp2)
	endif
c
c     dia
        elseif (lattce .eq. 'dia') then
        if(zs(i).ne.4) then
        zs(i)=4
        print *,' neighbor inconsistency - neighbors set to ',zs(i)
        endif
	omegas(i)=alats(i)**3/8
	res(i,i)=alats(i)*sqrt(3.)/4
	a2nn(i)=sqrt(8./3.)
	b2nn(i)=0.
	if(nn.and.cmin(i,i,i).lt.0.500001) then
	write(6,*) 'can not do 2NN MEAM for dia'
	stop
	endif
c
c     hcp (ideal c/a)
        elseif (lattce .eq. 'hcp') then
        if(zs(i).ne.12) then
        zs(i)=12
        print *,' neighbor inconsistency - neighbors set to ',zs(i)
        endif
	if(nn .and. ntypes.gt.1 .and. ibarr(i).lt.0) then
	print *,' can not do 2NN MEAM for alloys and hcp',
     1     '  reference structure'
	stop
	endif
        omegas(i)=alats(i)**3*sqrt(2.)
        res(i,i)=alats(i)
	a2nn(i)=sqrt(2.)
        if(noscr.eq.0) then
          arg=1.
          fac=(cmx-arg)/(cmx-cmn)
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
        b2nn(i)=-xfac**4*6./12.
        if (ibarr(i).lt.0) then
	rho0=zs(1)*(1-b2nn(i)*exp(-betas(1,i)*(a2nn(i)-1)))
	rh3=exp(-betas(1,i)*(a2nn(i)-1))
        arg=(1.0/3.0)*ts(4,i)*rh3**2
	z0s(i)=bar(rho0,arg,ibarr(i),zs(i)*rozros(i),temp1,temp2)
	endif
c
c     dimer
        elseif (lattce .eq. 'dim') then
        if(zs(i).ne.1) then
        zs(i)=1
        print *,' neighbor inconsistency - neighbors set to ',zs(i)
        endif
        omegas(i)=alats(i)**3
        res(i,i)=alats(i)
	b2nn(i)=0.
	a2nn(i)=1.
c
c     error 
	else
 235  write (6,240) lattce
 240  format(1x,'Error in subroutine meamdat, lattce=',a3)
         stop 135
	endif
	if(.not.nn .or. ibarr(i).ge.0) then
        z0s(i)=zbar(ibarr(i),zs(i),lattces(i),ts(2,i),ts(3,i),ts(4,i),
     1    nn,1.)
	endif
        write (6,140) i,alats(i),asubs(i),zs(i),z0s(i),rozros(i)
 140  format(1x,i3,5f10.5)
 250  continue
c
c     compute off diagonal elements
c
	do 310 j=1,ntypes
	do 310 i=2,ntypes
	do 310 k=1,i-1
	cmin(i,j,k)=cmin(k,j,i)
	cmax(i,j,k)=cmax(k,j,i)
 310	continue
      do 260 i=1,ntypes
      do 255 j=1,ntypes
       if(i.eq.j) then
	go to 255
       elseif(i.gt.j) then
	esubs(i,j)=esubs(j,i)
	all(i,j)=all(j,i)
	res(i,j)=res(j,i)
	alphas(i,j)=alphas(j,i)
	repuls(i,j)=repuls(j,i)
	attrac(i,j)=attrac(j,i)
       else
	 if(all(i,j).eq.'L12') then
	i1=min(i,j)
	i2=max(i,j)
	if(esubs(i,j).eq.0.) then
	esubs(i,j)=(3*esubs(i1,i1)+esubs(i2,i2))/4. - deltas(i,j)
	endif
        if(alphas(i,j).eq.0.) then
	alphas(i,j)=(3*alphas(i1,i1)+alphas(i2,i2))/4.
	endif
        if(res(i,j).eq.0.) then
	res(i,j)=(res(i,i)+res(j,j))/2.
	endif
	 else
	if(esubs(i,j).eq.0.) then
	esubs(i,j)=(esubs(i,i)+esubs(j,j))/2. - deltas(i,j)
	endif
        if(alphas(i,j).eq.0.) then
        alphas(i,j)=(alphas(i,i)+alphas(j,j))/2.
	endif
        if(res(i,j).eq.0.) then
        res(i,j)=(res(i,i)+res(j,j))/2.
	endif
	 endif
         if (rz0.ne.1.0) then
           res(i,j)=rz0*res(i,j)
         endif
         if (alpha0.ne.1.0) then
           alphas(i,j)=alpha0*alphas(i,j)
         endif
        endif
 255  continue
 260  continue
c	compute default repulsive term
cdo 1311 i=1,ntypes
cdo 1311 j=1,ntypes
cif(repuls(i,j).lt.0.) then
crepuls(i,j)=hartre*bohr*ielement(i)*ielement(j)/alphas(i,j)**3
c    1   /esubs(i,j)/res(i,j)*exp(-alphas(i,j))*(zs(i)+zs(j))/4
cendif
c311	continue
c
c     edit complete arrays
c
      write (6,261)
 261  format(1x,'Computed arrays')
      write (6,265) 
 265  format(/,2x,'i',5x,'esubs(i,j)')
      do 275 i=1,ntypes
      write (6,270) i,(esubs(i,j),j=1,ntypes)
 270  format(1x,i2,5x,4f10.5)
 275  continue
c
      write (6,280) 
 280  format(/,2x,'i',5x,'alphas(i,j)')
      do 290 i=1,ntypes
      write (6,285) i,(alphas(i,j),j=1,ntypes)
 285  format(1x,i2,5x,4f10.5)
 290  continue
c
      write (6,295) 
 295  format(/,2x,'i',7x,'res(i,j)')
      do 305 i=1,ntypes
      write (6,300) i,(res(i,j),j=1,ntypes)
 300  format(1x,i2,5x,4f10.5)
 305  continue
c
      write (6,297) 
      do 311 i=1,ntypes
      write (6,300) i,(repuls(i,j),j=1,ntypes)
 311  continue
      write (6,298) 
      do 312 i=1,ntypes
      write (6,300) i,(attrac(i,j),j=1,ntypes)
 312  continue
c
      write (6,296)
 296  format(/,2x,'i',7x,'all(i,j)')
      do 306 i=1,ntypes
      write (6,301) i,(all(i,j),j=1,ntypes)
 301  format(1x,i2,5x,4a10)
 306  continue 
c
      alat=alats(1)
c
	write(6,*) ' For screening:'
	do 307 i=1,ntypes
	write(6,*) ' type ',i
	write(6,*) ' cmax'
	do 308 j=1,ntypes
	write(6,9331) (cmax(j,i,k),k=1,ntypes)
  308	continue
	write(6,*) ' cmin'
	do 309 j=1,ntypes
	write(6,9331) (cmin(j,i,k),k=1,ntypes)
  309	continue
9331	format(10f8.2)
  307	continue
400	continue
      write (6,330)
 330  format(/,1x,'End of MEAM data')
c
	if(itable.gt.0) then
	write(6,*) 'setting up table for phi'
	do 407 i=1,ntypes
	if(nn) then
	nmax=10
	write(6,*) 'calculating 2NN MEAM with a2nn,b2nn=',a2nn(i),b2nn(i)
	else
	nmax=0
	endif
	do 407 k=1,ntable
	r=k*rcut/ntable
	phitab(k,i,i)=phif(r,i,i)
	do 408 n=1,nmax
cprint *,' phi ',r,phitab(k,i,i),n
	phitab(k,i,i)=phitab(k,i,i)+b2nn(i)**n*phif(r*a2nn(i)**n,i,i)
408	continue
407	continue
	itable=-1
	do 406 i=1,ntypes
	do 406 j=i+1,ntypes
	print *, 'calculating mixed phi',i,j
        do 406 k=1,ntable
        r=k*rcut/ntable
	phitab(k,i,j)=phif(r,i,j)
	phitab(k,j,i)=phitab(k,i,j)
406	continue
	itable=-2
	endif
      return
	end
c************************************************************************
c
c  this routine outputs the current configuration to file 20
c  as well as the instantaneous values of the temperature and
c  stress tensor.
c
      subroutine output(t,iclfor)
      include 'param.inc'
      include 'dyn.inc'
      include 'meam.inc'
      include 'wrkspc.inc'
	common/head/ header
	character *80 header
      data de/0.0/,eold/0.0/,ifirst/0/,ipcnt/-1/
c
c
c  compute the physical (unscaled) coordinates
c
	call uscale
c
c  compute the energy, temperature and pressure
c
      call calce(ekin,pot,etot,temp,press,iclfor,t)
c
c  determine the maximun change in total energy
c
      if(ifirst.eq.0)then
        ifirst=1
        eold=etot
        detops = 0.0
      else
        de=abs(etot-eold)
        detops = max(detops,de)
        eold=etot
      endif
c
c  print out energy, temperature and pressure
c
      ipcnt = mod(ipcnt + 1,ipitera)
      if(ipcnt.eq.0)then
c   write the restart file each output
      rewind(21)
      write(21,9501) header
9501  format(a)
      write(21,9502) natoms,ntypes,t,beta
9502  format(2i10,e15.8,f10.5)
      write(21,9503) (perub(i),i=1,3),(perlb(i),i=1,3)
9503  format(3e25.16)
      write(21,9504) (amass(i),ielement(i),i=1,ntypes)
9504  format(e25.16,i10)
      write(21,9505) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
9505  format(3e25.16/3e25.16/i10)
      write(21,9503) (bndvel(i),i=1,3)
c output the averages over the output period
	call average(2,t)
	endif
      call flush(21)
      call flush(6)
c
c  write to file 20
c
      if (iconf.eq.1) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,3),itype(j),j=1,natoms)
      elseif (iconf.eq.2) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
      elseif (abs(iconf).eq.3) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,3),itype(j),j=1,natoms)
         if(iconf.gt.0) then
	  write(20) (((avslocal(i,j,k),i=1,j),j=1,3),k=1,natoms)
	 else
          write(20) (((avslocal(i,j,k),i=1,3),j=1,3),k=1,natoms)
	 endif
      elseif (iconf.eq.4) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,6),itype(j),e(j),j=1,natoms)
      elseif (iconf.eq.5) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,3),itype(j),j=1,natoms)
         write(20) (rho(k),c(k),(rhsq(j,k),j=1,3),k=1,natoms)
      elseif (abs(iconf).eq.6) then
         write(20) t,temp,press,(perlen(i),i=1,3)
         write(20) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
         if(iconf.gt.0) then
          write(20) (((avslocal(i,j,k),i=1,j),j=1,3),k=1,natoms)
         else
          write(20) (((avslocal(i,j,k),i=1,3),j=1,3),k=1,natoms)
         endif
       end if
c reinitialize averages
	call average(0,t)
      return
      end
