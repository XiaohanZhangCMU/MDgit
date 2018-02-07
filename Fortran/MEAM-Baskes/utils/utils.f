c utils.f
c$Author: baskes $
c$Date: 2004/10/08 18:13:12 $
c$Source: /home/baskes/programs/EAM/dyn/utils/RCS/utils.f,v $
c$Revision: 1.33 $
c$Log: utils.f,v $
cRevision 1.33  2004/10/08 18:13:12  baskes
cin verlet replaced dyn.inc with utils.inc
c
cRevision 1.32  2003/12/05 22:03:12  baskes
cupdating
c
cRevision 1.31  2003/11/10 18:42:49  baskes
cchanged rcutsq to rctsqn
cmakes neighbor region larger
c
cRevision 1.30  2003/08/20 17:17:57  baskes
cnew .inc
c
cRevision 1.29  2003/07/01 03:37:37  baskes
cgetting ready to change dyn.inc
c
cRevision 1.28  2002/11/04 18:29:26  baskes
cchanged print
c
cRevision 1.27  2002/08/25 20:41:41  baskes
cfixed print formats
cremoved getcell print
c
cRevision 1.26  2002/05/08 16:25:53  baskes
cfixed iper
c
cRevision 1.25  2002/05/07 21:38:13  baskes
cas of 5/2/2002
c
cRevision 1.24  2002/04/30 20:49:49  baskes
caveraging fixed again
c
cRevision 1.23  2002/04/30 18:09:17  baskes
cfix averaging start
c
cRevision 1.22  2002/04/29 18:43:26  baskes
cadded RMS for force
c
cRevision 1.21  2002/04/07 00:47:24  baskes
cmove restart to utils
cadd uniform rate of periodic vector change
c
cRevision 1.20  2001/12/04 22:39:34  baskes
cCm fixed
c
cRevision 1.19  2001/12/04 18:53:22  baskes
cadd icm zero CM velocity
c
cRevision 1.18  2001/06/12 15:38:31  baskes
cadded acf
c
cRevision 1.17  2001/01/10 21:18:28  baskes
cfixed short period for nmeth=-2
c
cRevision 1.16  2001/01/08 23:19:10  baskes
cfixed more dimensions for timing
c
cRevision 1.15  2001/01/08 23:12:25  baskes
cfixed dimension for timing
c
	subroutine getnei(r,ncell,natc,natcmax,icmax,icmay,icmaz,
     1     nei,neigh,neimx,cmono)
c   get the possible neighbors of atom at r
        implicit real*8 (a-h,o-z)
        common /cell/ xmin(3),xcell(3),icell(3)
	dimension r(3),ijk(3),ncell(natcmax,icmax,icmay,icmaz),
     1     natc(icmax,icmay,icmaz),neigh(neimx)
	nei=0
	do 1 k=1,3
1	ijk(k)=(r(k)-xmin(k))/xcell(k)+1
	nmono=abs(cmono)/xcell(1)+0.5
	if(cmono.lt.0) nmono=-nmono
	if(icell(3).eq.1) then
	kl=0
	ku=0
	elseif(icell(3).eq.2) then
	kl=0
	ku=1
	else
	kl=1
	ku=1
	endif
	do 21 kx=ijk(3)-kl,ijk(3)+ku
	kxp=mod(kx-1+icell(3),icell(3))+1
        if(icell(2).eq.1) then
        jl=0
        ju=0
        elseif(icell(2).eq.2) then
        jl=0
        ju=1
        else
        jl=1
        ju=1
        endif
        do 22 jx=ijk(2)-jl,ijk(2)+ju
	jxp=mod(jx-1+icell(2),icell(2))+1
        if(icell(1).eq.1) then
        il=0
        iu=0
        elseif(icell(1).eq.2) then
        il=0
        iu=1
        else
        il=1
        iu=1
        endif
	ixmin=ijk(1)-il
	ixmax=ijk(1)+iu
	if(abs(cmono).gt.1e-3.and.icell(1).gt.3.and.
     1        abs(ijk(3)-kxp).gt.1) then
	 if(kxp.eq.1) then
	  if(icell(1).gt.4) then
	  ixmin=ijk(1)-2-nmono
	  ixmax=ijk(1)+2-nmono
	  else
	  ixmin=1
	  ixmax=icell(1)
	  endif
	 elseif(kxp.eq.icell(3)) then
	  if(icell(1).gt.4) then
	  ixmin=ijk(1)-2+nmono
	  ixmax=ijk(1)+2+nmono
	  else
	  ixmin=1
	  ixmax=icell(1)
	  endif
	 endif
	endif
        do 23 ix=ixmin,ixmax
	ixp=mod(ix-1+icell(1),icell(1))+1
	do 2 j=1,natc(ixp,jxp,kxp)
	nei=nei+1
	if(nei.gt.neimx) go to 100
2	neigh(nei)=ncell(j,ixp,jxp,kxp)
23	continue
22	continue
21	continue
	return
100	write(6,*) ' in getnei, too many atoms, increase neimx'
	stop
	end
	subroutine getcell(r,natoms,ncell,natc,natcmax,icmax,icmay,
     1      icmaz,rc)
c   put the atoms in cells
        implicit real*8 (a-h,o-z)
	common /cell/ xmin(3),xcell(3),icell(3)
	dimension xmax(3),ijk(3)
	dimension r(6,natoms),ncell(natcmax,icmax,icmay,icmaz),
     1     natc(icmax,icmay,icmaz)
	do 20 i=1,icmax
	do 20 j=1,icmay
	do 20 k=1,icmaz
20	natc(i,j,k)=0
	do 1 k=1,3
	xmax(k)=-1e10
 	xmin(k)=1e10
	do 1 i=1,natoms
	xmax(k)=max(xmax(k),r(k,i))
1	xmin(k)=min(xmin(k),r(k,i))
	do 2 k=1,3
	if(k.eq.1) icell(k)=min((xmax(k)-xmin(k))/rc,icmax*1.0D0)
	if(k.eq.2) icell(k)=min((xmax(k)-xmin(k))/rc,icmay*1.0D0)
	if(k.eq.3) icell(k)=min((xmax(k)-xmin(k))/rc,icmaz*1.0D0)
	icell(k)=max(icell(k),1)
2	xcell(k)=(xmax(k)-xmin(k))/icell(k)+rc/1000
cwrite(6,1000) xcell,icell
1000	format(' creating cells ',3f10.2,3i5)
	do 3 i=1,natoms
	do 4 k=1,3
4	ijk(k)=(r(k,i)-xmin(k))/xcell(k)+1
	natc(ijk(1),ijk(2),ijk(3))= natc(ijk(1),ijk(2),ijk(3))+1
	if(natc(ijk(1),ijk(2),ijk(3)).gt.natcmax) go to 100
3	ncell(natc(ijk(1),ijk(2),ijk(3)), ijk(1),ijk(2),ijk(3))=i
	return
100	write(6,*) ' in getcell, too many atoms, increase natcmax'
	stop
	end
	function zbl(r,z1,z2)
      implicit real*8(a-h,o-z)
	integer z1,z2
	dimension c(4),d(4)
	data c/0.028171,0.28022,0.50986,0.18175/
	data d/0.20162,0.40290,0.94229,3.1998/
	data azero/0.4685/
	data cc/14.3997/
c   azer0=(9pi^2/128)^1/3 (0.529) Angstroms
	a=azero/(z1**0.23+z2**0.23)
	zbl=0.0
	x=r/a
	do 1 i=1,4
1	zbl=zbl+c(i)*exp(-d(i)*x)
	zbl=zbl*z1*z2/r*cc
	return
	end
      function rgauss()
      implicit real*8(a-h,o-z)
      data twopi/6.283185308/
   10 continue
      a1 = ranl()
      if (a1.eq.0.) goto 10
      a2 = ranl()
      rgauss = sqrt(-2.*log(a1))*cos(twopi*a2)
      return
      end
c
c***********************************************************************
c
c  random number functions
c
c**********************************************************************
      subroutine lranst(iseed)
      implicit real*8 (a-h,o-z)
	integer time
      if (iseed.le.0) then
         call srand(time())
      else
         call srand(iseed)
      endif
      return
      end
      function ranl()
      implicit real*8(a-h,o-z)
      ranl = rand(0)
      return
      end
c***********************************************************************
c
c  vector manipulation routines
c
      function sdot(n,x,incx,y,incy)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n)
      tot = 0.0
      do 100 i = 1,n
      ix = 1 + (i-1)*incx
      iy = 1 + (i-1)*incy
      tot = tot + x(ix)*y(iy)
100   continue
      sdot = tot
      return
      end
      subroutine saxpy(n,a,x,incx,y,incy)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n)
      do 100 i = 1,n
      ix = 1 + (i-1)*incx
      iy = 1 + (i-1)*incy
      y(iy) = a*x(ix) + y(iy)
100   continue
      return
      end
      function snrm2(n,x,inc)
        implicit real*8 (a-h,o-z)
      dimension x(1)
      ilast = n*inc
      sum = 0.
      do 100 i = 1,ilast,inc
         sum = sum + x(i)**2
 100     continue
      snrm2 = sqrt(sum)
      return
      end
      subroutine scopy(n,x,incx,y,incy)
        implicit real*8 (a-h,o-z)
      dimension x(1),y(1)
      do 100 i = 0,n-1
         jx = 1 + i*incx
         jy = 1 + i*incy
         y(jy) = x(jx)
 100     continue
      return
      end
c*******************************************************************
      subroutine m01aaf(a,ii,jj,ip,ia,ifail)
c*****************************************************************
c
c  subroutine qsort is a quick sort algorithm designed to look
c  like the nag routine m01aaf
c
c  a is the input array to be sorted
c   (a is unchanged after the sort)
c  ii,jj are the bounds of the part of the array to be sorted
c  ip(i) = position of a(i) in the sorted list
c  ia(i) = index of the ith element in the sorted list
c    (i.e. the elements in order are a(ia(i))
c  ifail is a error message
c
c  the algorithm is taken from  Richard C. Singleton
c   Communications of the ACM, vol. 12 Number 3, March 1969, p. 185
c   (Algorithm 347)
c
        implicit real*8 (a-h,o-z)
      dimension a(jj),ia(jj),ip(jj),iu(16),il(16)
      ifail = 0
      if (ii.gt.jj) ifail = 1
      if (ii.lt.1) ifail = 2
      if (jj.lt.1) ifail = 3
      if (ifail.ne.0) return
      do 2 i = ii,jj
2     ia(i) = i
      m = 1
      i = ii
      j = jj
5     if (i.ge.j) goto 70
10    k = i
      ij = (i+j)/2
      t = a(ij)
      it = ia(ij)
      if (a(i).le.t) goto 20
      a(ij) = a(i)
      a(i) = t
      t = a(ij)
      ia(ij) = ia(i)
      ia(i) = it
      it = ia(ij)
20    l = j
      if (a(j) .ge. t) goto 40
      a(ij) = a(j)
      a(j) = t
      t = a(ij)
      ia(ij) = ia(j)
      ia(j) = it
      it = ia(ij)
      if (a(i) .le. t) goto 40
      a(ij) = a(i)
      a(i) = t
      t = a(ij)
      ia(ij) = ia(i)
      ia(i) = it
      it = ia(ij)
      goto 40
30    a(l) = a(k)
      a(k) = tt
      ia(l) = ia(k)
      ia(k) = itt
40    l = l - 1
      if (a(l) .gt. t) goto 40
      tt = a(l)
      itt = ia(l)
50    k = k + 1
      if (a(k) .lt. t) goto 50
      if (k .le. l) goto 30
      if (l-i .le. j-k) goto 60
      il(m) = i
      iu(m) = l
      i = k
      m = m + 1
      goto 80
60    il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      goto 80
70    m = m - 1
      if (m .eq. 0) goto 1000
      i = il(m)
      j = iu(m)
80    if (j-i .ge. ii) goto 10
      if (i .eq. ii) goto 5
      i = i - 1
90    i = i + 1
      if (i .eq. j) goto 70
      t = a(i+1)
      it = ia(i+1)
      if (a(i) .le. t) goto 90
      k = i
100   a(k+1) = a(k)
      ia(k+1) = ia(k)
      k = k - 1
      if (t .lt. a(k)) goto 100
      a(k+1) = t
      ia(k+1) = it
      goto 90
1000  continue
      do 1100 i = ii,jj
1100  ip(ia(i)) = i
c
c  now return the data to the original order
c
      do 2000 k = ii,jj
      if (ia(k).lt.0) goto 2000
      k1 = k
      ik1 = ia(k1)
      temp = a(k1)
2100  continue
      ia(k1) = -ia(k1)
      store = a(ik1)
      a(ik1) = temp
      temp = store
      k1 = ik1
      ik1 = ia(k1)
      if (ik1.gt.0) goto 2100
2000  continue
      do 2200 k = ii,jj
2200  ia(k) = -ia(k)
      return
      end
c************************************************************************
c
c   the following four functions provide timing information
c
c***********************************************************************
c**********************************************************************
c
c  timing routines
c
      subroutine initsec()
      implicit real*8(a-h,o-z)
      real*4 tarray(2),etime,start
      common /seccom/ start
      start = etime(tarray)
      return
      end
      function seconds()
      implicit real*8(a-h,o-z)
      real*4 tarray(2),etime,start,t
      common /seccom/ start
      t = etime(tarray) - start
      seconds = t
      return
      end
c
c  return the remaining time left for this job
c  not meaning for this version, just return a large number
c
      subroutine trmain(timlft)
        implicit real*8 (a-h,o-z)
      timlft = 1.e10
      return
      end
c
      subroutine mytime(timestr)
        implicit real*8 (a-h,o-z)
      character*24 timestr, ctime
      integer*4 time
      timestr = ctime(time())
      return
      end
	subroutine rdx3(x1,x2,per12,perlen,cmono,rd,rsq)
c   calculate the closest distances between two points, x1 and x2
c   rd contains the xyz distances and rsq the distance squared
      implicit real*8(a-h,o-z)
	dimension x1(3),x2(3),per12(3),perlen(3),rd(3)
	kz=0
	if(cmono.gt.10000000.0) then
	   if(x1(3).ge.perlen(3)) then
	   kz=kz+1
	   elseif (x1(3).lt.0.) then
	   kz=kz-1
	   endif
           if(x2(3).ge.perlen(3)) then
           kz=kz-1
           elseif (x2(3).lt.0.) then
           kz=kz+1
           endif
	endif
      rd(3) = x1(3)-x2(3)
      if (rd(3) .lt. -per12(3)) then
        rd(3) = rd(3) + perlen(3)
	kz=kz+1
      else if (rd(3) .gt. per12(3)) then
        rd(3) = rd(3) - perlen(3)
	kz=kz-1
	endif
      rd(1) = x1(1)-x2(1) + kz*cmono
      if (rd(1) .lt. -per12(1)) then
        rd(1) = rd(1) + perlen(1)
      else if (rd(1) .gt. per12(1)) then
        rd(1) = rd(1) - perlen(1)
	endif
      rd(2) = x1(2)-x2(2)
      if (rd(2) .lt. -per12(2)) then
        rd(2) = rd(2) + perlen(2)
      else if (rd(2) .gt. per12(2)) then
        rd(2) = rd(2) - perlen(2)
	endif
      rsq = rd(1)*rd(1)+rd(2)*rd(2)+rd(3)*rd(3)
	return
	end
c***************************************************************
c
      subroutine latgen(contin)
      include 'param.inc'
      include 'utils.inc'
      logical contin
      character*8 lattype
      dimension avec(3),bvec(3),cvec(3),xrot(3),yrot(3),zrot(3)
      dimension xold(3),yold(3),zold(3),avecp(3),bvecp(3),cvecp(3)
      dimension roter(3,3),aperub(3),aperlb(3)
      dimension axbnd(2),aybnd(2),azbnd(2)
      dimension perubs(3),perlbs(3)
      dimension rcell(3,10),ccell(nelmax,10)
c
c  default values for namelist
c
      data alat/3.52/
      data lattype/'none'/,xrot/3*0./,yrot/3*0./,zrot/3*0./
      data xold/1.,0.,0./,yold/0.,1.,0./,zold/0.,0.,1./
      data avec/3*0.0/,bvec/3*0.0/,cvec/3*0.0/
      data perub/3*9999.0/,perlb/3*-9999.0/
      data xbound/-9999.,9999./,ybound/-9999.,9999./
      data zbound/-9999.,9999./
      data aperlb/3*-9999./,aperub/3*9999./
      data axbnd/-9999.,9999./,aybnd/-9999.,9999./,azbnd/-9999.,9999./
      data ncell/1/
      data rcell/30*0./
      namelist /latcard/ lattype,avec,bvec,cvec,alat,xrot,yrot,zrot,
     1perub,perlb,xbound,ybound,zbound,aperub,aperlb,axbnd,aybnd,
     2azbnd,ncell,rcell,ccell
c
c	only generate orthogonal periodic vectors
	betar=90.
      do 1 i = 1,nelmax
         do 1 j = 1,10
            ccell(i,j) = 0.
 1         continue
      ccell(1,1) = 1.
c
c  start the time clock at 0 for a generated configuration
c
      t0 = 0.0
c
c  read the lattice parameters
c
      if (contin) then
         perlb(1) = perlbs(1)
         perub(1) = perubs(1)
         perlb(2) = perlbs(2)
         perub(2) = perubs(2)
         perlb(3) = perlbs(3)
         perub(3) = perubs(3)
       else
         read(5,latcard)
       endif
c
c  if no periodic bounds set, assume a 4 period cube and scale the
c  periodic lengths from units of lattice constants to real units
c
      do 3 i = 1,3
      if (perlb(i).le.-9998..and.aperlb(i).le.-9998.) aperlb(i) = -2.0
      if (perub(i).ge.9998..and.aperub(i).ge.9998.) aperub(i) = 2.0
      if (perlb(i).le.-9998.) perlb(i) = alat*aperlb(i)
      if (perub(i).ge.9998.) perub(i) = alat*aperub(i)
3     continue
c
c  save the input periodic lenths for use in a continue
c
         perlbs(1) = perlb(1)
         perubs(1) = perub(1)
         perlbs(2) = perlb(2)
         perubs(2) = perub(2)
         perlbs(3) = perlb(3)
         perubs(3) = perub(3)
c
c  set the x,y,and z bounds to fill the periodic bounds if they have not
c  been set by the user
c
      if (axbnd(1).gt.-9998.) xbound(1) = alat*axbnd(1)
      if (aybnd(1).gt.-9998.) ybound(1) = alat*aybnd(1)
      if (azbnd(1).gt.-9998.) zbound(1) = alat*azbnd(1)
      if (axbnd(2).lt.9998.) xbound(2) = alat*axbnd(2)
      if (aybnd(2).lt.9998.) ybound(2) = alat*aybnd(2)
      if (azbnd(2).lt.9998.) zbound(2) = alat*azbnd(2)
      if (xbound(1).le.-9998.) xbound(1) = perlb(1) + 0.01
      if (xbound(2).ge.9998.) xbound(2) = perub(1) + 0.01
      if (ybound(1).le.-9998.) ybound(1) = perlb(2) + 0.01
      if (ybound(2).ge.9998.) ybound(2) = perub(2) + 0.01
      if (zbound(1).le.-9998.) zbound(1) = perlb(3) + 0.01
      if (zbound(2).ge.9998.) zbound(2) = perub(3) + 0.01
c
c  see if primitive lattice vectors have been defined
c  if not, find them using lattype
c      if lattype is not defined from namelist, then use latty of
c      first atom type
      primdef = (avec(1)**2 + avec(2)**2 + avec(3)**2)*
     1(bvec(1)**2 + bvec(2)**2 + bvec(3)**2)*
     2(cvec(1)**2 + cvec(2)**2 + cvec(3)**2)
      if(primdef.ne.0.0)then
         lattype = 'none'
      else
         if(lattype.eq.'none')lattype = 'fcc'
         call storelat(lattype,avec,bvec,cvec)
      endif
      if(lattype.ne.'none') then
         write(6,9115)lattype
       end if
9115  format(//' ******  generating ',a8,' lattice ')
      write(6,9116)avec,bvec,cvec
9116  format(' primitive lattice vectors ',3f5.2,3x,3f5.2,3x,3f5.2)
      xnorm = snrm2(3,xrot,1)
      ynorm = snrm2(3,yrot,1)
      znorm = snrm2(3,zrot,1)
c
c       if two defined, find missing one from them
c
      if((xnorm.ne.0.0).and.(ynorm.ne.0.0).and.(znorm.eq.0.0))
     1      call cross(zrot,xrot,yrot)
      if((ynorm.ne.0.0).and.(znorm.ne.0.0).and.(xnorm.eq.0.0))
     1      call cross(xrot,yrot,zrot)
      if((znorm.ne.0.0).and.(xnorm.ne.0.0).and.(ynorm.eq.0.0))
     1      call cross(yrot,zrot,xrot)
c
c       if none or only one defined, set equal to no rotation
c
      if((xnorm*ynorm.eq.0.0).and.(xnorm*znorm.eq.0.0).
     1and.(ynorm*znorm.eq.0.0))then
            call scopy(3,xold,1,xrot,1)
            call scopy(3,yold,1,yrot,1)
            call scopy(3,zold,1,zrot,1)
      endif
c
c       normalize the unit vectors and define rotation matrix
c
      xnorm = snrm2(3,xrot,1)
      ynorm = snrm2(3,yrot,1)
      znorm = snrm2(3,zrot,1)
      do 5 i=1,3
      xrot(i) = xrot(i)/xnorm
      yrot(i) = yrot(i)/ynorm
      zrot(i) = zrot(i)/znorm
5     continue
      do 10 i=1,3
      do 10 j=1,3
      roter(i,j)= xold(i)*xrot(j) + yold(i)*yrot(j) + zold(i)*zrot(j)
10    continue
c
c       rotate avec,bvec,cvec
c
      do 20 i=1,3
      avecp(i)=0.0
      bvecp(i)=0.0
      cvecp(i)=0.0
      do 20 j=1,3
      avecp(i) = avecp(i) + roter(i,j)*avec(j)
      bvecp(i) = bvecp(i) + roter(i,j)*bvec(j)
 20   cvecp(i) = cvecp(i) + roter(i,j)*cvec(j)
      do 30 i=1,3
      avec(i) = avecp(i)
      bvec(i) = bvecp(i)
 30   cvec(i) = cvecp(i)
      write(6,9216)avec,bvec,cvec
9216  format(' rotated lattice vectors   ',3f5.2,3x,3f5.2,3x,3f5.2)
c
c  compute perlen
c
      do 50 i = 1,3
50     perlen(i) = perub(i) - perlb(i)
c
c       half-lattice constant (angstroms)
c
      ahlc = 0.5*alat
      write(6,9030)ahlc
 9030 format(' half-lattice constant = ',g15.5)
      write(6,9031)
 9031 format('                ',15x,'x',14x,'y',14x,'z')
      write(6,9032)xbound(1),ybound(1),zbound(1)
 9032 format(' lower cell bound      ',3g15.5)
      write(6,9033)xbound(2),ybound(2),zbound(2)
 9033 format(' upper cell bound      ',3g15.5)
      write(6,9034)perlb
 9034 format(' lower periodic bound  ',3g15.5)
      write(6,9035)perub
 9035 format(' upper periodic bound  ',3g15.5)
      write(6,9036)perlen
 9036 format(' length        ',8x,3g15.5)
      xlen = max(abs(avec(1)),abs(bvec(1)),abs(cvec(1)))*ahlc
      ylen = max(abs(avec(2)),abs(bvec(2)),abs(cvec(2)))*ahlc
      zlen = max(abs(avec(3)),abs(bvec(3)),abs(cvec(3)))*ahlc
      nsidex = int((xbound(2)-xbound(1))/xlen) + 1
      nsidey = int((ybound(2)-ybound(1))/ylen) + 1
      nsidex = max(nsidex,nsidey)
      nsidez = int((zbound(2)-zbound(1))/zlen) + 1
      nsidex = max(nsidex,nsidez)
      nsidex = 2*nsidex
      ncentx = int(0.6 + 0.5*float(nsidex))
      do 190 l = 1,ncell
      do 190 m = 2,ntypes
      ccell(m,l) = min(1.,ccell(m-1,l) + ccell(m,l))
190   continue
      natoms = 0
      do 200 i=1,nsidex
      ii = i - ncentx
      do 200 j=1,nsidex
      jj = j - ncentx
      do 200 k=1,nsidex
      kk = k - ncentx
      xc = ahlc*(ii*avec(1) + jj*bvec(1) + kk*cvec(1))
      yc = ahlc*(ii*avec(2) + jj*bvec(2) + kk*cvec(2))
      zc = ahlc*(ii*avec(3) + jj*bvec(3) + kk*cvec(3))
      do 200 l = 1,ncell
      xt = xc + rcell(1,l)
      yt = yc + rcell(2,l)
      zt = zc + rcell(3,l)
      if(xt.ge.xbound(2))go to 200
      if(xt.lt.xbound(1))go to 200
      if(yt.ge.ybound(2))go to 200
      if(yt.lt.ybound(1))go to 200
      if(zt.ge.zbound(2))go to 200
      if(zt.lt.zbound(1))go to 200
      tst = ranl()
      it = 0
      do 201 m=ntypes,1,-1
201   if (tst.lt.ccell(m,l)) it = m
      if (it.le.0.or.it.gt.ntypes) goto 200
      natoms = natoms + 1
      if(natoms.gt.natmax)go to 300
      itype(natoms) = it
      rv(1,natoms) = xt
      rv(2,natoms) = yt
      rv(3,natoms) = zt
 200  continue
      return
c
 300  write(6,9601)natmax
 9601 format('   number of atoms generated exceeds maximum dimension:'
     1 ,i10)
      return
      end
c*******************************************************************
c
c  dummy program to read the latcard namelist entry
c
      subroutine latdum(contin)
        implicit real*8 (a-h,o-z)
      logical contin
      dimension avec(3),bvec(3),cvec(3),xrot(3),yrot(3),zrot(3)
      dimension perub(3),perlb(3),xbound(2),ybound(2),zbound(2)
      dimension aperub(3),aperlb(3),axbnd(2),aybnd(2),azbnd(2)
      namelist /latcard/ lattype,avec,bvec,cvec,alat,xrot,yrot,zrot,
     1perub,perlb,xbound,ybound,zbound,aperub,aperlb,axbnd,aybnd,
     2azbnd
c
c  read the namelist latcard.  note that the values are discarded.
c
      if (contin) return
      read(5,latcard)
      return
      end
      subroutine cross(z,x,y)
        implicit real*8 (a-h,o-z)
      dimension z(3),x(3),y(3)
      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
      return
      end
c***************************************************************
c
      subroutine storelat(latty,avec,bvec,cvec)
        implicit real*8 (a-h,o-z)
      character*8 lstored(10),latty
      dimension avecs(3,10),bvecs(3,10),cvecs(3,10)
      dimension avec(3),bvec(3),cvec(3)
      data lstored/'fcc     ','bcc     ','sc      ',7*'        '/
      data avecs/1.,1.,0., 1.,1.,-1., 2.,0.,0., 21*0.0/
      data bvecs/0.,1.,1., -1.,1.,1., 0.,2.,0., 21*0.0/
      data cvecs/1.,0.,1., 1.,-1.,1., 0.,0.,2., 21*0.0/
c       primitive lattice vectors given in terms of half-lattice constants
      imatch = 0
      do 10 i=1,10
      if(latty.ne.lstored(i))go to 10
        imatch = 1
        do 5 j=1,3
        avec(j) = avecs(j,i)
        bvec(j) = bvecs(j,i)
        cvec(j) = cvecs(j,i)
 5      continue
 10   continue
      if(imatch.eq.1)return
      write(6,15)latty
 15   format('   could not find lattice type ',a8,
     1'. assume default type.')
c  assume default
      i = 1
      latty = lstored(i)
      do 20 j=1,3
      avec(j) = avecs(j,i)
      bvec(j) = bvecs(j,i)
      cvec(j) = cvecs(j,i)
 20   continue
      return
      end
c*****************************************************************************
c
c  the subroutine readef reads in the defcards
c
      subroutine readef(ndef)
      include 'param.inc'
      include 'utils.inc'
      integer type,oldtype,newtype
      dimension num(2),pos(3),delpos(3),vel(3)
      namelist /defcard/ num,type,oldtype,pos,delpos,vel,newtype,
     1                   xmin,xmax,ymin,ymax,zmin,zmax
c
c  initialize the count of the number of real defcards
c
      ndef = 0
100   continue
c
c   set all the values of the namelist variables to their default values
c
      newtype = 99
      oldtype = 99
      type = 99
      num(1) = -1
      num(2) = -1
      do 10 i = 1,3
      pos(i) = 9999.
      delpos(i) = 0.
      vel(i) = 9999.
10    continue
      xmin = -9999.
      xmax = 9999.
      ymin = -9999.
      ymax = 9999
      zmin = -9999.
      zmax = 9999.
c
c  read in the next defcard and check to see if anyting is changed
c
      read(5,defcard)
      if (newtype.ne.99) type=newtype
      if (num(1).ne.-1) goto 1000
      if (num(2).ne.-1) goto 1000
      if (type.ne.99) goto 1000
      if (oldtype.ne.99) goto 1000
      if (xmin.ne.-9999..or.xmax.ne.9999.) goto 1000
      if (ymin.ne.-9999..or.ymax.ne.9999.) goto 1000
      if (zmin.ne.-9999..or.zmax.ne.9999.) goto 1000
      if (delpos(1).ne.0.) goto 1000
      if (delpos(2).ne.0.) goto 1000
      if (delpos(3).ne.0.) goto 1000
c
c  nothing changed so return
c
      return
c
1000  continue
c
c  there is a nontrivial defcard
c
      ndef = ndef + 1
      write(6,defcard)
c
c  branch depending on whether the defects are specified by number or position
c
      if (num(1).ne.-1) goto 1500
c
c  specified by position
c
      do 1100 i = 1,natoms
      if (rv(1,i).lt.xmin.or.rv(1,i).gt.xmax) goto 1100
      if (rv(2,i).lt.ymin.or.rv(2,i).gt.ymax) goto 1100
      if (rv(3,i).lt.zmin.or.rv(3,i).gt.zmax) goto 1100
      if (oldtype.ne.99.and.itype(i).ne.oldtype) goto 1100
      if (type.ne.99) itype(i) = type
      do 1110 j = 1,3
      if (pos(j).ne.9999.) then
         rv(j,i) = pos(j)
       else
         rv(j,i) = rv(j,i) + delpos(j)
       endif
      if (vel(j).ne.9999.) rv(3+j,i) = vel(j)
1110  continue
1100  continue
      goto 100
c
c  specified by number
c
1500  continue
      if (num(2).eq.-1) num(2) = num(1)
c
c  creating a new atom
c
      if (num(1).eq.0) then
         natoms = natoms + 1
         n = natoms
         if (type.ne.99) then
            itype(n) = type
          else
            itype(n) = 1
          endif
         do 1605 j = 1,3
            rv(j,n) = pos(j)
            if (vel(j).eq.9999.) then
               rv(3+j,n) = 0.
             else
               rv(3+j,n) = vel(j)
             endif
 1605       continue
       else
c
c  modify existing atoms
c
         do 1600 i = num(1),num(2)
         if (oldtype.ne.99.and.itype(i).ne.oldtype) goto 1600
         if (type.ne.99) itype(i) = type
         do 1610 j = 1,3
         if (pos(j).ne.9999.) then
            rv(j,i) = pos(j)
          else
            rv(j,i) = rv(j,i) + delpos(j)
          endif
         if (vel(j).ne.9999.) then
            rv(3+j,i) = vel(j)
          else
            rv(3+j,i) = 0.
          endif
1610     continue
1600     continue
       endif
      goto 100
      end
c***************************************************************
c       deletes vacancies
c
      subroutine delvac
      include 'param.inc'
      include 'utils.inc'
c
c
c       delete vacancies
      j = 0
      do 30 i=1,natoms
      if(itype(i).eq.0)go to 30
      j = j + 1
      itype(j) = itype(i)
      do 20 k=1,6
 20   rv(k,j) = rv(k,i)
 30   continue
      ndel = natoms - j
      natoms = j
      write(6,9000)ndel,natoms
9000  format(//' deleted ',i10,' vacancies'/
     1,i10,' atoms left')
      return
      end
c***************************************************************
c       sorts atoms
c
      subroutine sorter
      include 'param.inc'
      include 'utils.inc'
      dimension idum1(natmax),idum2(natmax),dum1(natmax),
     1                idum3(natmax)
      dimension rmult(3)
c
c
c  set newlst = 1 to insure that if the second neighbor method
c  is in effect that a new neighbor list will be determined
c
      newlst = 1
c
c       order according to non-periodic directions first
c
      ifail=0
      permax = max(perlen(1),perlen(2),perlen(3))
      if(permax.ge.1000.)then
         call m01aaf(perlen,1,3,idum1,idum2,ifail)
         rmult(1) = 10**(4*idum1(1)-3)
         rmult(2) = 10**(4*idum1(2)-3)
         rmult(3) = 10**(4*idum1(3)-3)
         do 40 i=1,natoms
         dum1(i) = rmult(1)*(rv(1,i)-xbound(1)) +
     1      rmult(2)*(rv(2,i)-ybound(1)) + rmult(3)*(rv(3,i)-zbound(1))
 40      continue
      write(6,9010)idum1(1),idum1(2),idum1(3)
9010  format(//' ****** sorting according to following hierarchy:',
     1/'       x      y      z',
     2/1x,3i6)
      else
         do 45 i=1,natoms
         dum1(i) = rv(1,i)**2 + rv(2,i)**2 + rv(3,i)**2
 45      continue
      write(6,9015)
9015  format(//' ****** sorting according to distance from origin')
      endif
      call m01aaf(dum1,1,natoms,idum1,idum2,ifail)
      if(ifail.ne.0)write(6,9020)ifail
 9020 format('    trouble in sorting routine   ',i5)
      do 50 i = 1,natoms
50    idum3(idum1(i)) =  itype(i)
      do 52 i = 1,natoms
52    itype(i) = idum3(i)
      do 60 k = 1,6
      do 62 i = 1,natoms
62    dum1(idum1(i)) = rv(k,i)
      do 64 i = 1,natoms
64    rv(k,i) = dum1(i)
60    continue
      return
      end
c************************************************************************
c
c  chktyp checks that all types used have been defined
c
      subroutine chktyp
      include 'param.inc'
      include 'utils.inc'
      do 90 i = 1,nelmax
90    netype(i) = 0
      do 100 i = 1,natoms
      ityp = itype(i)
      if (ityp.le.0.or.ityp.gt.ntypes) then
         write(6,9001) ityp
9001     format(1x,'** undefined type found. itype = ',i5,' **')
         stop
       end if
      netype(ityp) = netype(ityp) + 1
100   continue
      write(6,9101) (netype(i),i=1,ntypes)
9101  format(1x,'number of atoms of each type:',10i10)
      return
      end
c***************************************************************
      subroutine chkper
      include 'param.inc'
      include 'utils.inc'
      data onlyoi/.true./,twoi/.false./,threei/.false./,fouri/.false./
      data ncalls/0/
      data nx0/0/,ny0/0/,nz0/0/
c
c       check periodicity:
c               allow only one periodic image in the x and y directions
c               allow up to four in the z direction
c
c       x and y:
c       if perlen is lt 2*sqrt(rcutsq), then have troubles with periodicity
c
      rcut = sqrt(rcutsq)
      permin = 2.*rcut
        if(cbeta.ne.0.0.and. perlen(3).lt.permin) then
        write(6,1000)
1000    format (' perlen(3) must be at least twice rcut to run shear')
        stop
        endif
      istop = 0
      do 100 i=1,2
      if(perlen(i).ge.permin)go to 100
      write(6,9230)i
 9230 format('   periodicity is too short in the ',i2,'  direction')
      write(6,9240)permin,perlen(i)
 9240 format('   permin = ',e15.5,'  periodicity = ',e15.5)
      istop = 1
 100   continue
      if(istop.eq.1)stop
c
c       z direction:
c       four levels
c
      if(perlen(3).lt.2.*rcut)then
          onlyoi=.false.
          twoi=.true.
      endif
      if(onlyoi)go to 500
      if(perlen(3).lt.rcut)threei=.true.
      if(perlen(3).lt.2.*rcut/3.)fouri=.true.
      if(perlen(3).lt.0.5*rcut)then
         i=3
         permin = 0.5*rcut
         write(6,9230)i
         write(6,9240)permin,perlen(i)
         stop
      endif
 500  continue
      nx = 1
      ny = 1
      if(fouri)nz = 4
      if(threei.and..not.fouri)nz = 3
      if(twoi.and..not.threei)nz = 2
      if(onlyoi)nz = 1
c
c  print out the number of images used on the first call only
c
      if (ncalls.eq.0)then
      write(6,9035)nx,ny,nz
 9035 format(' # of periodic images   ',i10,5x,i10,5x,i10)
      nx0 = nx
      ny0 = ny
      nz0 = nz
      endif
      ncalls = 1
c
c  print out notice if # of periodic images changes
      if(nx.ne.nx0.or.ny.ne.ny0.or.nz.ne.nz0)then
        nx0 = nx
        ny0 = ny
        nz0 = nz
        write(6,9040)
 9040   format(//' ****** changing # of periodic images ')
        write(6,9035)nx,ny,nz
      endif
      return
      end
c***************************************************************
c
c  this routine produces random velocities with a boltzman
c  distribution appropriate to temp
c    (the center of mass velocity is set to zero)
c
      subroutine velgen(contin)
      include 'param.inc'
      include 'utils.inc'
      logical contin
      integer time
      dimension vcm(3)
      data boltz/8.617e-5/
      data iseed0/87654321/,ngen/1/
      data iseed/12345678/
      namelist /velcard/ temp,iseed,icm,ngen
c
c
c       temperature in degrees kelvin
c       generate velocities in angstroms/picosecond
c
      	if (.not.contin) then
      	temp = 0.0
      	iseed = 0
	icm=.false.
	read(5,velcard)
	temp0=temp
	endif
      if (iseed.lt.0) then
         iseed = time()
       else if (iseed.eq.0) then
         iseed = iseed0
       endif
      call lranst(iseed)
      write(6,9015)temp0,iseed
 9015 format(//' ******  generating velocities according to',
     1' temperature =',g15.5,' with seed',i20)
	iseed=-1
c
c       convert temp to eV
c
      tempin = temp0 * boltz
c
c       mass is in eV-psec**2/Angstroms**2
c
      tmass = 0.0
        do 109 jc=1,3
109     vcm(jc)= 0.0
      do 110 i=1,natoms
      it = itype(i)
      if (amass(it).le.0.0) then
      write(*,*) it,i,amass(it)
      stop
      endif
      vnorm = sqrt(tempin/amass(it))
      do 105 jc=4,6
      rv(jc,i) = vnorm*rgauss()
      vcm(jc-3) = vcm(jc-3) + amass(it)*rv(jc,i)
105   continue
      tmass = tmass + amass(it)
110   continue
c
c       also ensure that the center of mass is not drifting
c
      do 200 jc=1,3
 200  vcm(jc)=vcm(jc)/tmass
      do 300 i=1,natoms
      do 301 jc=4,6
	rv(jc,i)=rv(jc,i)-vcm(jc-3)
 301  continue
 300  continue
c
c  compute the kinetic energy
c
      ekin=0.0
      do 320 i=1,natoms
      v2 = rv(4,i)**2 + rv(5,i)**2  + rv(6,i)**2
 320  ekin = ekin + amass(itype(i))*v2
      ekin = ekin*0.5
c
c  compute the actual temperature
c
      tempact = 7736.6*ekin/float(natoms)
c
c  scale the velocities to get the exact temperature
c
      if (tempact.gt.0.) then
         scalev = sqrt(temp/tempact)
       else
         if (temp.le.0.0) then
            scalev = 0.
           else
            write(6,*)'error in velgen'
            stop
           endif
       endif
      do 350 jc=4,6
      do 350 i=1,natoms
350   rv(jc,i) = rv(jc,i)*scalev
c
c  finally set the initial boundary velocities to zero
c
c     do 400 i = 1,3
c400   bndvel(i) = 0.0
c      return
      end
c***************************************************************
c
c  dummy program to read the namelist card velcard
c
      subroutine veldum(contin)
c       implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'utils.inc'
      logical contin
      namelist /velcard/ temp,iseed,icm,ngen
      if (contin) return
	icm=.false.
      read(5,velcard)
      return
      end
c***************************************************************
c
c  this routine produces random displacements to move the
c  atom positions off of symmetric saddple points
c
      subroutine disgen(dismag)
      include 'param.inc'
      include 'utils.inc'
c
      call lranst(iseed)
c
c  generate small random displacements to
c  push system off of symmetric extremum
c
      write(6,9020)
 9020 format(//' ******  generating small random displacements' )
      do 50 i=1,natoms
      do 50 ic=1,3
50    rv(ic,i) = rv(ic,i) + 2.*dismag*(ranl()-0.5)
      return
      end
c***************************************************************
c
c  this routine checks the namelist intcard to determine the
c  integration scheme to be used and the corresponding parameters.
c  it then calls the appropriate integrator.
c
      subroutine newton(contin)
      include 'param.inc'
      include 'utils.inc'
      logical contin
c
c  define namelist defaults
c
      data tottim,outtim/0.1,0.05/
      data nfmax/0/,stmax/0.0/
      namelist /intcard/ tottim,outtim,iaccur,dt,tol,inte,nfmax,stmax,
     1                   accbnd
c
c  read the data from the namelist
c
      if (.not.contin) read(5,intcard)
        write(6,intcard)
c
c  inte=-1  minimization
c  inte=0   dynamics with slatec routine
c  inte=1   dynamics with nordsieck routine
c  inte=2   just calculate the energy once
c
      if(inte.eq.2)return
      write(6,9010)
9010  format(/,/,' ******  solving problem')
      nsteps = tottim/dt/ngen
      if(inte.gt.0)ipitera = iabs(ipitera)
      if(inte.eq.-1)then
         nsteps = nfmax
         dxmax = stmax
      endif
      if(outtim.gt.tottim)then
         write(6,9020)outtim,tottim
 9020    format('  output time =',g15.5,
     1    ' is greater than tottim=',g15.5,' will set outtim=tottim')
         outtim = tottim
      endif
      noutp = nint(outtim/dt)
c
c  determine the number of outputs before averaging
c
      nequil = eqtim/outtim
c
c  determine the number of degrees of freedom
c
      if (ibdtyp.eq.1) then
          ndegfr = natoms
       else if (ibdtyp.eq.2) then
          ndegfr = natoms + 1
          y(4,ndegfr) = 0.0
          y(5,ndegfr) = 0.0
          y(6,ndegfr) = 0.0
       else
          write(6,9030)ibdtyp
9030      format(1x,'undefined boundary type: ibdtyp =',i3)
          stop
       end if
c
c  check that the number of degrees of freedom does not exceed natmax
c
      if (ndegfr.gt.natmax) then
         write(6,9901) ndegfr,natmax
9901     format(1x,'ndegfr=',i6,' is greater than natmax=',i6)
         stop
       end if
c
c  compute the scaled coordinates to pass to solvers
c
      call scale_sub
c
c  call the requested solver
c
      if (inte.lt.0)then
         write(6,9100)
9100     format(//' minimizing energy')
      else
         write(6,9110)
9110     format(//' solving newton''s equation')
      endif
      if (inte.eq.-1)call minimize
      if (inte.eq. 0)call slatint
      if (inte.eq. 1)call nordint
c
      write(6,9120)
9120  format(' ****** solution finished')
c
c  unscale the coordinates
c
      call uscale
      return
      end
c**************************************************************************
c
c  this routine minimizes the potential energy using
c  va08a (conjugate gradients)
c
      subroutine minimize
      include 'param.inc'
      include 'utils.inc'
      external vafunc
      data df/-.1/
c       df = order of magnitude of change of f
c       eps = convergence criterion (smallest change in variables)
c       dxmax = largest stepsize allowed
c
      eps = -tol
c
c  zero all velocities
c
      do 5 i=1,ndegfr
      do 5 j=4,6
5     y(j,i) = 0.0
c
c      if dxmax has not been set by user, then set dxmax to be a fraction of
c      the nearest neighbour distance
c
      if(dxmax.eq.0.0)then
        distnn = alat/2.
        dxmax = 0.1*distnn
      endif
c
      n = 3*ndegfr
c
c  call minimizer once if nsteps was set by the user or until time
c  is up if it has not been set to a non-zero value
c
      if (nsteps.gt.0) then
         write(6,9110)nsteps,eps,dxmax
         call va08a(vafunc,n,ener,df,eps,nsteps,ipitera,
     1              dxmax,iterm)
c
c  set mnterm to 1 if the minimization was terminated by nfmax
c
         mnterm = iabs(iterm)
       else
c
c  determine the maximum number of steps that the time limit allows
c
         iterm = 0
 25      call trmain(timlft)
         nsteps = (timlft-2.0*frctmx)/(1.10*frctmx)
c
c  check to see if there is time for (more) iterations
c
         if(nsteps.lt.100.and.iterm.ne.0)goto 26
         write(6,9110)nsteps,eps,dxmax
9110     format('  calling minimizer (va08a)  nfmax = ',i10,
     1          '  eps = ',g15.5,'  dxmax = ',g15.5)
         call va08a(vafunc,n,ener,df,eps,nsteps,ipitera,
     1              dxmax,iterm)
c
c  determine if va08a terminated because of time limit
c
         if(iterm.eq.-1)go to 25
26       continue
         mnterm = iterm
       end if
      return
      end
c**************************************************************************
c
c       this is a slightly modified version of va08a
c
      subroutine va08a (funct,n,fl,dfn,epps,maxfn,iprint,dsmax,
     1    iterm)
      include 'param.inc'
      include 'utils.inc'
      dimension x(6*natmax)
      equivalence (x(1),y(1,1))
      external funct
c--     va08a  has been modified from the original harwell routine.
c**      the original array x,g,and s are now passed to and from this
c**      routine via common/pointer variables
c--          eps has been replaced in argument list by epps and
c--     dsmax has been added to argument list.
c--     the former  eps=abs(epps).  dsmax will not be referenced
c--     unless  epps.lt.0.  in that case dsmax will be the magnitude
c--     (in euclidean norm sense) of the largest step size that will
c--     be taken during iteration in moving from any currently best
c--     point to another point.
c--     iterm = 0 if normal termination
c--     iterm = -1 if terminated due to time limit
c--          this modification is brought to you by
c--          thomas jefferson.
      common /mincom/ g(3*natmax),s(3*natmax)
c

      ntry = 0
      nmax = 0
      iterm=0
        if( idynper(1).lt.0.and.idynper(2).lt.0
     1     .and.idynper(3).lt.0) then
        iper=4
        elseif(idynper(1).lt.0.and.idynper(2).lt.0) then
        iper=3
        elseif(idynper(1).lt.0.and.idynper(3).lt.0) then
        iper=2
        elseif(idynper(3).lt.0.and.idynper(2).lt.0) then
        iper=1
        else
        iper=0
        endif
      if (iprint.ne.0) write (6,250) n,fl,dfn,epps,maxfn
      eps=abs(epps)
      z=fl
      itn=0
      call funct (n,fl)
      flast=fl
      smax=0.
      fy=1.0e+6
      imax=0
      igrad=1
c--         igrad=1  means that array g(*) contains gradient of f at x()
      ifn=1
      df=dfn
      if (dfn.eq.0.) df=fl-z
      if (dfn.lt.0.) df=abs(df*fl)
      if (df.le.0.) df=1.
110   continue
      do 120 i=1,n
120   s(i)=0.
      smaxmx=0.
      gg=1.
      do 220 icyc=1,n
      itn=itn+1
c               gg = 0. + (g,g)
      gglast=gg
130   gg=sdot(n,g,1,g,1)
      z=gg/gglast
      if (iprint.eq.0) go to 140
      if (mod(itn,iabs(iprint)).ne.0) go to 140
      write (6,260) icyc,itn,ifn,fl,gg,imax,smax
      if (iprint.lt.0) go to 140
      write (6,270) (x(i),i=1,n)
      write (6,270) (g(i),i=1,n)
140   continue
         if (z.ne.0.) go to 150
c--         if g is gradient at current  x we are through.  otherwise,
c--evaluate gradient at current x and make test again.
         if (igrad.eq.1) go to 240
      call funct (n,fl)
      igrad=1
      ifn=ifn+1
      go to 130
c--
150   imax=1
      do 160 i=1,n
      s(i)=z*s(i)-g(i)
      if (abs(s(i)) .gt. abs(s(imax))) imax=i
160   continue
c--           now  imax  is index of the component of s with largest mag
      gs=sdot(n,g,1,s,1)
c               gs = 0. + (g,s)
      if (gs.ge.0.) go to 110
c--     if stepsize is to restricted in length to dsmax, then find the
c--     maximum alpha corresponding to the direction given by s vector.
      if (epps.lt.0.) then
c--                                 ss = 0. + (s,s)
          ss=sdot(n,s,1,s,1)
          alphmx=dsmax/sqrt(ss)
      endif
      gs0=gs
      if ((flast-fl).gt.0.) df=flast-fl
      flast=fl
      alpha=-2.*df/gs
      alpha0=0.
c--     alpha0  is the scalar factor multiplying the  vector s(*)
c--               thereby giving the step increment from x(i) at this pt
170   continue
         if (ifn.ge.maxfn) go to 240
c--     if maximum stepsize is to be restricted, then restrict it by
c--               limiting alpha.
      ntry = ntry + 1
      if (alpha.ge.alphmx)nmax = nmax + 1
      if (epps.lt.0.) alpha=min(alphmx,alpha)
      alpha0=alpha0+alpha
c       the following is modified to handle the rv array which has both
c       positions and velocities in it ---------msd 11/16/83
c       do x,y, and z coordinates separately
c       s is incremented by 3, but x is incremented by 6
      call saxpy (ndegfr,alpha,s(1),3,x(1),6)
      call saxpy (ndegfr,alpha,s(2),3,x(2),6)
      call saxpy (ndegfr,alpha,s(3),3,x(3),6)
      call funct (n,fy)
      ifn=ifn+1
      igrad=1
      icon=0
      if (abs(alpha*s(imax)).gt.eps) icon=1
      gys=sdot(n,g,1,s,1)
         if (fy.ge.fl) go to 180
c--     gys/gs0  is equal to the ratio of the projection of the current
c--              gradient in the s-vector direction  to the
c--              projection in the s direction of the gradient at the
c--              start of the minimization in the s direction.
         if (abs(gys/gs0).le..1) go to 200
         if (gys.gt.0.) go to 180
      z=4.
      if (gs.lt.gys) z=gys/(gs-gys)
      if (z.gt.4.) z=4.
      alpha=alpha*z
      fl=fy
      gs=gys
         go to 170
c--               come here if we have gone too far in the s-vector
c--          direction.  i.e. either function value is larger or
c--          gradient dot s is now positive.
180   iset=1
         if (icon.eq.0.or.(alpha.eq.0.)) go to 190
      z=3.*(fl-fy)/alpha+gys+gs
         if (z.eq.0.) go to 190
      w=sqrt(1.-gs/z*gys/z)*abs(z)
      z=1.-(gys+w-z)/(2.*w+gys-gs)
      iset=2
c--      if iset=2 we will continue minimizing in  s direction.
c-- otherwise(iset=1), set x-vector back to its previous value only if
c-- the current function value fy  is no improvement.
190      if (iset.eq.1.and.fy.lt.fl) go to 200
c
c       the following is modified to handle the rv array which has both
c       positions and velocities in it ---------msd 11/16/83
c       do x,y, and z coordinates separately
c       s is incremented by 3, but x is incremented by 6
      call saxpy (ndegfr,-alpha,s(1),3,x(1),6)
      call saxpy (ndegfr,-alpha,s(2),3,x(2),6)
      call saxpy (ndegfr,-alpha,s(3),3,x(3),6)
      alpha0=alpha0-alpha
      igrad=2
      alpha=alpha*z
         if (iset.eq.2) go to 170
         go to 210
c--               come here if new f is .lt. previous best f and
c--     projection of gradient in s-vector direction has been
c--     sufficiently reduced.
c--            also come here if we have gone too far in s-direction,
c--      yet the change in the unknown parameters is too small to be
c--      of significance.
200   fl=fy
c--       if parameters have changed sufficiently during this
c--  minimization in this s-direction, then go on to next direction.
c--  otherwise, either start again in gradient direction (icyc=1) or
c--  we are finished
210   smax=alpha0*s(imax)
      smaxmx=max(smaxmx,abs(smax))
         if (abs(smax).le..05*smaxmx) go to 230
         if (abs(smax).le.eps) go to 230
220   continue
         go to 110
c--   start again in gradient direction unless latest minimization was
c--         in gradient direction.
230      if (icyc.ne.1) go to 110
c--                   the end
240   fl=min(fl,fy)
      if (ifn.ge.maxfn)then
         write(6,9240)
9240     format(' ********warning: minimization terminated due to',
     1     'time limit *********')
         iterm=-1
      endif
      if (igrad.eq.2) call funct (n,fl)
      if (igrad.eq.2) ifn=ifn+1
c
      if (iprint.eq.0) return
      gg=sdot(n,g,1,g,1)
      icyc=n+1
      pmax = float(nmax)/float(ntry)
      write (6,260) icyc,itn,ifn,fl,gg,imax,smax
      write (6,261) pmax
      if (iprint.le.0) return
      write (6,270) fl
      write (6,270) (x(i),i=1,n)
      write (6,270) (g(i),i=1,n)
      return
250   format (20h entry to va08a,  n=,i5,3h f=,e12.4,5h dfn=,e12.4,6h ep
     1ps=,e12.4,7h maxfn=,i5)
260   format (17h  in va08a--icyc=,i5,5h itn=,i5,5h ifn=,i5,3h f=,e19.11
     1,10h g dot g =,e14.6,6h imax=,i5,6h smax=,e12.4)
261   format(' in va08a--fraction of steps that were limited by',
     1  ' maximum step size=',/,1x,e12.4)
270   format ((8e15.7))
      end
c**************************************************************************
c   constrain.f
c$Log: utils.f,v $
cRevision 1.33  2004/10/08 18:13:12  baskes
cin verlet replaced dyn.inc with utils.inc
c
cRevision 1.32  2003/12/05 22:03:12  baskes
cupdating
c
cRevision 1.31  2003/11/10 18:42:49  baskes
cchanged rcutsq to rctsqn
cmakes neighbor region larger
c
cRevision 1.30  2003/08/20 17:17:57  baskes
cnew .inc
c
cRevision 1.29  2003/07/01 03:37:37  baskes
cgetting ready to change utils.inc
c
cRevision 1.28  2002/11/04 18:29:26  baskes
cchanged print
c
cRevision 1.27  2002/08/25 20:41:41  baskes
cfixed print formats
cremoved getcell print
c
cRevision 1.26  2002/05/08 16:25:53  baskes
cfixed iper
c
cRevision 1.25  2002/05/07 21:38:13  baskes
cas of 5/2/2002
c
cRevision 1.24  2002/04/30 20:49:49  baskes
caveraging fixed again
c
cRevision 1.23  2002/04/30 18:09:17  baskes
cfix averaging start
c
cRevision 1.22  2002/04/29 18:43:26  baskes
cadded RMS for force
c
cRevision 1.21  2002/04/07 00:47:24  baskes
cmove restart to utils
cadd uniform rate of periodic vector change
c
cRevision 1.20  2001/12/04 22:39:34  baskes
cCm fixed
c
cRevision 1.19  2001/12/04 18:53:22  baskes
cadd icm zero CM velocity
c
cRevision 1.18  2001/06/12 15:38:31  baskes
cadded acf
c
cRevision 1.17  2001/01/10 21:18:28  baskes
cfixed short period for nmeth=-2
c
cRevision 1.16  2001/01/08 23:19:10  baskes
cfixed more dimensions for timing
c
cRevision 1.15  2001/01/08 23:12:25  baskes
cfixed dimension for timing
c
cRevision 1.14  2001/01/08 23:11:55  baskes
cfixed dimension for timing
c
cRevision 1.13  2001/01/08 22:44:38  baskes
cput in header
c
cRevision 1.12  2001/01/08 22:43:21  baskes
cmoved data statements out
c
cRevision 1.11  2000/11/21 17:20:13  baskes
cdouble rand
c
cRevision 1.10  2000/11/20 20:47:15  baskes
c*** empty log message ***
c
cRevision 1.9  2000/11/20 17:31:03  baskes
cmoved nmeth3
c
cRevision 1.8  2000/11/20 17:20:58  baskes
cgetnei working for mono
c
cRevision 1.7  2000/11/17 20:13:29  baskes
cmoved gneigh
c
cRevision 1.6  2000/11/16 23:21:34  baskes
caccel moved to utils
c
cRevision 1.5  2000/11/16 21:37:36  baskes
cmoved initlat to nords to utils
c
cRevision 1.4  2000/11/16 16:43:23  baskes
croutines moved from dyn
c
cRevision 1.3  2000/11/15 21:40:55  baskes
cnow contains lutils and constrain
c
cRevision 1.7  2000/11/09 21:39:23  baskes
cmoved utils.inc to top
c
cRevision 1.6  2000/09/25 21:35:44  baskes
cimplement circumferential velocity by setting forces
c
cRevision 1.5  2000/09/21 18:35:23  baskes
csign fixed
c
cRevision 1.4  2000/09/21 16:10:10  baskes
c*** empty log message ***
c
cRevision 1.3  2000/09/21 16:09:20  baskes
c*** empty log message ***
c
cRevision 1.2  2000/09/21 16:07:19  baskes
cadd mode 7 for circumferential strain rate
c
c*****************************************************************************
c
c  subroutine setfix determines the constraints on each atom
c
      subroutine setfix(contin,lastconf,nfixcard)
      include 'param.inc'
      include 'utils.inc'
      logical contin,lastconf
      integer type
      dimension num(2),vector(3),dvecdt(3)
      namelist /fixcard/ num,type,mode,vector,dvecdt,
     1                   xmin,xmax,ymin,ymax,zmin,zmax
c
c  set the last time fixfor was called (timefx) to t0
c
      timefx = t0
      nsize1=natoms
c
c  do not clear on continuation
c
      if (contin) goto 1000
c
C  clear the flag for the existence of constraints
c
      nfixes = 0
      ifxal1 = 0
      ifxal2 = 0
c
c  clear the fixcard data arrays
c
      do 100 i = 0,100
      modefx(i) = 0
      vectfx(1,i) = 0.
      vectfx(2,i) = 0.
      vectfx(3,i) = 0.
      dvctfx(1,i) = 0.
      dvctfx(2,i) = 0.
      dvctfx(3,i) = 0.
100   continue
1000  continue
c
c     if using the last configuration for start, preserve the work by
c     external forces
c
      if (.not.lastconf.and.ifxal2.ne.0) then
         do 200 i = 1,natoms
200      workfx(i) = 0.0
       endif
c
c  clear the namelist variables to read the next card
c
      num(1) = 0
      num(2) = 0
      type = 0
      mode = -9999
      xmin = -9999.
      xmax = 9999.
      ymin = -9999.
      ymax = 9999.
      zmin = -9999.
      zmax = 9999.
      vector(1) = 9999.
      vector(2) = 9999.
      vector(3) = 9999.
      dvecdt(1) = 9999.
      dvecdt(2) = 9999.
      dvecdt(3) = 9999.
c
c  read the fix card
c
      read(5,fixcard)
c
c  see if anything is to be fixed
c
      if (mode.eq.-9999) goto 2000
      nfixes = nfixes + 1
c
c  if not done already
c
      if (ifxal1.eq.0) then
         ifxal1 = 1
         do 300 i = 1,nsize1
300      ipntfx(i) = 0
       endif
      if (ifxal2.eq.0.and.mode.eq.4) then
         ifxal2 = 1
         do 400 i = 1,nsize1
         workfx(i) = 0.
         posfx(1,i) = 0.
         posfx(2,i) = 0.
         posfx(3,i) = 0.
400      continue
       endif
c
c  write out the fixcard
c
      write(6,fixcard)
c
c  if vector refers to a direction, normalize it.
c
      if (mode.eq.1.or.mode.eq.2) then
         vnorm = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)
         vector(1) = vector(1)/vnorm
         vector(2) = vector(2)/vnorm
         vector(3) = vector(3)/vnorm
       endif
c
c  store the fixcard data
c
      modefx(nfixes) = mode
      if(vector(1).ne.9999.) vectfx(1,nfixes) = vector(1)
      if(vector(2).ne.9999.) vectfx(2,nfixes) = vector(2)
      if(vector(3).ne.9999.) vectfx(3,nfixes) = vector(3)
      if(dvecdt(1).ne.9999.) dvctfx(1,nfixes) = dvecdt(1)
      if(dvecdt(2).ne.9999.) dvctfx(2,nfixes) = dvecdt(2)
      if(dvecdt(3).ne.9999.) dvctfx(3,nfixes) = dvecdt(3)
c
c  determine to which atoms these constraints apply and set the ipntfx
c  values for those atoms to this set of constraints
c
c  branch depending on whether the atoms are specified by number or position
c
      if (num(1).ne.0) goto 1500
c
c  specified by position
c
      do 1100 i = 1,natoms
      if (rv(1,i).lt.xmin.or.rv(1,i).gt.xmax) goto 1100
      if (rv(2,i).lt.ymin.or.rv(2,i).gt.ymax) goto 1100
      if (rv(3,i).lt.zmin.or.rv(3,i).gt.zmax) goto 1100
      if (type.ne.0.and.itype(i).ne.type) goto 1100
      ipntfx(i) = nfixes
      if (mode.eq.4) then
         posfx(1,i) = rv(1,i)
         posfx(2,i) = rv(2,i)
         posfx(3,i) = rv(3,i)
       endif
1100  continue
      goto 1000
c
c  specified by number
c
1500  continue
      if (num(2).eq.0) num(2) = num(1)
      do 1600 i = num(1),num(2)
      ipntfx(i) = nfixes
      if (mode.eq.4) then
         posfx(1,i) = rv(1,i)
         posfx(2,i) = rv(2,i)
         posfx(3,i) = rv(3,i)
       endif
1600  continue
      goto 1000
2000  continue
      if (nfixes.eq.0) return
c
c  set the velocities of any fixed atoms to 0 in the appropriate
c  direction
c
      do 2100 i = 1,natoms
      if (modefx(ipntfx(i)).eq.1) then
      vdvec = rv(4,i)*vectfx(1,ipntfx(i)) + rv(5,i)*vectfx(2,ipntfx(i))
     1         + rv(6,i)*vectfx(3,ipntfx(i))
         rv(4,i) = rv(4,i) - vdvec*vectfx(1,ipntfx(i))
         rv(5,i) = rv(5,i) - vdvec*vectfx(2,ipntfx(i))
         rv(6,i) = rv(6,i) - vdvec*vectfx(3,ipntfx(i))
      elseif (modefx(ipntfx(i)).eq.2) then
      vdvec = rv(4,i)*vectfx(1,ipntfx(i)) + rv(5,i)*vectfx(2,ipntfx(i))
     1         + rv(6,i)*vectfx(3,ipntfx(i))
         rv(4,i) = vdvec*vectfx(1,ipntfx(i))
         rv(5,i) = vdvec*vectfx(2,ipntfx(i))
         rv(6,i) = vdvec*vectfx(3,ipntfx(i))
      elseif (modefx(ipntfx(i)).eq.3) then
         rv(4,i) = 0.0
         rv(5,i) = 0.0
         rv(6,i) = 0.0
       elseif (modefx(ipntfx(i)).eq.5) then
         if(dvctfx(1,ipntfx(i)).eq.0.)  rv(4,i) = vectfx(1,ipntfx(i))
         if(dvctfx(2,ipntfx(i)).eq.0.)  rv(5,i) = vectfx(2,ipntfx(i))
         if(dvctfx(3,ipntfx(i)).eq.0.)  rv(6,i) = vectfx(3,ipntfx(i))
       elseif (modefx(ipntfx(i)).eq.6) then
         if(dvctfx(1,ipntfx(i)).eq.0.)
     1	  rv(4,i) = rv(1,i)*vectfx(1,ipntfx(i))
         if(dvctfx(2,ipntfx(i)).eq.0.)
     1	  rv(5,i) = rv(2,i)*vectfx(2,ipntfx(i))
         if(dvctfx(3,ipntfx(i)).eq.0.)
     1	  rv(6,i) = rv(3,i)*vectfx(3,ipntfx(i))
       elseif (modefx(ipntfx(i)).eq.7) then
         if(dvctfx(1,ipntfx(i)).eq.0.)
     1   rv(4,i) = rv(2,i)*vectfx(1,ipntfx(i))
         if(dvctfx(2,ipntfx(i)).eq.0.)
     1   rv(5,i) = rv(1,i)*vectfx(2,ipntfx(i))
	if(dvctfx(3,ipntfx(i)).eq.0.)  rv(6,i)=vectfx(3,ipntfx(i))
       elseif (modefx(ipntfx(i)).eq.8) then
         if(dvctfx(1,ipntfx(i)).eq.0.)
     1   rv(4,i) = rv(3,i)*vectfx(1,ipntfx(i))
         if(dvctfx(3,ipntfx(i)).eq.0.)
     1   rv(6,i) = rv(1,i)*vectfx(3,ipntfx(i))
        if(dvctfx(2,ipntfx(i)).eq.0.)  rv(5,i)=vectfx(2,ipntfx(i))
        endif
2100  continue
c
c  last card read  -  print out current constraints
c
      if(iconst.eq.1)then
        write(6,9010)
9010    format(1x,/,1x,'******  the following constraints apply')
        write(6,9020)
9020    format(1x,'   num  type  mode',15x,'vector',35x,'dvecdt')
        do 2200 i = 1,natoms
        if (modefx(ipntfx(i)).eq.0) goto 2200
        write(6,9030) i,itype(i),modefx(ipntfx(i)),
     1              (vectfx(j,ipntfx(i)),j=1,3),
     2              (dvctfx(j,ipntfx(i)),j=1,3)
9030    format(1x,3i6,3g12.4,5x,3g12.4)
2200    continue
        write(6,9035) nfixes,nsize1
9035    format(1x,' number of fixes ',2i6)
        write(6,9040)
9040    format(1x,'end of constraints',/,1x)
      endif
      nfixcard=1
      return
      end
c***************************************************************************
c
c  fixfor implements the force fixing defined by setfix
c
      subroutine fixfor(time)
      include 'param.inc'
      include 'utils.inc'
      dimension delta(3)
c
c  if no constraints, return
c
      if (nfixes.eq.0) return
c
c  fix vector for helping to calculate external work
c
c     per12(1) = 0.5*perlen(1)
c     per12(2) = 0.5*perlen(2)
c     per12(3) = 0.5*perlen(3)
c
c  increment the vectors with time
c
      do 100 i = 1,nfixes
      vectfx(1,i) = vectfx(1,i) + (time-timefx)*dvctfx(1,i)
      vectfx(2,i) = vectfx(2,i) + (time-timefx)*dvctfx(2,i)
      vectfx(3,i) = vectfx(3,i) + (time-timefx)*dvctfx(3,i)
100   continue
c
c  loop through all the atoms and implement the constraints
c
      do 1000 i = 1,natoms
c
c  skip if no constraints on this atom
c
      if (ipntfx(i).eq.0) goto 1000
c
c  branch to the appropriate contraint mode
c
      if (modefx(ipntfx(i)).eq.0) then
	goto 1000
c
c  constrain to a plane
c
      elseif (modefx(ipntfx(i)).eq.1) then
      dotprd = f(1,i)*vectfx(1,ipntfx(i)) + f(2,i)*vectfx(2,ipntfx(i))
     1       + f(3,i)*vectfx(3,ipntfx(i))
      f(1,i) = f(1,i) - dotprd*vectfx(1,ipntfx(i))
      f(2,i) = f(2,i) - dotprd*vectfx(2,ipntfx(i))
      f(3,i) = f(3,i) - dotprd*vectfx(3,ipntfx(i))
c
c  constrain to a line
c
      elseif (modefx(ipntfx(i)).eq.2) then
 	if(ibdtyp.eq.1) then
      dotprd = f(1,i)*vectfx(1,ipntfx(i)) + f(2,i)*vectfx(2,ipntfx(i))
     1       + f(3,i)*vectfx(3,ipntfx(i))
      f(1,i) = dotprd*vectfx(1,ipntfx(i))
      f(2,i) = dotprd*vectfx(2,ipntfx(i))
      f(3,i) = dotprd*vectfx(3,ipntfx(i))
	else
      dotprd = f(1,i)*vectfx(1,ipntfx(i))*perlen(1)**2
     1	     + f(2,i)*vectfx(2,ipntfx(i))*perlen(2)**2
     1       + f(3,i)*vectfx(3,ipntfx(i))*perlen(3)**2
      f(1,i) = dotprd*vectfx(1,ipntfx(i))/perlen(1)**2
      f(2,i) = dotprd*vectfx(2,ipntfx(i))/perlen(2)**2
      f(3,i) = dotprd*vectfx(3,ipntfx(i))/perlen(3)**2
	endif
c
c  constrain to a point
c
      elseif (modefx(ipntfx(i)).eq.3) then
      f(1,i) = 0.0
      f(2,i) = 0.0
      f(3,i) = 0.0
c
c  add a force
c
      elseif (modefx(ipntfx(i)).eq.4) then
      f(1,i) = f(1,i) + vectfx(1,ipntfx(i))
      f(2,i) = f(2,i) + vectfx(2,ipntfx(i))
      f(3,i) = f(3,i) + vectfx(3,ipntfx(i))
      delta(1) = rv(1,i)-posfx(1,i)
      delta(2) = rv(2,i)-posfx(2,i)
      delta(3) = rv(3,i)-posfx(3,i)
      delta(1) = delta(1) - perlen(1)*nint(delta(1)/perlen(1))
      delta(2) = delta(2) - perlen(2)*nint(delta(2)/perlen(2))
      delta(3) = delta(3) - perlen(3)*nint(delta(3)/perlen(3))
      workfx(i) = workfx(i) +
     1            vectfx(1,ipntfx(i))*delta(1) +
     2            vectfx(2,ipntfx(i))*delta(2) +
     3            vectfx(3,ipntfx(i))*delta(3)
      e(i) = e(i) - workfx(i)
      posfx(1,i) = rv(1,i)
      posfx(2,i) = rv(2,i)
      posfx(3,i) = rv(3,i) 
c
c     fix velocity in a direction
c
      elseif (modefx(ipntfx(i)).eq.5.or.modefx(ipntfx(i)).eq.6) then
      if(dvctfx(1,ipntfx(i)).eq.0.) f(1,i)=0.0
      if(dvctfx(2,ipntfx(i)).eq.0.) f(2,i)=0.0
      if(dvctfx(3,ipntfx(i)).eq.0.) f(3,i)=0.0
c
c     fix velocity to be circumferential by setting forces
c	v_x=y*vect_x
c	v_y=x*vect_y
c	v_z=  vect_z
c
c
      elseif (modefx(ipntfx(i)).eq.7) then
	if(dvctfx(1,ipntfx(i)).eq.0.) then
	f(1,i)=rv(5,i)*vectfx(1,ipntfx(i))*amass(itype(i))
	endif
	if(dvctfx(2,ipntfx(i)).eq.0.) then
	f(2,i)=rv(4,i)*vectfx(2,ipntfx(i))*amass(itype(i))
	endif
	if(dvctfx(3,ipntfx(i)).eq.0.) f(3,i)=0.0
c
c     fix velocity to be circumferential by setting forces
c       v_x=z*vect_x
c       v_z=x*vect_z
c       v_y=  vect_y
c
c
      elseif (modefx(ipntfx(i)).eq.8) then
        if(dvctfx(1,ipntfx(i)).eq.0.) then
        f(1,i)=rv(6,i)*vectfx(1,ipntfx(i))*amass(itype(i))
        endif
        if(dvctfx(3,ipntfx(i)).eq.0.) then
        f(3,i)=rv(4,i)*vectfx(3,ipntfx(i))*amass(itype(i))
        endif
        if(dvctfx(2,ipntfx(i)).eq.0.) f(2,i)=0.0
      else
      write(6,9101) modefx(i)
9101  format(1x,'undefined modefx',i6)
      stop
	endif
1000  continue
c
c  reset the last time called variable
c
      timefx = time
      return
      end
c**************************************************
c
c   settmp determines the temperature control parameters
c
      subroutine settmp(contin)
c	modified 10/28/98 by MIB to include option of T control in xyz
      include 'param.inc'
      include 'utils.inc'
      integer time
      logical contin,first
      dimension qdefault(4),num(2),itemp(3),vcm(3)
      data qdefault/1730.,40.,0.1,10000./,destmp/300./,tmptim/1.0/
	data itemp/3*1/,vcm/3*0./
c     data ialloc/0/
      namelist /tmpcard/ ifxtmp,follow,iseed
      namelist /regcard/ q,tmptim,destmp,xmin,xmax,ymin,ymax,zmin,zmax,
     1   num,itemp,vcm,e0,v0,pzz
	nsize1=natoms
c
c
c  read the temperature control card
c
      read(5,tmpcard)
c
      write(6,9229)
9229   format(//,' ******  energy/temperature conditions ')
c
      if (ifxtmp.eq.0)then
         write(6,9230)
 9230    format(/,'  constant energy',/)
c
c  read in regcards (to be ignored)
c
100      destmp = -9999.
         read(5,regcard)
         if(destmp.eq.-9999.)return
         go to 100
      endif
c
c  set time variable for fixtmp
c
      timetmp = t0
      tmpreg(0) = 0.
	itempg(0,1)=1
	itempg(0,2)=1
	itempg(0,3)=1
c     ialloc = 1
      do 300 i = 1,nsize1
300   ipntrg(i) = 0
c
c  set representative mass for determination of q values from tmptim
c
      repmass = amass(1)
      do 305 ityp=1,ntypes
305   repmass = min(repmass,amass(ityp))
c
c  read in cards to define regions
c
      if(.not.contin)nregs = 0
      first = .true.
c
c  read in parameters
c
1000  q = qdefault(ifxtmp)
      tmptim = 9999.
      destmp = -9999.
      xmin = -9999.
      xmax = 9999.
      ymin = -9999.
      ymax = 9999.
      zmin = -9999.
      zmax = 9999.
      num(1)=0
      num(2)=0
	itemp(1)=1
	itemp(2)=1
	itemp(3)=1
c
      read(5,regcard)
      if(destmp.eq.-9999.)goto 2000
      if(contin.and.first)nregs = 0
      first = .false.
      nregs = nregs + 1
c
      write(6,9233)nregs
9233  format(/' temperature control region ',i5)
      write(6,regcard)
c
c  store region information
c
      if(tmptim.ne.9999.)then
          if(ifxtmp.eq.1)qreg(nregs) = destmp*tmptim/(2.2*repmass)
          if(ifxtmp.eq.2)qreg(nregs) = destmp*tmptim**2/(188.*repmass)
          if(ifxtmp.eq.3)qreg(nregs) = 3.*repmass/tmptim
          if(ifxtmp.eq.4)qreg(nregs) = 600./tmptim**2
      else
          qreg(nregs) = q
      endif
      tmpreg(nregs) = destmp
	do 777 kkk=1,3
	vcmrg(nregs,kkk)=vcm(kkk)
	itempg(nregs,kkk)=itemp(kkk)
777	continue
	if(num(1).eq.0) then
      xcenter(nregs) = (xmin + xmax)/2.
      xwidth(nregs) =  (xmax - xmin)/2.
      ycenter(nregs) = (ymin + ymax)/2.
      ywidth(nregs) =  (ymax - ymin)/2.
      zcenter(nregs) = (zmin + zmax)/2.
      zwidth(nregs) =  (zmax - zmin)/2.
c
c  if following by atom number, set up indirection pointer
c
      if(follow) then
        do 2500 i=1,natoms
          if(abs(rv(1,i)-xcenter(nregs)).gt.xwidth(nregs))goto 2501
          if(abs(rv(2,i)-ycenter(nregs)).gt.ywidth(nregs))goto 2501
          if(abs(rv(3,i)-zcenter(nregs)).gt.zwidth(nregs))goto 2501
          ipntrg(i) = nregs
2501    continue
2500    continue
      endif
	else
      follow=.true.
      if (num(2).eq.0) num(2) = num(1)
      do 1600 i = num(1),num(2)
          ipntrg(i) = nregs
1600  continue
	endif
c
c  set initial value (used only for ifxtmp=2,4)
c
      drag(nregs) = 0.0
c
	go to 1000
c
2000  continue
c
      if(nregs.eq.0)then
        write(6,9555)
9555    format(' no regions defined:',/,
     .        ' must at least specify a temperature ')
        stop
      endif
c
c  if following by atom number, set up indirection pointer
c
c     if(follow) then
c       do 2500 i=1,natoms
c       ipntrg(i) = 0
c       do 2501 ireg=1,nregs
c         if(abs(rv(1,i)-xcenter(ireg)).gt.xwidth(ireg))goto 2501
c         if(abs(rv(2,i)-ycenter(ireg)).gt.ywidth(ireg))goto 2501
c         if(abs(rv(3,i)-zcenter(ireg)).gt.zwidth(ireg))goto 2501
c         ipntrg(i) = ireg
c501    continue
c500    continue
c     endif
c
c  set random number seed for langevin constraint
c
      if(ifxtmp.eq.3) then
        if (iseed.lt.0) then
           iseed = time()
        else if (iseed.eq.0) then
           iseed = iseed0
        endif
        call lranst(iseed)
        write(6,9015)iseed
 9015   format(//' ****** stochastic force seed set at ',i20)
      endif
	if(ifxtmp.eq.4) then
	ireg=1
	write(6,9017)
9017	format('  region ','initial energy    initial pressure    initial',
     1       ' volume')
	write(6,9016) ireg,e0,pzz,v0
9016	format(i10,3f20.5)
	endif
        if(icm) then
        do 299  ireg=0,nregs
      tmass = 0.0
        do 109 jc=1,3
109     vcm(jc)= 0.0
      do 110 i=1,natoms
        if(ipntrg(i).ne.ireg) go to 110
      it = itype(i)
      do 105 jc=4,6
      vcm(jc-3) = vcm(jc-3) + amass(it)*rv(jc,i)
105   continue
      tmass = tmass + amass(it)
110   continue
              do 201 jc=1,3
        vcm(jc)=vcm(jc)/tmass
        vcm(jc)=vcm(jc)-vcmrg(ireg,jc)
201     continue
      do 302 i=1,natoms
        if(ipntrg(i).ne.ireg) go to 302
      do 301 jc=1,3
301     rv(jc+3,i)=rv(jc+3,i)-vcm(jc)
302     continue
299     continue
        endif
      return
      end
 
c************************************************************************
c
c  this routine applies temperature control to the equations of motion
c  by adding a drag term
c  two methods implemented:  both have f = f - drag*v
c      ifxtmp = 1:   drag = (actual temp - desired tmp)/q
c      ifxtmp = 2:   d(drag)/dt = (actual temp - desired tmp)/q  (Hoover)
c      ifxtmp = 3:   langevin (Biswas & Hamann)
c      ifxtmp = 4:   hugoniostat (Maillet PRE 63, 016121 (2000)
c
      subroutine fixtmp(time)
      include 'param.inc'
      include 'utils.inc'
      data boltz/8.617e-5/
c
c  modify the forces to enforce temperature control if appropriate
c
      if (ifxtmp.eq.0) return
c
c  if following particle numbers then update indirection pointer
      if(.not.follow)then
        do 200 i=1,natoms
        ipntrg(i) = 0
        do 201 ireg=1,nregs
          if(abs(rv(1,i)-xcenter(ireg)).gt.xwidth(ireg))goto 201
          if(abs(rv(2,i)-ycenter(ireg)).gt.ywidth(ireg))goto 201
          if(abs(rv(3,i)-zcenter(ireg)).gt.zwidth(ireg))goto 201
          ipntrg(i) = ireg
201    continue
200    continue
      endif
c
c  calculate drag for each region and apply to atoms
c
      deltat = time-timetmp
      drag(0) = 0.
c
c  standard drag
c
      if(ifxtmp.eq.1)then
c
c  calculate the temperature of each region
c
        do 299 ireg=0,nregs
        natrg(ireg) = 0
299     acttmp(ireg) = 0.
ccdir$ novector
        do 300 i=1,natoms
        ireg = ipntrg(i)
        natrg(ireg) = natrg(ireg) + 1
300     acttmp(ireg) = acttmp(ireg) +
     .     amass(itype(i))*(itempg(ireg,1)*rv(4,i)**2
     .    + itempg(ireg,2)*rv(5,i)**2 + itempg(ireg,3)*rv(6,i)**2)
        do 305 ireg=0,nregs
	three=itempg(ireg,1)+itempg(ireg,2)+itempg(ireg,3)
        if(natrg(ireg).ne.0)acttmp(ireg) = acttmp(ireg)/
     .      (three*boltz*float(natrg(ireg)))
305     continue
c
c  add drag term to force
c
        do 306 ireg=1,nregs
        drag(ireg) = (acttmp(ireg)-tmpreg(ireg))/qreg(ireg)
306     continue
        do 307 i=1,natoms
        ireg = ipntrg(i)
        do 307 j = 1,3
307     f(j,i) = f(j,i) - drag(ireg)*rv(j+3,i)*itempg(ireg,j)
ccdir$ vector
c
      elseif(ifxtmp.eq.2)then
c
c  hoover drag
c
c
c  calculate the temperature of each region
c
ccm=0
csumf=0.
        do 399 ireg=0,nregs
        natrg(ireg) = 0
399     acttmp(ireg) = 0.
ccdir$ novector
        do 400 i=1,natoms
        ireg = ipntrg(i)
        natrg(ireg) = natrg(ireg) + 1
ccm=cm+amass(itype(i))*rv(4,i)
csumf=sumf+f(1,i)
cprint *,' in hoover',i, cm,sumf
400     acttmp(ireg) = acttmp(ireg) +
     .     amass(itype(i))*(itempg(ireg,1)*rv(4,i)**2
     .    + itempg(ireg,2)*rv(5,i)**2 + itempg(ireg,3)*rv(6,i)**2)
        do 405 ireg=0,nregs
        three=itempg(ireg,1)+itempg(ireg,2)+itempg(ireg,3)
        if(natrg(ireg).ne.0)acttmp(ireg) = acttmp(ireg)/
     .      (three*boltz*float(natrg(ireg)))
405     continue
c
c  add drag term to force
c
        do 420 ireg=1,nregs
        dragdot = (acttmp(ireg)-tmpreg(ireg))/qreg(ireg)
cprint *,'deltat',deltat
        drag(ireg) = drag(ireg) + dragdot*deltat
420     continue
        do 421 i=1,natoms
        ireg = ipntrg(i)
        do 421 j = 1,3
421     f(j,i) = f(j,i) - drag(ireg)*rv(j+3,i)*itempg(ireg,j)
ccdir$ vector
      elseif( (ifxtmp.eq.3).and.(deltat.ne.0.0) ) then
c
c  langevin
c
c  note:  do not calculate the temperature
c
ccdir$ novector
        do 530 ireg=1,nregs
530     drag(ireg) = qreg(ireg)
        do 540 i=1,natoms
        ireg = ipntrg(i)
        do 540 j=1,3
	if(itempg(ireg,j).eq.0) go to 540
        stoch = sqrt( 2.*drag(ireg)*boltz*tmpreg(ireg)/
     .                deltat ) * rgauss()
        f(j,i) = f(j,i) - 
     1       (drag(ireg)*rv(j+3,i) - stoch)*itempg(ireg,j)
540     continue
ccdir$ vector
	elseif( ifxtmp.eq.4) then
c	Hugoniostat
	ireg=1
	acte=0.0
	do 522 i=1,natoms
	acte=acte+e(i)
	do 522 j=1,3
	acte=acte+0.5*amass(itype(i))*rv(j+3,i)**2
522	continue
	acte=acte/natoms
	actv=perlen(1)*perlen(2)*perlen(3)/natoms
	actpzz=stresst(3,3)/actv/natoms
        dragdot = (acte-e0-0.5*(actpzz+pzz)*(v0-actv))*qreg(ireg)
        drag(ireg) = drag(ireg) + dragdot*deltat
cprint *,'in fixtmp',dragdot,drag(ireg)
c
c  add drag term to force
c
	do 520 i=1,natoms
	do 520 j=1,3
	f(j,i) = f(j,i) - drag(ireg)*rv(j+3,i)*amass(itype(i))
520	continue
	
      endif
c
c  update last time
      timetmp = time
c
      return
      end
c*******************************************************************
c
c
      subroutine scale_sub
      include 'param.inc'
      include 'utils.inc'
c
c  determine the type of variable scaling used
c
      goto (1000,2000) ibdtyp
      write(6,9001) ibdtyp
9001  format(1x,'undefined boundary type:  ibdtyp=',i3)
      stop
c
1000  continue
c
c  ibdtyp = 1
c    fixed boundaries  no scaling
c
      do 1100 j = 1,6
      do 1100 i = 1,natoms
1100  y(j,i) = rv(j,i)
      return
c
c  ibdtyp = 2
c    dynamic periodic lengths  rectangular orientation
c
2000  continue
      do 2100 j = 1,3
	if (idynper(j).eq.0) then
	scfact=1.0
        origin =  0.0
	else
	scfact = 1.0/perlen(j)
        origin = perlb(j)
	endif
      do 2100 i = 1,natoms
      y(j,i) = (rv(j,i) - origin)*scfact
      y(3+j,i) = rv(3+j,i)*scfact
2100  continue
      do 2200 j = 1,3
      y(j,natoms+1) = perlen(j)
      y(j+3,natoms+1) = bndvel(j)
2200  continue
      return
      end
      subroutine uscale
      include 'param.inc'
      include 'utils.inc'
c
c  determine the type of variable scaling used
c
      goto (1000,2000) ibdtyp
      write(6,9001) ibdtyp
9001  format(1x,'undefined boundary type:  ibdtyp=',i3)
      stop
c
1000  continue
c
c  ibdtyp = 1
c    no scaling
c
      do 1100 j = 1,6
      do 1100 i = 1,natoms
1100  rv(j,i) = y(j,i)
      return
c
c  ibdtyp = 2
c    dynamic periodic lengths  rectangular orientation
c
2000  continue
      do 2100 j = 1,3
      perlen(j) = y(j,natoms+1)
      perub(j) = perlb(j) + y(j,natoms+1)
      bndvel(j) = y(j+3,natoms+1)
2100  continue
      do 2200 j = 1,3
        if (idynper(j).eq.0) then
        scfact=1.0
        origin = 0.0
        else
	scfact = perlen(j)
        origin = perlb(j)
        endif
      do 2200 i = 1,natoms
      rv(j,i) = y(j,i)*scfact + origin
      rv(3+j,i) = y(3+j,i)*scfact
2200  continue
      return
      end
c********************************************************************
c
c  this subroutine puts all the particles back into the basic
c  periodic cell.
c
      subroutine colect
      include 'param.inc'
      include 'utils.inc'
c
c  determine the type of variable scaling used
c
      goto (1000,2000) ibdtyp
      write(6,9001) ibdtyp
9001  format(1x,'undefined scaling method ibdtyp=',i3)
      stop
c
1000  continue
c
c  ibdtyp = 1
c    no scaling
c
      do 1100 i = 1,natoms
      if (y(2,i).gt.perub(2)) y(2,i) = y(2,i) - perlen(2)
      if (y(2,i).lt.perlb(2)) y(2,i) = y(2,i) + perlen(2)
      if (y(3,i).gt.perub(3)) then
	y(3,i) = y(3,i) - perlen(3)
	y(1,i) = y(1,i) - ctbeta*perlen(3)
	endif
      if (y(3,i).lt.perlb(3)) then
	y(3,i) = y(3,i) + perlen(3)
	y(1,i) = y(1,i) + ctbeta*perlen(3)
	endif
      if (y(1,i).gt.perub(1)) y(1,i) = y(1,i) - perlen(1)
      if (y(1,i).lt.perlb(1)) y(1,i) = y(1,i) + perlen(1)
1100  continue
      return
c
c  ibdtyp = 2
c    dynamic periodic lengths 
c
2000  continue
      do 2100 j = 1,3
	if(idynper(j).ne.0) then
	pub=1.0
	plb=0.0
	scfact=1.0
	ct=ctbeta*perlen(3)/perlen(1)
	else
	pub=perub(j)
	plb=perlb(j)
	scfact=perlen(j)
	ct=ctbeta*perlen(3)
	endif
      do 2100 i = 1,natoms
     	if (y(j,i).gt.pub) then
	y(j,i) = y(j,i) - scfact
	if(j.eq.3) y(1,i)=y(1,i)-ct
	endif
     	if (y(j,i).lt.plb) then
	y(j,i) = y(j,i) + scfact
	if(j.eq.3) y(1,i)=y(1,i)+ct
	endif
2100  continue
c     do 2200 j = 1,3
c     y(j,natoms+1) = perlen(j)
c     y(j+3,natoms+1) = bndvel(j)
2200  continue
      return
      end
c*******************************************************************
      subroutine modlat
      include 'param.inc'
      include 'utils.inc'
      return
      end
c*******************************************************************
      subroutine modvel
      include 'param.inc'
      include 'utils.inc'
      return
      end
c*******************************************************************
      subroutine modint
      include 'param.inc'
      include 'utils.inc'
      return
      end
c*******************************************************************
      subroutine modfor
      include 'param.inc'
      include 'utils.inc'
      return
      end
c************************************************************************
c
c  this routine computes the instantaneous energies, temperature
c  and pressure for the current configuration.
c
      subroutine calce(ekin,pot,etot,tempx,press,iclfor,time)
      include 'param.inc'
      include 'utils.inc'
	common/switch1/swval,potref,eref(natmax)
	common/free/ifree
      data boltz/8.617e-5/
c
c
c  compute the kinetic energy
c
cprint *,'entering calce ',ifree,iclfor
      ekin=0.0
      do 20 i=1,natoms
      v2 = rv(4,i)**2 + rv(5,i)**2 + rv(6,i)**2
 20   ekin = ekin + amass(itype(i))*v2
      ekin = ekin*0.5
c
c  compute the temperature
c
      tempx = 7.7366e+3*ekin/float(natoms)
c
c  compute the occupation and temperature of the regions
c
      if(ifxtmp.ne.0)then
        do 40 ireg=0,nregs
        natrg(ireg) = 0
40      acttmp(ireg) = 0.
ccdir$ novector
        do 50 i=1,natoms
        ireg = ipntrg(i)
        natrg(ireg) = natrg(ireg) + 1
50      acttmp(ireg) = acttmp(ireg) +
     .     amass(itype(i))*(itempg(ireg,1)*rv(4,i)**2 
     .   + itempg(ireg,2)*rv(5,i)**2 + itempg(ireg,3)*rv(6,i)**2)
        do 60 ireg=0,nregs
	three=itempg(ireg,1)+itempg(ireg,2)+itempg(ireg,3)
        if(natrg(ireg).ne.0)acttmp(ireg) = acttmp(ireg)/
     .     (three*boltz*float(natrg(ireg)))
60      continue
ccdir$ vector
      endif
c
c  call force to obtain the potential energy and pressure information
c
	if(ifree.eq.0) then
         if (iclfor.eq.1)  call force(time)
cprint *,'called force ', e(1)
c
c  compute the potential energy and pressure
c
         pot=0.0
         do 300 i=1,natoms
300      pot=pot+e(i)
	 pot=pot/natoms
	 ekin=ekin/natoms
         etot=ekin+pot
cprint *,'called force ', e(1)
	else
         if (iclfor.eq.1)  call xforce(time)
         potref=0.0
         do 301 i=1,natoms
	 pot=pot+e(i)
301      potref=potref+eref(i)
	 pot=pot/natoms
	 potref=potref/natoms
	 ekin=ekin/natoms
	 etot=ekin+(1.0d0-swval)*pot+swval*potref
	endif
      volume = perlen(1)*perlen(2)*perlen(3)
      press = (stresst(1,1)+stresst(2,2)+stresst(3,3))/(3.0*volume)
c
c  convert the pressure to bar
c
      press = 1.602e6*press
cprint *,'in calce ', natoms,pot,etot,press,swval
c  all energies  are per atom
      return
      end
c*******************************************************************
c
c  initlat initializes the positions and velocities of the atoms
c
      subroutine initlat(contin,lastconf)
      include 'param.inc'
      include 'utils.inc'
      logical contin,lastconf
      logical rstart,genlat,genvel,gendis,sort
      character*80 initf
      save initf
      dimension scale(3),origin(3),scalev(3)
        data beta,betar,cbeta,sbeta,ctbeta/90.0,-9999.,0.0,1.0,0.0/
	data shear,sheardot/2*0.0/
	data pi/3.1415926/
      data rstart,genlat,genvel,gendis,sort/5*.false./
      data scale/3*-9999./
      data scalev/3*-9999./
      data dismag/0.1/
      namelist /initcard/ genlat,genvel,gendis,dismag,sort,initf,scale,
     1     beta,shear,scalev,sheardot
c
c  determine how the program is to be initialized
c
	pi=2*acos(0.d0)
      if (.not.contin) read(5,initcard)
c
c     read in the start-up file (22) from a previous run if available
c     initf of 'none' indicates that no restart file is available
c
      if (contin.and.lastconf) then
         genvel = .false.
         genlat = .false.
         sort = .false.
         t0 = tend
       else
         if (initf(1:4).eq.'NONE'.or.initf(1:4).eq.'none') then
c
c  if no restart file, generate position and velocity and also sort
c
            sort = .true.
            genvel = .true.
            genlat = .true.
          else
            rstart = .true.
            open(unit=22,file=initf,form='FORMATTED',status='OLD')
            write(6,9135)initf
9135        format(/,/,1x,'******   using restart',/,
     1         ' reading in initial file named ',a80)
            call restart(rstart)
	    tend=t0
            close(22)
          end if
       endif
c
c  generate a lattice if requested
c    if lattice generation not requested, read in latcard and ignore
c    its contents
c
      if (genlat) then
         call latgen(contin)
       else
         call latdum(contin)
       end if
        beta=beta*pi/180.
        if(shear.ne.0.0) then
c ignore beta if shear specified
        cbeta=cos(betar*pi/180.)-2*shear
        beta=acos(cbeta)
        elseif (betar.gt.90.) then
c ignore beta if restart file has beta not equal to 90
        beta=pi/180.*betar
        endif
        cbeta=cos(beta)
        sbeta=sin(beta)
        ctbeta=cbeta/sbeta
        beta=180./pi*beta
c
c  check that the number of atoms is not too large
c
      if (natoms.gt.natmax) then
         write(6,9160) natoms,natmax
 9160    format(' natoms=',i10,' is greater than natmax=',i10)
         stop
       end if
c
c  sort lattice if requested
c
      if (sort) call sorter
c
c  scale the lattice if requested
c
         if (scale(1).eq.-9999.) scale(1) = 1.
         if (scale(2).eq.-9999.) scale(2) = scale(1)
         if (scale(3).eq.-9999.) scale(3) = scale(1)
      scale(3)=scale(3)*sbeta/sin(betar*pi/180.)
            write(6,9010) (scale(j),j=1,3),beta
 9010       format(1x,'***************',/,'  scaling lattice by ',3g13.5,
     1           /,'***************',/,'  shear angle ',f10.5)
c
c  scale the velocities if requested
c
         if (scalev(1).eq.-9999.) scalev(1) = 1.
         if (scalev(2).eq.-9999.) scalev(2) = scalev(1)
         if (scalev(3).eq.-9999.) scalev(3) = scalev(1)
	 write(6,9011) scalev
 9011	 format(1x,'***************',/,'  scaling velocities by ',3g13.5,
     1           /,'***************')
        call colect
        br=pi/180*betar
        factor=ctbeta-cos(br)/sbeta
         do 100 i = 1,3
         perlen(i) = scale(i)*(perub(i)-perlb(i))
         origin(i) = 0.5*(perlb(i)+perub(i))
         perub(i) = origin(i) + 0.5*perlen(i)
         perlb(i) = perub(i) - perlen(i)
         do 100 j = 1,natoms
         rv(i,j) = origin(i) + scale(i)*(rv(i,j)-origin(i))
        if(i.eq.3) rv(1,j)=rv(1,j)+rv(3,j)*factor
	 rv(i+3,j)=scalev(i)*rv(i+3,j)
100      continue
            
c      call modlat
c
c  check periodicity
c
      call chkper
c
c  generate random displacements if requested
c
      if (gendis) call disgen(dismag)
c
c  generate the velocities if requested
c    if velocities not requested, read velcard and ignore its contents
c
      if (genvel) then
         call velgen(contin)
       else
         call veldum(contin)
       end if
c      call modvel
c
c  check that all types have been defined
c
      call chktyp
c
      return
      end
c****************************************************************************
c
c  this subroutine reads in the parameters required to restart
c  a run from fortran file 22
c
      subroutine restart(rstart)
      include 'param.inc'
      include 'utils.inc'
      logical rstart
	character *80 header
c
c  note that the old header is not passed to the rest of the
c  program
c
      read(22,9501,end=8000) header
9501  format(a)
      write(6,9001)header
9001  format(/,' previous header:',/,3x,a,/)
      rstart=.true.
      read(22,9502) natoms,ntyp,t0,betar
9502  format(2i10,e15.8,f10.5)
        if(betar.eq.0) then
        betar=90.
        else
        write(6,9511) betar
9511    format(' ***monoclinic structure*** using beta',
     1     ' from restart file ',f10.5)
        endif
c
c  check that the number of atoms is not too large
c
      if (natoms.gt.natmax) then
         write(6,9901) natoms,natmax
9901     format(1x,'natoms=',i6,' is greater than natmax=',i6)
         stop
       end if
c
      write(6,9002)natoms,ntyp
9002  format(1x,i10,' atoms',i10,' particle types')
      write(6,9003) t0
9003  format(1x,'final time of previous run: ',f10.3)
      read(22,9503) (perub(i),i=1,3),(perlb(i),i=1,3)
9503  format(3e25.16)
      read(22,9504) (amass(i),ielement(i),i=1,ntyp)
9504  format(e25.16,i10)
      read(22,9505) ((rv(i,j),i=1,6),itype(j),j=1,natoms)
9505  format(3e25.16/3e25.16/i10)
      read(22,9503,end=200) (bndvel(i),i=1,3)
	write(6,9506) bndvel
9506	format('boundary velocities = ',3f15.8)
      goto 1000
200   continue
      do 210 i = 1,3
210   bndvel(i) = 0.0
      write(6,9201)
9201  format(1x,'no bndvel available, zero assumed')
c
c       print out types
c
1000  continue
      write(6,9100)
 9100 format('  restart uses these types: ')
      write(6,9102)
 9102 format('   type  element      amass  ',/,
     1       '   ----  -------    ---------')
      write(6,9103)(i,ielement(i),amass(i),i=1,ntypes)
 9103 format(1x,i4,i9,2x,g11.4)
c
c  compute perlen
c
      do 100 i = 1,3
100   perlen(i) = perub(i) - perlb(i)
      write(6,9031)
 9031 format('  periodicity   ',15x,'x',14x,'y',14x,'z')
      write(6,9032)perlb
 9032 format('  lower periodic bound  ',3g15.5)
      write(6,9033)perub
 9033 format('  upper periodic bound  ',3g15.5)
      write(6,9034)perlen
 9034 format('  length        ',8x,3g15.5)
c
 8000 return
      end
c**************************************************************************
c
c  this routine sets up the gradient of the energy for va08a
c
      subroutine vafunc(n,ener)
      include 'param.inc'
      include 'utils.inc'
      dimension save(3)
      common /mincom/ grad(3*natmax),s(3*natmax)
c
c
c  unscale the coordinates (compute physical coordinates)
c
      call colect
      call uscale
c
c  determine the number of periodic images needed if the periodic
c  lengths are dynamic
c
      if (ibdtyp.ne.1) call chkper
c
c  compute the forces for this configuration
c
      call force(t0)
c
c  determine the energy and the gradients
c
      if (ibdtyp.eq.1) then
         ener = 0.0
         do 10 i = 1,natoms
         ener = ener + e(i)
         grad(1+3*(i-1)) = -f(1,i)
         grad(2+3*(i-1)) = -f(2,i)
         grad(3+3*(i-1)) = -f(3,i)
10       continue
       else
         ener = 0.0
         do 20 i = 1,natoms
         ener = ener + e(i)
         do 20 j=1,3
         if(idynper(j).eq.0) then
         grad(j+3*(i-1)) = -f(j,i)
         else
         grad(j+3*(i-1)) = -perlen(j)*f(j,i)
         endif
20       continue
         ener = ener + e(natoms+1)
         volume = perlen(1)*perlen(2)*perlen(3)
	do 21 k=1,3
	i1=mod(k+1,3)+1
	i2=mod(k+3,3)+1
         save(k) =((dpress-stresst(k,k)/volume)*perlen(i1)*perlen(i2)
     2               + perlen(k)*dstress(k))/perlen(k)
21	continue
        if(iper.eq.4) then
         sav=(save(1)+save(2)+save(3))/3
         grad(3*natoms+1)=sav*perlen(1)
         grad(3*natoms+2)=sav*perlen(2)
         grad(3*natoms+3)=sav*perlen(3)
        elseif (iper.eq.3) then
         sav=(save(1)+save(2))/2
         grad(3*natoms+1)=sav*perlen(1)
         grad(3*natoms+2)=sav*perlen(2)
         grad(3*natoms+3)=save(3)*abs(idynper(3))*perlen(3)
        elseif (iper.eq.2) then
         sav=(save(1)+save(3))/2
         grad(3*natoms+1)=sav*perlen(1)
         grad(3*natoms+3)=sav*perlen(3)
         grad(3*natoms+2)=save(2)*abs(idynper(2))*perlen(2)
        elseif (iper.eq.1) then
         sav=(save(2)+save(3))/2
         grad(3*natoms+2)=sav*perlen(2)
         grad(3*natoms+3)=sav*perlen(3)
         grad(3*natoms+1)=save(1)*abs(idynper(1))*perlen(1)
        else
         grad(3*natoms+1)=save(1)*abs(idynper(1))*perlen(1)
         grad(3*natoms+2)=save(2)*abs(idynper(2))*perlen(2)
         grad(3*natoms+3)=save(3)*abs(idynper(3))*perlen(3)
        endif
       endif
      return
      end
c**************************************************************************
c
c  this routine integrates newton's equations using the slatec
c  integration routine deabm
c
      subroutine slatint
        implicit real*8 (a-h,o-z)
      print *,'slatint not implemented'
      stop
      end
c***************************************************************************
c
c  this routine calls the nordsieck integrator
c
      subroutine nordint
      include 'param.inc'
      include 'utils.inc'
	logical contin
c
c
      t=t0
      write(6,9001)dt
 9001 format(' calling integrator (nordsieck):  dt =',e15.5)
      x1 = seconds()
c
c  nords produces an output file every noutp time steps
c
      otime = dt*noutp
	do 777 ii=1,ngen
	contin=.true.
	if(ii.gt.1) call velgen(contin)
      call nords(t,ndegfr,otime)
777	continue
c
c  ensure that all the particles are inside the periodic box at
c  end of the integration
c
      call colect
c
      x2 = seconds() - x1
      write(6,9101)x2
 9101 format(' time for call to integrator ',e15.5)
      write(6,9102)nforce
 9102 format(' number of calls to subroutine force  ',i20)
      return
      end
c**************************************************************************
c
c  this routine implements the nodsieck integration scheme
c   (see j.r. beeler, "radiation effects computer experiments",
c    p. 100 (1983))
c
      subroutine nords(t,n,otime)
      include 'param.inc'
      include 'utils.inc'
      dimension u1(3,natmax),u2(3,natmax),u3(3,natmax),u4(3,natmax),
     .          u5(3,natmax),phi(3,natmax)
	dimension tmass(0:6),du1(0:6,3),avu2(0:6,3),
     1            avu3(0:6,3),avu4(0:6,3),avu5(0:6,3)
      data c0,c1,c2,c3,c4,c5/0.1875,0.697222222,1.,0.61111111111,
     1                      0.16666666667,0.016666666667/
c
      dt2 = 0.5*dt**2
      dtrec = 1.0/dt
      dt1 = dt
c
c  initialize the higher derivatives to zero at the start
c
      do 100 j = 1,3
      do 100 i = 1,n
      u3(j,i) = 0.0
      u4(j,i) = 0.0
      u5(j,i) = 0.0
100   continue
c
c  set the velocity and acceleration variables
c
      call accel(t)
      do 200 j = 1,3
      do 200 i = 1,n
      u1(j,i) = dt*y(j+3,i)
      u2(j,i) = dt2*acc(j,i)
200   continue
c
c  the program is now ready to integrate the equations of motion
c
c  first initialize the timing value
c
      tinit = t
      tfinal = tinit + nsteps*dt - dt/100.
      dt1 = dt
      timstr = seconds()
1000  continue
      t = t + dt1
c
c  first compute the estimated values at the next time step.
c
      do 2000 j = 1,3
      do 2000 i = 1,n
      y(j,i) = y(j,i) + u1(j,i) + u2(j,i) + u3(j,i) +
     $         u4(j,i) + u5(j,i)
      u1(j,i) = u1(j,i) + 2.0*u2(j,i) + 3.0*u3(j,i) +
     $           4.0*u4(j,i) + 5.0*u5(j,i)
      y(j+3,i) = dtrec*u1(j,i)
      u2(j,i) = u2(j,i) + 3.0*u3(j,i) + 6.0*u4(j,i) + 10.*u5(j,i)
      u3(j,i) = u3(j,i) + 4.0*u4(j,i) + 10.0*u5(j,i)
      u4(j,i) = u4(j,i) + 5.0*u5(j,i)
2000  continue
c
c  compute the actual accelerations at this point
c
      call accel(t)
c
c  Evaluate the ACFS if needed.
c
      if(ncalcacf.EQ.0)then
        dens=natoms/(perlen(1)*perlen(2)*perlen(3))
       call evalvacf(natoms,rv,acc,slocal,dt,amass,itype,
     $  nbuff,nval,ntypes,netype,nlimitacf,dens,e)
       endif
c
c  determine the maximum acceleration on any particle
c
      amax = 0.
c      do 2100 i = 1,n
c      amax1 = max(amax,acc(1,i)**2 + acc(2,i)**2 + acc(3,i)**2)
c2100  continue
c      amax1 = dt1*sqrt(amax)
c
c  check if the largest acceleration is too large,
c   if so, cut time step in half
c
      if (amax.lt.accbnd(1)) goto 2900
c
c  undo the previous estimates for the next time step
c
      do 2300 j = 1,3
      do 2300 i = 1,n
      u4(j,i) = u4(j,i) - 5.0*u5(j,i)
      u3(j,i) = u3(j,i) - 4.0*u4(j,i) - 10.0*u5(j,i)
      u2(j,i) = u2(j,i) - 3.0*u3(j,i) - 6.0*u4(j,i) - 10.0*u5(j,i)
      u1(j,i) = u1(j,i) - 2.0*u2(j,i) - 3.0*u3(j,i) - 4.0*u4(j,i)
     1        - 5.0*u5(j,i)
      y(j+3,i) = u1(j,i)*dtrec
      y(j,i) = y(j,i) - u1(j,i) - u2(j,i) - u3(j,i) - u4(j,i) - u5(j,i)
2300  continue
c
c  change the dt values
c
      t = t - dt1
      dt1 = dt1/2.
      dt2 = dt2/4.
      dtrec = 2.*dtrec
c
c  rescale the u arrays
c
      do 2400 j = 1,3
      do 2400 i = 1,n
      u1(j,i) = u1(j,i)/2.
      u2(j,i) = u2(j,i)/4.
      u3(j,i) = u3(j,i)/8.
      u4(j,i) = u4(j,i)/16.
      u5(j,i) = u5(j,i)/32.
2400  continue
      write(6,9291) dt1,t
9291  format(1x,'** time step changed to ',g12.5,' at t =',f8.4)
      goto 1000
2900  continue
c
c  now compute the displacement function phi
c
      do 3000 i = 1,n
      do 3000 j = 1,3
      phi(j,i) = dt2*acc(j,i) - u2(j,i)
3000  continue
c
c  now correct the predicted values of the functions
c
      do 4000 j = 1,3
      do 4000 i = 1,n
      y(j,i) = y(j,i) + c0*phi(j,i)
      u1(j,i) = u1(j,i) + c1*phi(j,i)
      y(j+3,i) = dtrec*u1(j,i)
      u2(j,i) = u2(j,i) + c2*phi(j,i)
      u3(j,i) = u3(j,i) + c3*phi(j,i)
      u4(j,i) = u4(j,i) + c4*phi(j,i)
      u5(j,i) = u5(j,i) + c5*phi(j,i)
4000  continue
c set CM velocity if required
       if(icm) then
      do 299  ireg=0,nregs
        tmass(ireg) = 0.0
      do 109 jc=1,3
        actforce(ireg,jc) = 0.0
	avu2(ireg,jc)= 0.0
	avu3(ireg,jc)= 0.0
	avu4(ireg,jc)= 0.0
	avu5(ireg,jc)= 0.0
109     actvcm(ireg,jc)= 0.0
299   continue
      do 110 i=1,natoms
	ireg=ipntrg(i)
        it = itype(i)
        tmass(ireg) = tmass(ireg) + amass(it)
      do 105 jc=1,3
        actforce(ireg,jc) = actforce(ireg,jc) + amass(it)*acc(jc,i)
        actvcm(ireg,jc) = actvcm(ireg,jc) + amass(it)*y(jc+3,i)
	avu2(ireg,jc)=avu2(ireg,jc)+u2(jc,i)
	avu3(ireg,jc)=avu3(ireg,jc)+u3(jc,i)
	avu4(ireg,jc)=avu4(ireg,jc)+u4(jc,i)
	avu5(ireg,jc)=avu5(ireg,jc)+u5(jc,i)
105   continue
110   continue
      do 298 ireg=0,nregs
	if(tmass(ireg).gt.0.0) then
      do 201 jc=1,3
         actvcm(ireg,jc)=actvcm(ireg,jc)/tmass(ireg)
         actforce(ireg,jc)=actforce(ireg,jc)/tmass(ireg)
	 du1(ireg,jc)=actvcm(ireg,jc)-vcmrg(ireg,jc)
	 avu2(ireg,jc)=avu2(ireg,jc)/natrg(ireg)
	 avu3(ireg,jc)=avu3(ireg,jc)/natrg(ireg)
	 avu4(ireg,jc)=avu4(ireg,jc)/natrg(ireg)
	 avu5(ireg,jc)=avu5(ireg,jc)/natrg(ireg)
201   continue
 	endif
298   continue
      do 300 i=1,natoms
	ireg=ipntrg(i)
      do 301 jc=1,3
	u2(jc,i)=u2(jc,i)-avu2(ireg,jc)
	u3(jc,i)=u3(jc,i)-avu3(ireg,jc)
	u4(jc,i)=u4(jc,i)-avu4(ireg,jc)
	u5(jc,i)=u5(jc,i)-avu5(ireg,jc)
301     u1(jc,i)=u1(jc,i)-du1(ireg,jc)*dt
300     continue
       endif
	call average(1,t)
c
c  output the results if appropriate
c
      if (abs(mod(t-tinit+0.5*dt1,otime)-0.5*dt1).lt.(0.5*dt1)) then
         call output(t,iaccur)
         tend = t
c
c  determine if there is enough time left to go to next output
c
         timend = seconds()
         timned = 1.2*timend
         call trmain(timlft)
         if (timlft.le.timned) then
            write(6,9201)
9201        format(1x,'job terminated early due to time limit')
            return
          end if
         timstr = timend
       end if
c
c  check to see if the time step can now be increased
c
      if (dt1.lt.dt.and.amax.lt.accbnd(2).and.
     1    (t-tinit-nint((t-tinit)/dt)*dt).le..001*dt) then
         dt1 = 2.*dt1
         dt2 = 4.*dt2
         dtrec = dtrec/2.
         do 5000 j = 1,3
         do 5000 i = 1,n
         u1(j,i) = 2.*u1(j,i)
         u2(j,i) = 4.*u2(j,i)
         u3(j,i) = 8.*u3(j,i)
         u4(j,i) = 16.*u4(j,i)
         u5(j,i) = 32.*u5(j,i)
5000     continue
         write(6,9591) dt1,t
9591     format(1x,'** time step changed to ',g12.5,' at t =',f8.4)
       endif
c
c  see if done with the integration
c
      if (t.lt.tfinal) goto 1000
      return
      end
c*******************************************************************
c
c  this routine computes the accleration on the potentially scaled
c  variables describing the particle positions and optionally the
c  periodic boundaries
c
      subroutine accel(time)
      include 'param.inc'
      include 'utils.inc'
        common /free/ ifree
      common /mincom/ grad(3*natmax),s(3*natmax)
        dimension save(3)
c
c  unscale the coordinates (compute physical coordinates)
c
      call colect
      call uscale
c
c  determine the number of periodic images needed if the periodic
c  lengths are dynamic
c
      if (ibdtyp.ne.1) call chkper
c
c  compute the forces for this configuration
c
	if(ifree.eq.0) then
      call force(time)
	else
      call xforce(time)
	endif
ccall scale_sub
c
c  convert the forces into accelerations on the scaled coordinates
c
c  branch to the desired boundary type
c
      goto (1000,2000) ibdtyp
      write(6,9010) ibdtyp
9010  format(1x,'undefined boundary type: ibdtyp =',i3)
      stop
c
c  ibdtyp = 1
c    fixed boundaries, no scaling of variables
c
1000  continue
      do 1100 j = 1,natoms
      do 1100 i = 1,3
1100  acc(i,j) = f(i,j)/amass(itype(j))
cprint *, 'in accel', f(1,1),acc(1,1)
      return
c
c  ibdtyp = 2
c    lengths of x, y and z sides of the periodic cell dynamic
c    (remains orthogonal)
c       (i.e., transform forces to scale coordinates)
c
2000  continue
      volume = y(1,ndegfr)*y(2,ndegfr)*y(3,ndegfr)
      do 2200 i = 1,3
      save(i) = (((stresst(i,i)-dpress*volume)/y(i,ndegfr)
     2 - y(i,ndegfr)*dstress(i))/bndmas(i)
     3 - bnddrg(i)*y(3+i,ndegfr))/perlen(i)
2200  continue
        if(iper.eq.4) then
         sav=(save(1)+save(2)+save(3))/3
         acc(1,ndegfr)=sav*perlen(1)
         acc(2,ndegfr)=sav*perlen(2)
         acc(3,ndegfr)=sav*perlen(3)
        elseif (iper.eq.3) then
         sav=(save(1)+save(2))/2
         acc(1,ndegfr)=sav*perlen(1)
         acc(2,ndegfr)=sav*perlen(2)
         acc(3,ndegfr)=save(3)*abs(idynper(3))*perlen(3)
        elseif (iper.eq.2) then
         sav=(save(1)+save(3))/2
         acc(1,ndegfr)=sav*perlen(1)
         acc(3,ndegfr)=sav*perlen(3)
         acc(2,ndegfr)=save(2)*abs(idynper(2))*perlen(2)
        elseif (iper.eq.1) then
         sav=(save(2)+save(3))/2
         acc(2,ndegfr)=sav*perlen(2)
         acc(3,ndegfr)=sav*perlen(3)
         acc(1,ndegfr)=save(1)*abs(idynper(1))*perlen(1)
        else
         acc(1,ndegfr)=save(1)*abs(idynper(1))*perlen(1)
         acc(2,ndegfr)=save(2)*abs(idynper(2))*perlen(2)
         acc(3,ndegfr)=save(3)*abs(idynper(3))*perlen(3)
        endif
	do 2201 k=1,3
	if(idynper(k).eq.2) acc(k,ndegfr)=0.0
2201	continue
      do 2100 j = 1,natoms
      do 2100 i = 1,3
        if(idynper(i).eq.0) then
        acc(i,j) = f(i,j)/amass(itype(j))
        else
c     acc(i,j) = (f(i,j)/amass(itype(j))-y(i,ndegfr)*acc(1,ndegfr)
      acc(i,j) = (f(i,j)/amass(itype(j))
     1         - 2.*y(i+3,ndegfr)*y(i+3,j))/y(i,ndegfr)
        endif
2100  continue
      return
      end
c**********************************************************************
c
c  nm3add adds atom i to the neighbor list of atom j
c
      subroutine nm3add(i,j)
      include 'param.inc'
      include 'utils.inc'
c
c  increment the number of neighbors for atom j
c
      nnindx(j) = nnindx(j) + 1
c  check that the number of neighbors has not exceeded neimax
c
      if (nnindx(j).gt.neimax) then
         write(6,9901) nnindx(j),j
9901     format(1x,'number of neighbors,',i5,'for atom,',i5)
         stop
      endif
c
c  skip over the elements in the list before i
c
      k = 0
100   continue
      k = k + 1
      if (nlist3(k,j).lt.i.and.k.lt.nnindx(j)) goto 100
c
c  shift the rest of the list up
c
      do 200 l = nnindx(j),k+1,-1
200   nlist3(l,j) = nlist3(l-1,j)
c
c  add the atom to the list
c
      nlist3(k,j) = i
      return
      end
c**********************************************************************
c
c  nm3del deletes atom i to the neighbor list of atom j
c
      subroutine nm3del(i,j)
      include 'param.inc'
      include 'utils.inc'
c
c  decrement the number of neighbors for atom j
c
      nnindx(j) = nnindx(j) - 1
c
c  skip over the elements in the list before i
c
      k = 0
100   continue
      k = k + 1
      if (nlist3(k,j).ne.i) goto 100
c
c  shift the rest of the list down
c
      do 200 l = k,nnindx(j)
200   nlist3(l,j) = nlist3(l+1,j)
      return
      end
c************************************************************************
c
c  this routine determines the neighbors of a given atom
c  for full=.t. all neighbors are obtained
c  for full=.f. only neighbors with index less than i are obtained
c
      subroutine gneigh(i,full)
      include 'param.inc'
      include 'utils.inc'
      include 'wrkspc.inc'
	logical full
      dimension per12(3)
c
c  determine the timing
c
c      call timeused(icpu,io,isys,imem)
c      timin = (1.e-6)*float(icpu)
	cmono=ctbeta*perlen(3)
c
c  define the constants needed to find the nearest periodic image
c
      do 10 kcoord = 1,3
10    per12(kcoord) = 0.50*perlen(kcoord)
      ity = itype(i)
c
c  branch to the appropriate neighbor finding method
c
      goto (1000,2000,3000,90) abs(nmeth)
90    write(6,9001) nmeth
9001  format(1x,'undefined neighbor method requested; nmeth=',i3)
      stop
c
1000  continue
c
c  nmeth = 1
c    this is the order n**2 method when there is no scaling of the
c    variables
c
c       first do off diagonal terms (i.ne.j)
c       fortran 77 convention: if i=1 then j=1,i-1 loop is skipped
c
	if(full) then
	jmax=natoms
	else
	jmax=i-1
	endif
      do 1100 j = 1,jmax
c
c     compute the square of the distance to the closest periodic image
c
	call rdx3(rv(1,i),rv(1,j),per12,perlen,cmono,dis(1,j),r(j))
1100  continue
c
c     determine which pairs are separated by less than rcut
c     and store the needed information about these pairs
c
      nneigh = 0

      do 1200 j = 1,jmax
c
c     if nearest periodic image is out of range, then all images
c     will be
c
      if (r(j).gt.rctsqn.or.j.eq.i) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = r(j)
      jneigh(nneigh) = j
      do 1250 kcoord = 1,3
1250  dneigh(kcoord,nneigh) = dis(kcoord,j)
      fpn(nneigh) = fp(j)
c     check periodic images in z direction
      if(onlyoi)go to 1200
      disz = -sign(perlen(3),dis(3,j))
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
c     if next nearest image is out of range, subsequent ones will also be
      if (rim.gt.rctsqn) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      fpn(nneigh) = fp(j)
      if(.not.threei)go to 1200
      disz = -disz
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
      if (rim.gt.rctsqn) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      fpn(nneigh) = fp(j)
      if(.not.fouri)go to 1200
      disz = -2.*sign(perlen(3),dis(3,j))
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
      if (rim.gt.rctsqn) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      fpn(nneigh) = fp(j)
1200  continue
c
c       now do diagonal (i=j) term for three or four images of self
c       both cases produce two images
c
      nneips = nneigh
      if(threei)then
c       first image of self
         nneips = nneips + 1
         disz = perlen(3)
         rim = disz**2
c     don't need to check against range here
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
c       second image of self
         nneips = nneips + 1
         disz = -disz
c        rim = disz**2    done above
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
      endif
      if(nneips.gt.neimax)then
          write(6,9012)nneigh,neimax
9012      format(' number of neighbors',i5,' exceeds array bound ',i5)
          stop
      endif
c
c  compute the current maximum number of neighbors seen
c
      nneimx = max0(nneimx,nneips)
      mxnnei = max0(mxnnei,nneimx)
c
c  compute the timing
c
c       call timeused(icpu,io,isys,imem)
c       gnetim = (1.e-6)*float(icpu) - timin
c       gnetmx = max(gnetmx,gnetim)
c       gnetmn = min(gnetmn,gnetim)
      return
c
c  end of nmeth = 1
c
c------------------------------
2000  continue
c
c  beginning of nmeth = 2
c    storage of neighbor indices
c
c newlst specifies whether the neighbor list is being updated
c
      if (newlst.ne.0) goto 2500
c
c  if here then use the neighbors found earlier
c
      jend = nnindx(i) - nnindx(i-1)
      do 2108 jtmp = 1,jend
      kx(jtmp) = nnlst(jtmp+nnindx(i-1))
      dis(1,jtmp) = rv(1,kx(jtmp))
      dis(2,jtmp) = rv(2,kx(jtmp))
      dis(3,jtmp) = rv(3,kx(jtmp))
2108  continue
      do 2110 jtmp = 1,jend
	call rdx3(rv(1,i),dis(1,jtmp),per12,perlen,cmono,
     1	dis(1,jtmp),r(jtmp))
2110  continue
      nneigh = 0
      do 2100 jtmp = 1,jend
c
c     determine which pairs are separated by less than rcut
c     and store the needed information about these pairs
c
c     if nearest periodic image is out of range, then all images
c     will be
c
      if (r(jtmp).gt.rctsqn) go to 2100
      j = kx(jtmp)
      nneigh = nneigh + 1
      rneigh(nneigh) = r(jtmp)
      jneigh(nneigh) = j
      do 2120 kcoord = 1,3
2120  dneigh(kcoord,nneigh) = dis(kcoord,jtmp)
      fpn(nneigh) = fp(j)
c     check periodic images in z direction
      if(onlyoi)go to 2100
      disz = -sign(perlen(3),dis(3,jtmp))
      rim = r(jtmp) + 2.*dis(3,jtmp)*disz + disz**2
c     if next nearest image is out of range, subsequent ones will also b
      if (rim.gt.rctsqn) go to 2100
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,jtmp)
      dneigh(2,nneigh) = dis(2,jtmp)
      dneigh(3,nneigh) = dis(3,jtmp) + disz
      fpn(nneigh) = fp(j)
      if(.not.threei)go to 2100
      disz = -disz
      rim = r(jtmp) + 2.*dis(3,jtmp)*disz + disz**2
      if (rim.gt.rctsqn) go to 2100
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,jtmp)
      dneigh(2,nneigh) = dis(2,jtmp)
      dneigh(3,nneigh) = dis(3,jtmp) + disz
      fpn(nneigh) = fp(j)
      if(.not.fouri)go to 2100
      disz = -2.*sign(perlen(3),dis(3,jtmp))
      rim = r(jtmp) + 2.*dis(3,jtmp)*disz + disz**2
      if (rim.gt.rctsqn) go to 2100
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,jtmp)
      dneigh(2,nneigh) = dis(2,jtmp)
      dneigh(3,nneigh) = dis(3,jtmp) + disz
      fpn(nneigh) = fp(j)
2100  continue
c
c       now do diagonal (i=j) term for three or four images of self
c       both cases produce two images
c
      nneips = nneigh
      if(threei)then
c       first image of self
         nneips = nneips + 1
         disz = perlen(3)
         rim = disz**2
c     don't need to check against range here
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
c       second image of self
         nneips = nneips + 1
         disz = -disz
c        rim = disz**2    done above
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
      endif
      if(nneips.gt.neimax)then
          write(6,9212)nneips,neimax
9212      format(' old neighbor list '/
     1      ' number of neighbors',i10,' exceeds array bound ',i5)
	write (6,*) 'at atom number ',i
	write (6,9999) (jneigh(kkk),rneigh(kkk),(dneigh(kk,kkk),kk=1,3)
     1       ,kkk=1,nneips)
9999	format(i6,4f10.5)
          stop
      endif
c
c  compute the maximum number of neighbors seen so far
c
      nneimx = max0(nneimx,nneips)
c
c  compute the timing
c
c       call timeused(icpu,io,isys,imem)
c       gnetim = (1.e-6)*float(icpu) - timin
c       gnetmx = max(gnetmx,gnetim)
c       gnetmn = min(gnetmn,gnetim)
      return
c
c  end of neighbor finding when using old neighbor list
c
c-------------------------------------------------------------
2500  continue
c
c  determine new neighbor list while getting the neighbors
c  nmax is the maximum number of atoms to look at
c  jmax is the maximum atom number to look at
c
	if(nmeth.gt.0) then
        if(full) then
        nmax=natoms
        else
        nmax=i-1
        jmax=i-1
        endif
	else
	call getnei(rv(1,i),mcell,natc,natcmax,icmax,icmay,icmaz,
     1    nei,neigh,neimx,cmono)
	nmax=nei
	 if(full) then
	 jmax=natoms
	 else
	 jmax=i-1
	 endif
	endif
      do 2600 n = 1,nmax
	if(nmeth.gt.0) then
	j=n
	else
	j=neigh(n)
	if(j.gt.jmax) go to 2600
	endif
c
c     compute the square of the distance to the closest periodic image
c
	call rdx3(rv(1,i),rv(1,j),per12,perlen,cmono,dis(1,j),r(j))
2600  continue
c
      nnindx(i) = nnindx(i-1)
      nneigh = 0
      do 2700 n = 1,nmax
        if(nmeth.gt.0) then
        j=n
        else
        j=neigh(n)
        if(j.gt.jmax) go to 2700
        endif

c
c  determine if these particles are within the storage distance
c
      if (r(j).gt.rctsqn.or.j.eq.i) goto 2700
c
c  store the index of the particle
c
      nnindx(i) = nnindx(i) + 1
      nnlst(nnindx(i)) = j
c
c     determine which pairs are separated by less than rcut
c     and store the needed information about these pairs
c
c     if nearest periodic image is out of range, then all images
c     will be
c
      if (r(j).gt.rctsqn) go to 2700
      nneigh = nneigh + 1
	if(nneigh.gt.neimax) then
	nneips=nneigh
	go to 7777
	endif
      rneigh(nneigh) = r(j)
      jneigh(nneigh) = j
      do 2750 kcoord = 1,3
2750  dneigh(kcoord,nneigh) = dis(kcoord,j)
      fpn(nneigh) = fp(j)
c     check periodic images in z direction
      if(onlyoi)go to 2700
      disz = -sign(perlen(3),dis(3,j))
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
c     if next nearest image is out of range, subsequent ones will also be
      if (rim.gt.rctsqn) go to 2700
      nneigh = nneigh + 1
        if(nneigh.gt.neimax) then
        nneips=nneigh
        go to 7777
        endif
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      fpn(nneigh) = fp(j)
      if(.not.threei)go to 2700
      disz = -disz
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
      if (rim.gt.rctsqn) go to 2700
      nneigh = nneigh + 1
        if(nneigh.gt.neimax) then
        nneips=nneigh
        go to 7777
        endif
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      fpn(nneigh) = fp(j)
      if(.not.fouri)go to 2700
      disz = -2.*sign(perlen(3),dis(3,j))
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
      if (rim.gt.rctsqn) go to 2700
      nneigh = nneigh + 1
        if(nneigh.gt.neimax) then
        nneips=nneigh
        go to 7777
        endif
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      fpn(nneigh) = fp(j)
2700  continue
c
c       now do diagonal (i=j) term for three or four images of self
c       both cases produce two images
c
      nneips = nneigh
      if(threei)then
c       first image of self
         nneips = nneips + 1
        if(nneips.gt.neimax) go to 7777
         disz = perlen(3)
         rim = disz**2
c     don't need to check against range here
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
c       second image of self
         nneips = nneips + 1
        if(nneips.gt.neimax) go to 7777
         disz = -disz
c        rim = disz**2    done above
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
      endif
7777     if(nneips.gt.neimax)then
          write(6,9222)nneips,neimax
9222      format(' new neighbor list '/
     1      ' number of neighbors',i10,' exceeds array bound ',i5)
        write (6,*) 'at atom number ',i
        write (6,9999) (jneigh(kkk),rneigh(kkk),(dneigh(kk,kkk),kk=1,3)
     1       ,kkk=1,nneips-1)
          stop
      endif
c
c  compute the maximum number of neighbors seen so far
c
      nneimx = max0(nneimx,nneips)
      mxnnei = max0(mxnnei,nnindx(i)-nnindx(i-1))
      mxlstu = max0(mxlstu,nnindx(i))
c
c  if this is the last call to gneigh, then set newlst=0 to indicate
c  that the current list can be used and also increment ngtlst
c
      if (i.eq.natoms) then
         newlst = 0
         ngtlst = ngtlst + 1
       end if
c
c  compute the timing
c
c       call timeused(icpu,io,isys,imem)
c       gnetim = (1.e-6)*float(icpu) - timin
c       gnetmx = max(gnetmx,gnetim)
c       gnetmn = min(gnetmn,gnetim)
      return
c
c  end of nmeth = 2
c
3000  continue
c
c  beginning of nmeth = 3
c    storage of neighbor indices with incremental update of
c    neighbor list
c
      jtmp = 1
      do 3105 jlst = 1,nnindx(i)
      if (nlist3(jlst,i).ge.i) goto 3106
      dis(1,jtmp) = rv(1,nlist3(jlst,i))
      dis(2,jtmp) = rv(2,nlist3(jlst,i))
      dis(3,jtmp) = rv(3,nlist3(jlst,i))
      kx(jtmp) = nlist3(jlst,i)
      jtmp = jtmp + 1
3105  continue
3106  continue
      jend = jtmp - 1
      do 3110 jtmp = 1,jend
	call rdx3(rv(1,i),dis(1,jtmp),per12,perlen,cmono,
     1	dis(1,jtmp),r(jtmp))
3110  continue
      nneigh = 0
      do 3100 jtmp = 1,jend
c
c     determine which pairs are separated by less than rcut
c     and store the needed information about these pairs
c
c     if nearest periodic image is out of range, then all images
c     will be
c
      if (r(jtmp).gt.rctsqn) go to 3100
      j = kx(jtmp)
      nneigh = nneigh + 1
      rneigh(nneigh) = r(jtmp)
      jneigh(nneigh) = j
      do 3120 kcoord = 1,3
3120  dneigh(kcoord,nneigh) = dis(kcoord,jtmp)
      fpn(nneigh) = fp(j)
c     check periodic images in z direction
      if(onlyoi)go to 3100
      disz = -sign(perlen(3),dis(3,jtmp))
      rim = r(jtmp) + 2.*dis(3,jtmp)*disz + disz**2
c     if next nearest image is out of range, subsequent ones will also be
      if (rim.gt.rctsqn) go to 3100
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,jtmp)
      dneigh(2,nneigh) = dis(2,jtmp)
      dneigh(3,nneigh) = dis(3,jtmp) + disz
      fpn(nneigh) = fp(j)
      if(.not.threei)go to 3100
      disz = -disz
      rim = r(jtmp) + 2.*dis(3,jtmp)*disz + disz**2
      if (rim.gt.rctsqn) go to 3100
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,jtmp)
      dneigh(2,nneigh) = dis(2,jtmp)
      dneigh(3,nneigh) = dis(3,jtmp) + disz
      fpn(nneigh) = fp(j)
      if(.not.fouri)go to 3100
      disz = -2.*sign(perlen(3),dis(3,jtmp))
      rim = r(jtmp) + 2.*dis(3,jtmp)*disz + disz**2
      if (rim.gt.rctsqn) go to 3100
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,jtmp)
      dneigh(2,nneigh) = dis(2,jtmp)
      dneigh(3,nneigh) = dis(3,jtmp) + disz
      fpn(nneigh) = fp(j)
3100  continue
c
c       now do diagonal (i=j) term for three or four images of self
c       both cases produce two images
c
      nneips = nneigh
      if(threei)then
c       first image of self
         nneips = nneips + 1
         disz = perlen(3)
         rim = disz**2
c     don't need to check against range here
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
c       second image of self
         nneips = nneips + 1
         disz = -disz
c        rim = disz**2    done above
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
         fpn(nneips) = fp(i)
      endif
      if(nneips.gt.neimax)then
          write(6,9312)nneips,neimax
9312      format(' number of neighbors',i5,' exceeds array bound ',i5)
          stop
      endif
c
c  compute the maximum number of neighbors seen so far
c
      nneimx = max0(nneimx,nneips)
c
c  compute the timing
c
c       call timeused(icpu,io,isys,imem)
c       gnetim = (1.e-6)*float(icpu) - timin
c       gnetmx = max(gnetmx,gnetim)
c       gnetmn = min(gnetmn,gnetim)
      return
c
c  end of neighbor finding when using nmeth=3
c
      end
c*********************************************************************
c
c  this subroutine computes the displacement of each particle since
c  the last update of the neighbor list.
c  if the maximum displacement is more than 1/2 of dradn, then newlst
c  is set to flag the creation of a new neighbor list by gneigh
c
      subroutine chkdis
      include 'param.inc'
      include 'utils.inc'
      include 'wrkspc.inc'
      common /nm3rold/ rold(3,natmax)
      dimension per12(3)
      dimension lstold(neimax)
	cmono=ctbeta*perlen(3)
      do 10 kcoord = 1,3
10    per12(kcoord) = 0.50*perlen(kcoord)

c	write(*,*) 'per12=',per12(1),per12(2),per12(3)
c	write(*,*) 'rctsqn=',rctsqn
c
c  branch to appropriate neighbor method
c
      goto (1000,2000,3000) abs(nmeth)
1000  return
c
c  nmeth = 2
c
2000  continue
c
c  treat the second call the same as the first in case defects have been
c  added
c
c  also if newlst has been set to 1, do not override even if the particles
c  have not moved
c
      if (newlst.eq.1) goto 2500
c
c  compare the new positions with the old ones
c
c     per12x = 0.5*perlen(1)
c     per12y = 0.5*perlen(2)
c     per12z = 0.5*perlen(3)
      do 2100 i = 1,natoms
	call rdx3(rold(1,i),rv(1,i),per12,perlen,cmono,
     1	dis(1,i),r(i))
2100  continue
c
c  determine the maximum displacement and compare with  dradn
c
      drmax1 = 0.0
      drmax2 = 0.0
      do 2200 i = 1,natoms
      tmp = min(drmax1,r(i))
      drmax1 = max(drmax1,r(i))
      drmax2 = max(drmax2,tmp)
2200  continue
      drmax = sqrt(drmax1) + sqrt(drmax2)
      if (drmax.gt.dradn) goto 2500
c
c  if here the old nighbor list can be used
c
      newlst = 0
      return
c-----------------
2500  continue
c
c  if here, a new neighbor list is needed so store the current coordinates
c
      newlst = 1
	if(nmeth.lt.0) then
	rcell=sqrt(rcutsq)+dradn
     	call getcell(rv,natoms,mcell,natc,natcmax,icmax,icmay,
     1    icmaz,rcell)
	endif
      do 2600 j = 1,3
      do 2600 i = 1,natoms
2600  rold(j,i) = rv(j,i)
      return
c
c  nmeth = 3
c
3000  continue
c
c  for the first two calls (the second is to allow for defect creation)
c  or if newlst has been set, create the neighbor lists from scratch
c
      if (newlst.eq.1) then
         call inmth3
         newlst = 0
         return
       endif
c
c  determine the displacement of all of the particles from there
c  last update position
c
c     per12x = 0.5*perlen(1)
c     per12y = 0.5*perlen(2)
c     per12z = 0.5*perlen(3)
      do 3100 i = 1,natoms
	call rdx3(rold(1,i),rv(1,i),per12,perlen,cmono,
     1	dis(1,i),r(i))
3100  continue
c
c  for each distance determine if the list for that atom
c  should be updated
c
      tst = (0.5*dradn)**2
      do 3200 i = 1,natoms
      if (r(i).le.tst) goto 3200
c
c  atom i's list needs to be updated
c
c  reset the old position
c
      rold(1,i) = rv(1,i)
      rold(2,i) = rv(2,i)
      rold(3,i) = rv(3,i)
c
c  find the new set of neighbors
c
      do 3300 j = 1,natoms
	call rdx3(rold(1,j),rv(1,i),per12,perlen,cmono,
     1	dis(1,j),p(j))
3300  continue
c
c store the old neighbor list for this atom
c
      do 3350 j = 1,nnindx(i)
3350  lstold(j) = nlist3(j,i)
      numold = nnindx(i)
c
c  create the new list for this atom
c
      nnindx(i) = 0
      do 3400 j = 1,natoms
      if (p(j).ge.rctsqn.or.j.eq.i) goto 3400
      nnindx(i) = nnindx(i) + 1
      nlist3(nnindx(i),i) = j
3400  continue
c
c  check that the number of neighbors has not exceeded neimax
c
      if (nnindx(i).gt.neimax) then
         write(6,9901) nnindx(i),i
9901     format(1x,'number of neighbors,',i5,'for atom,',i5)
         stop
       endif
c
c  now correct all the other lists
c
c  the algorithm uses the fact that both lists are in ascending order
c
      kold = 1
      knew = 1
3500  continue
      if (lstold(kold)-nlist3(knew,i)) 3510,3520,3530
3510  continue
c
c nlist3(knew,i) > lstold(kold)
c
      call nm3del(i,lstold(kold))
      kold = kold + 1
      if (kold.le.numold) goto 3500
      goto 3540
3520  continue
c
c nlist(knew,i) = lstold(kold)
c
      kold = kold + 1
      knew = knew + 1
      if (knew.le.nnindx(i).and.kold.le.numold) goto 3500
      goto 3540
3530  continue
c
c  lstold(kold) > nlist3(knew,i)
c
      call nm3add(i,nlist3(knew,i))
      knew = knew + 1
      if (knew.le.nnindx(i)) goto 3500
3540  continue
c
c  take care of the ends of the lists
c
      do 3550 kount = kold,numold
3550  call nm3del(i,lstold(kount))
      do 3560 kount = knew,nnindx(i)
3560  call nm3add(i,nlist3(kount,i))
3200  continue
      return
      end
c**********************************************************************
c
c  subroutine inmth3 initializes the neighbor lists for neighbor
c  method 3
c
      subroutine inmth3
      include 'param.inc'
      include 'utils.inc'
      include 'wrkspc.inc'
      common /nm3rold/ rold(3,natmax)
      dimension per12(3)
	cmono=ctbeta*perlen(3)
      do 10 kcoord = 1,3
10    per12(kcoord) = 0.50*perlen(kcoord)
c
c  initialize the rold array
c
      do 100 i=1,natoms
      rold(1,i) = rv(1,i)
      rold(2,i) = rv(2,i)
      rold(3,i) = rv(3,i)
100   continue
c
c  clear the nnindx array
c
      do 200 i = 1,natoms
200   nnindx(i) = 0
c
c  determine the periodicity data
c
c     per12x = 0.5*perlen(1)
c     per12y = 0.5*perlen(2)
c     per12z = 0.5*perlen(3)
c
c  loop over all the atom pairs and store the close neighbors
c
      do 1000 j2 = 1,natoms
      do 1100 j1 = 1,j2-1
	call rdx3(rold(1,j1),rold(1,j2),per12,perlen,cmono,
     1	dis(1,j1),r(j1))
1100  continue
c
c  store the indices for the neighbor pairs
c
      do 1200 j1 = 1,j2-1
      if (r(j1).ge.rctsqn) goto 1200
      nnindx(j1) = nnindx(j1) + 1
      nlist3(nnindx(j1),j1) = j2
      nnindx(j2) = nnindx(j2) + 1
      nlist3(nnindx(j2),j2) = j1
1200  continue
1000  continue
c
c  check that the number of neighbors has not exceeded neimax
c
      do 2000 j1=1,natoms
      if (nnindx(j1).gt.neimax) then
         write(6,9901) nnindx(j1),j1
9901     format(1x,'number of neighbors,',i5,'for atom,',i5)
         stop
      endif
2000  continue
      return
      end
c**********************************************************************
	subroutine finalpr(ekin1,pot1,etot1,eperf,volint)
      include 'param.inc'
      include 'utils.inc'
      if(abs(ipatoms).ge.2) call patoms
      call calce(ekin2,pot2,etot2,temper,press,1,tend)
      volume = perlen(1)*perlen(2)*perlen(3)/natoms
	write(6,9301) natoms,' total number of atoms'
9301    format(i20,A)
9300    format(g20.12,' final ',A)
        write(6,9300) temper,'temperature'
        write(6,9300) press,'pressure'
        write(6,9300) volume,'volume'
        write(6,9300) ekin2,'kinetic energy'
        write(6,9300) pot2,'potential energy'
        write(6,9300) etot2,'total energy'
        write(6,9300) gg,'g dot g'
        write(6,9300) perlen(1),'x period'
        write(6,9300) perlen(2),'y period'
        write(6,9300) perlen(3),'z period'
c
c  for a minimization, print out a warning if the minimization
c  terminated for any reason other than convergence to tol
c
      if (inte.eq.-1.and.mnterm.eq.1) then
         write(6,9232)
9232     format(1x,'warning************************************',/,
     1          1x,'minimization terminated due to nfmax',/,
     2          1x,'warning************************************')
       endif
      if (inte.eq.-1.and.mnterm.eq.-1) then
         write(6,9234)
9234     format(1x,'warning************************************',/,
     1          1x,'minimization terminated due to time limit',/,
     2          1x,'warning************************************')
       endif
c
c  determine the total work done by the external forces
c
      work = 0.0
      if (ifxal2.ne.0) then
         do 90 i = 1,natoms
90       work = work + workfx(i)
         if (work.ne.0.0) write(6,9242) work
9242     format(1x,'work done by external forces:',g15.6)
       endif
c
c  determine the changes in the energies
c
      dekin=ekin2-ekin1
      dpot=pot2-pot1
      detot=etot2-etot1
      erel=pot2-eperf
      write(6,9245)dekin,dpot,detot
 9245 format('  changes in energies: kinetic, potential, and total ',
     1/1x,3e17.7)
      write(6,9250)detops
 9250 format('  largest change in total e reported by output: ',e17.7)
      write(6,9260)erel
 9260 format('  energy of final lattice with defects',
     1' relative to initial lattice without defects =',e17.7)
      write(6,9270)((stresst(i,j),i=1,3),j=1,3)
 9270 format(/'  stress tensor'/' sx',3g15.5/' sy',3g15.5/' sz',3g15.5)
      dvperc = 100.0*(volume-volint)/(volint)
      write(6,9290)dvperc
 9290 format('   percent volume change:',g12.5)
	return
	end
c**********************************************************************
c
c   calculate the averages
	subroutine average(init,t)
      include 'param.inc'
      include 'utils.inc'
	dimension cper(3),cper2(3),cstrs(3,3),cstrs2(3,3),ctmprg(0:6),
     1     cvcm(0:6,3),cforce(0:6,3)
	dimension ctmprg2(0:6),cvcm2(0:6,3),cforce2(0:6,3)
	save cper,cper2,cstrs,cstrs2,ctemp,ctemp2,cpres,cpres2
	save cke,cke2,cpe,cpe2,cvol,cvol2,ctmprg,cvcm,cforce
	save ctmprg2,cvcm2,cforce2
c   initialize run averages
	if(init.eq.-1) then
         noutput=0
         ctemp=0
         ctemp2=0
         cpres=0
         cpres2=0
         cke=0
         cke2=0
         cpe=0
         cpe2=0
         cvol=0
         cvol2=0
	 do 501 ireg=0,6
         do 16 j=1,3
         cvcm2(ireg,j)=0.0
         cvcm(ireg,j)=0.0
         cforce2(ireg,j)=0.0
16       cforce(ireg,j)=0.0
	 ctmprg2(ireg)=0
501	 ctmprg(ireg)=0
         do 500 i=1,3
         cper(i)=0
         cper2(i)=0
         do 500 j=1,3
         cstrs(i,j)=0
500      cstrs2(i,j)=0
	elseif (init.eq.0) then
c   initialize averages
	 nave = 0
         avpe = 0.0
         avvol=0
         avke = 0.0
         avpres = 0.0
         avtemp = 0.0
         do 10 i = 1,3
         avper(i)=0.
         do 10 j = 1,3
         do 1011 k=1,natoms
1011     avslocal(i,j,k)=0.
10       avstrs(i,j) = 0.0
         do 11 ireg = 0,6
	 do 12 j=1,3
	 avvcm(ireg,j)=0.0
12	 avforce(ireg,j)=0.0
11       avtmprg(ireg) = 0.0
c        do 891 j=1,3
c        do 891 i=1,3
c        do 891 k=1,natoms
c891      avslocal(i,j,k)=0.
	elseif (init.eq.1) then
c   increment averages
         call calce(ekin,pot,etot,temp,press,0,0.0d0)
         nave = nave + 1
         avpres = avpres + press
         avtemp = avtemp + temp
         avvol = avvol + perlen(1)*perlen(2)*perlen(3)
         avke = avke + ekin
         avpe = avpe + pot
         do 100 i = 1,3
         avper(i) = avper(i) + perlen(i)
         do 100 j = 1,3
100      avstrs(i,j) = avstrs(i,j) + stresst(i,j)
         do 101 ireg = 0,nregs
	 if(natrg(ireg).gt.0) then
          do 13 j=1,3
          avvcm(ireg,j)= avvcm(ireg,j) + actvcm(ireg,j)
13        avforce(ireg,j)= avforce(ireg,j) +actforce(ireg,j)
	  avtmprg(ireg) = avtmprg(ireg) + acttmp(ireg)
	 endif
101	 continue
	elseif(init.eq.2) then
         nskipd = nskipd + 1
c   calculate interim final values
      if (nave.gt.0) then
	avtemp = avtemp/nave
        avpres = avpres/nave
        avke = avke/nave
        avpe = avpe/nave
        avvol = avvol/nave
        avden = natoms/avvol
        do 300 i = 1,3
        avper(i)=avper(i)/nave
        do 300 j = 1,3
300     avstrs(i,j) = avstrs(i,j)/nave
        do 102 ireg = 0,nregs
	if(natrg(ireg).gt.0) then
         do 14 j=1,3
         avvcm(ireg,j)=avvcm(ireg,j)/nave
14       avforce(ireg,j)=avforce(ireg,j)/nave
	endif
102     avtmprg(ireg) = avtmprg(ireg)/nave
        do 880 j=1,3
        do 880 i=1,3
        do 880 k=1,natoms
880     avslocal(i,j,k)=avslocal(i,j,k)/nave
        etot=avke+avpe
        write(6,9001)t,avke,avpe,etot
 9001   format(' output: time,kin,pot,tot',f9.3,3e13.5)
        write(6,9002)avtemp,avpres
 9002   format(9x,'temperature:',g12.5,5x,'pressure:',g12.5)
        do 1235 ireg=0,nregs
        if (ifxtmp.ne.0.and.natrg(ireg).gt.0) then
            write(6,9004) ireg,natrg(ireg),avtmprg(ireg)
9004        format('   region:',i5,'  number of atoms:',i10,
     .       '  temperature:',g12.5)
            if(icm) then
             write(6,9044) ireg,(avvcm(ireg,j),j=1,3)
9044         format('   region:',i5,' vcm:', 3g14.5)
             write(6,9045) ireg,(avforce(ireg,j),j=1,3)
9045         format('   region:',i5,' force:', 3g14.5)
            endif
         endif
1235     continue
         if (ibdtyp.eq.2) then
            write(6,9003)(avper(i),i=1,3)
9003        format(9x,'perlen:',3g17.8)
          end if
c   increment the averages for the entire run
	if(nskipd.gt.nequil) then
	 noutput=noutput+1
         ctemp=ctemp+avtemp
         ctemp2=ctemp2+avtemp**2
         cpres=cpres+avpres
         cpres2=cpres2+avpres**2
         cke=cke+avke
         cke2=cke2+avke**2
         cpe=cpe+avpe
         cpe2=cpe2+avpe**2
         cvol=cvol+avvol
         cvol2=cvol2+avvol**2
	 do 201 ireg=0,nregs
	 if(natrg(ireg).gt.0) then
          do 17 j=1,3
          cvcm(ireg,j)= cvcm(ireg,j) + avvcm(ireg,j)
          cvcm2(ireg,j)= cvcm2(ireg,j) + avvcm(ireg,j)**2
          cforce2(ireg,j)= cforce2(ireg,j) +avforce(ireg,j)**2
17        cforce(ireg,j)= cforce(ireg,j) +avforce(ireg,j)
	  ctmprg2(ireg)=ctmprg2(ireg)+avtmprg(ireg)**2
	  ctmprg(ireg)=ctmprg(ireg)+avtmprg(ireg)
	 endif
201	 continue
	 do 200 i=1,3
	 cper(i)=cper(i)+avper(i)
	 cper2(i)=cper2(i)+avper(i)**2
	 do 200 j=1,3
	 cstrs(i,j)=cstrs(i,j)+avstrs(i,j)
200	 cstrs2(i,j)=cstrs2(i,j)+avstrs(i,j)**2
	endif
      endif
        elseif(init.eq.3) then
c   calculate final average values
      if (noutput.gt.0) then
        avtemp = ctemp/noutput
        avtemp2 = sqrt(ctemp2/noutput-avtemp**2)
        avpres = cpres/noutput
        avpres2 = sqrt(cpres2/noutput-avpres**2)
        avke = cke/noutput
        avke2 = sqrt(cke2/noutput-avke**2)
        avpe = cpe/noutput
        avpe2 = sqrt(cpe2/noutput-avpe**2)
        avvol = cvol/noutput
        avvol2 = sqrt(abs(cvol2/noutput-avvol**2))
        avden = natoms/avvol
        do 400 i = 1,3
        avper(i)=cper(i)/noutput
        avper2(i)=sqrt(abs(cper2(i)/noutput-avper(i)**2))
        do 400 j = 1,3
        avstrs(i,j) = cstrs(i,j)/noutput
400     avstrs2(i,j) = sqrt(cstrs2(i,j)/noutput-avstrs(i,j)**2)
	do 401 ireg=0,nregs
	avtmprg(ireg)=ctmprg(ireg)/noutput
        ctmprg2(ireg)=sqrt(abs(ctmprg2(ireg)/noutput-avtmprg(ireg)**2))
        do 15 j=1,3
        avvcm(ireg,j)=cvcm(ireg,j)/noutput
        cvcm2(ireg,j)=sqrt(abs(cvcm2(ireg,j)/noutput-avvcm(ireg,j)**2))
        avforce(ireg,j)=cforce(ireg,j)/noutput
        cforce2(ireg,j)=sqrt(abs(cforce2(ireg,j)/noutput
     1       -avforce(ireg,j)**2))
15      continue
401	continue
9300    format(g20.7,' average ',A)
	write(6,9300) avtemp,'temperature'
	write(6,9300) avpres,'pressure'
	avvol=avvol/natoms
	write(6,9300) avvol,'volume'
	write(6,9300) avden,'density'
	write(6,9300) avke,'kinetic energy'
	write(6,9300) avpe,'potential energy'
	write(6,9300) avper(1),'x period'
	write(6,9300) avper(2),'y period'
	write(6,9300) avper(3),'z period'
        write(6,9301) ,' x ',(avstrs(1,j),j=1,3)
        write(6,9301) ,' y ',(avstrs(2,j),j=1,3)
        write(6,9301) ,' z ',(avstrs(3,j),j=1,3)
9301    format('average stress tensor',A,3g15.5)
9305    format(g20.7,' RMS ',A)
        write(6,9305) avtemp2,'temperature'
        write(6,9305) avpres2,'pressure'
	avvol2=avvol2/natoms
        write(6,9305) avvol2,'volume'
        write(6,9305) avke2,'kinetic energy'
        write(6,9305) avpe2,'potential energy'
        write(6,9305) avper2(1),'x period'
        write(6,9305) avper2(2),'y period'
        write(6,9305) avper2(3),'z period'
        write(6,9303) ,' x ',(avstrs2(1,j),j=1,3)
        write(6,9303) ,' y ',(avstrs2(2,j),j=1,3)
        write(6,9303) ,' z ',(avstrs2(3,j),j=1,3)
9303    format('RMS stress tensor',A,3g15.5)
        do 1234 ireg=0,nregs
        if (ifxtmp.ne.0.and.natrg(ireg).gt.0) then
         write(6,9304) ireg,avtmprg(ireg),ctmprg2(ireg)
9304     format(' temperature of region:',i5,' average:',g14.7,
     1     ' RMS:',g14.7)
         if(icm) then
           write(6,9046) ireg,(avvcm(ireg,j),j=1,3)
9046       format(' region:',i5,' average vcm:', 3g14.5)
           write(6,9047) ireg,(avforce(ireg,j),j=1,3)
9047       format(' region:',i5,' average force:', 3g14.5)
           write(6,9048) ireg,(cvcm2(ireg,j),j=1,3)
9048       format(' region:',i5,' RMS vcm:', 3g14.5)
           write(6,9049) ireg,(cforce2(ireg,j),j=1,3)
9049       format(' region:',i5,' RMS force:', 3g14.5)
         endif
        endif
1234    continue
	endif
	endif
	return
	end
****************************************************************
	subroutine getbnd
      include 'param.inc'
      include 'utils.inc'
      namelist /bndcard/ ibdtyp,dpress,bndmas,bnddrg,dstress,idynper
      data ibdtyp/1/,dpress/0.0/,bndmas/3*100.0/,bnddrg/3*-1.0/
      data idynper/3*1/,dstress/3*0./
        data iper /0/
c
c  read in the type of boundaries to be used
c
      read(5,bndcard)
c
c  determine the boundary drag if it is not specified by the user
c  this determination is for critical damping assuming a
c  bulk modulus of 1mbar
c
      if (ibdtyp.eq.2.and.bnddrg(1).lt.0.0) then
         bnddrg(1) = sqrt(7.5*perlen(1)/bndmas(1))
         bnddrg(2) = sqrt(7.5*perlen(2)/bndmas(2))
         bnddrg(3) = sqrt(7.5*perlen(3)/bndmas(3))
       end if
        do 777 k=1,3
        if(idynper(k).eq.2) bnddrg(k)=0
777     continue
      if (ibdtyp.eq.2) then
         trace = (dstress(1)+dstress(2)+dstress(3))/3.
         dpress = dpress + trace
         dstress(1) = dstress(1) - trace
         dstress(2) = dstress(2) - trace
         dstress(3) = dstress(3) - trace
       endif
      write(6,9210)
9210  format(/,/' ******  boundary conditions')
      if(ibdtyp.eq.1)write(6,9215)
9215  format(/,'  constant volume (fixed periodic vectors) ')
      if(ibdtyp.eq.2)then
        write(6,9220)
9220    format(/,'  dynamic periodic lengths (fixed directions) ')
        write(6,9225)dpress,dstress,bndmas,bnddrg
9225    format('  desired pressure:',g12.5,/,
     1         '  desired stresses:',3g14.5,/,
     2'  boundary mass:',3g14.5,/,'  boundary drag:',3g14.5,/)
        write(6,9227) idynper
9227    format('  boundaries allowed to move:',3i2)
        if( idynper(1).lt.0.and.idynper(2).lt.0
     1     .and.idynper(3).lt.0) then
        iper=4
        elseif(idynper(1).lt.0.and.idynper(2).lt.0) then
        iper=3
        elseif(idynper(1).lt.0.and.idynper(3).lt.0) then
        iper=2
        elseif(idynper(3).lt.0.and.idynper(2).lt.0) then
        iper=1
        else
        iper=0
        endif
      endif
c  convert pressure and stresses from bar to internal units
      dpress = dpress/1.602e+6
      dstress(1) = dstress(1)/1.602e+06
      dstress(2) = dstress(2)/1.602e+06
      dstress(3) = dstress(3)/1.602e+06
c
c  multiply dstress by the appropriate combinations of perlens
c  for use elsewhere
c
      dstress(1) = dstress(1)*perlen(2)*perlen(3)/perlen(1)
      dstress(2) = dstress(2)*perlen(1)*perlen(3)/perlen(2)
      dstress(3) = dstress(3)*perlen(1)*perlen(2)/perlen(3)
	return
	end
c**************************************************************************
      subroutine verlet(t,n,otime)
      
c  this routine implements the verlet integration scheme
c  (Ramon Ravelo 2/95
c  
      include 'param.inc'
      include 'utils.inc'
      dimension u1(3,natmax),u2(3,natmax)
c     
c  sets the acceleration and velocities at time t
c  
        print *,'in verlet'
      call scale_sub
      
      do l=1,3
         do i=1,n
            u1(l,i)=rv(l+3,i)
         enddo
      enddo
      
      call accel(t)
      
c
c  integrates the equations of motion
c  
c  first initialize the timing value
c     
      dt1=dt
      tinit = t
      tfinal = tinit + nsteps*dt
      tequil=tinit+eqtim
      delt=(tfinal-tequil)
      timstr = seconds()
      swval=0.0d0
1000  continue
      t = t + dt1
      
      if (iswitch.eq.0) then
         swval=0.0d0
      elseif(iswitch.eq.-1) then
         swval=1.0d0
      else
         if (t.le.tequil) then
            swval=swi
         else
            swval=swi+(t-tequil)*(swf-swi)/delt
         endif
      endif

c  Calculates velocities, positions and forces at t=t+dt1/2

      do l=1,3
         do i=1,n
            y(l+3,i)=u1(l,i) + 0.5d0*dt1*acc(l,i)
            y(l,i)=y(l,i) + dt1*y(l+3,i)
            u2(l,i)=y(l+3,i)
            y(l+3,i)=2*u2(l,i)-u1(l,i)
         enddo
      enddo

      call accel(t)

c     Updates the velocities to t=t+dt1


      do l=1,3
         do i=1,n
            y(l+3,i)=u2(l,i)+0.5d0*dt1*acc(l,i)
            u1(l,i)=y(l+3,i)
         enddo
      enddo
c
c  output the results if appropriate
c    
      if (abs(mod(t-tinit+0.5*dt1,otime)-0.5*dt1).lt.(0.5*dt1)) then
         call froutput(t,iaccur)
         tend = t
c
c     determine if there is enough time left to go to next output
c
         timend = seconds()
         timned = 1.2*timend
         call trmain(timlft)
         if (timlft.le.timned) then
            write(6,9201)
 9201       format(1x,'job terminated early due to time limit')
            return
         end if
         timstr = timend
      end if
c
c  see if done with the integration
c
      if (t.lt.tfinal) goto 1000
      return
      end
c************************************************************************
