

      implicit real*8 (a-h,o-z)
******************************************************************
*  Develop interpolation coefficients for interatomic potentials *
******************************************************************

      parameter(ngrid = 6000,icomp=2) 
      parameter(nspmax	  = 4)
      real rho(icomp,ngrid),rhop(icomp,ngrid)
      real phi(icomp,ngrid),phip(icomp,ngrid)
      real frho(icomp,ngrid),frhop(icomp,ngrid)
      real phix(ngrid),phixp(ngrid)
      character*4 element(nspmax)

      common /chemical/ cutoff2,xmass(nspmax),np(nspmax,nspmax),
     &                  numat,npairs,npr,element 

      call INTERPOL

      cutoff = sqrt(cutoff2)

      print*,'CUTOFF =',cutoff
      dr = cutoff/ngrid
      drar = dr
      drhoar = dr

      open(8,file='pal.output',status='unknown')
      open(9,file='palti.output',status='unknown')
      open(10,file='pti.output',status='unknown')
      open(11,file='for022.output',status='unknown')
      open(12,file='fti.output',status='unknown')
      open(14,file='fral.output',status='unknown')
      open(15,file='frti.output',status='unknown')

      r = 0.0D0
      do i=1,ngrid
         r = r + dr
         phi(1,i) = SPL(1,r,0)
         phix(i) = SPL(2,r,0)
         phi(2,i) = SPL(3,r,0)
         rho(1,i) = SPL(4,r,0)
         rho(2,i) = SPL(5,r,0)
         frho(1,i) = SPL(6,r,0)
         frho(2,i) = SPL(7,r,0)
         phip(1,i) = SPL(1,r,1)
         phixp(i) = SPL(2,r,1)
         phip(2,i) = SPL(3,r,1)
         rhop(1,i) = SPL(4,r,1)
         rhop(2,i) = SPL(5,r,1)
         frhop(1,i) = SPL(6,r,1) 
         frhop(2,i) = SPL(7,r,1) 

         write(8,*)r,phi(1,i),phip(1,i)
         write(9,*)r,phix(i),phixp(i)
         write(10,*)r,phi(2,i),phip(2,i)
         write(11,*)r,rho(1,i),rhop(1,i)
         write(12,*)r,rho(2,i),rhop(2,i)
         write(14,*)r,frho(1,i),frhop(1,i)
         write(15,*)r,frho(2,i),frhop(2,i)

      end do 

      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(14)
      close(15) 

      open(16,file='eamdata',status='unknown')

      write(16,*)'TiAl from Diana Farkas Al (1) Ti (2)'
      write(16,*)13, 26.98, 4.05
      write(16,*)drar,drhoar,cutoff

* Write rho for components 1 and 2

      do i=1,ngrid
         write(16,*)(rho(icom,i),rhop(icom,i),icom=1,2)
      end do

* Write phi for components 1 and 2 (phi11 and phi22)   

      do i=1,ngrid
         write(16,*)(phi(icom,i),phip(icom,i),icom=1,2)
      end do

* Write cross terms phi12

       do i=1,ngrid
           write(16,*)phix(i), phixp(i)
       end do 


* Write F(rho) for components 1 and 2 

      do i=1,ngrid
         write(16,*)(frho(icom,i),frhop(icom,i),icom=1,2)
      end do

      close(16)

      stop
      end

      SUBROUTINE INTERPOL
*-----------------------------------------
* Yuri Mishin, Feb. 1996; revised in June 1998
*===================================================================
* This subroutine reads in data from potential files and develops
* cubic spline coefficients that are subsequently used by function
* SPL to retrieve interpolated values of interatomic potentials and
* their derivatives if requested. The input data and the spline 
* coefficients are stored in arrays a, b, c, and d in common 
* block /pot/. 
*
* Subroutines used: SPLINE
*===================================================================

      implicit real*8 (a-h,o-z)
      parameter (nspmax   = 4,    ! max number of chem. species
     &           maxpnt   = 3000, ! max number of tabulated values
     &           maxf     = nspmax*(nspmax+1)/2 + 2*nspmax  
     &                            ! max number of potential files
     &) 
      character*4 element(nspmax)
      character*72 filename(maxf)
      common /chemical/ cutoff2,xmass(nspmax),np(nspmax,nspmax),
     &                  numat,npairs,npr,element 
      common/pot/a(maxf,maxpnt),b(maxf,maxpnt),c(maxf,maxpnt),
     &           d(maxf,maxpnt),x0(maxf),xn(maxf),
     &           step(maxf),n(maxf)
     &
      dimension x(maxpnt),y(maxpnt),w1(maxpnt),w2(maxpnt),w3(maxpnt)

* ==== Input the number of chemical species ==========

      open (1, file = 'dianasfiles', status= 'old')
      read (1,*) 
      read(1,*) numat

      if ((numat.lt.1).or.(numat.gt.4)) then
         write(6,100)
         stop
      endif

***********************************************************
*         Input names of potential files.                 *
* The first numat*(numat+1)/2 files represent pair        *
* interaction (e.g. for numat = 3, p11, p12, p13, p22,    *
* p23, p33). These are followed by numat electron-density *
* files and numat embedding-energy files.                 *
***********************************************************

      npairs = numat*(numat+1)/2
      npr=npairs + numat
      nfiles = npr  + numat

      write(6,103)
      do i=1,numat
         read(1,*) element(i), xmass(i)
         write(6, 104) i, element(i)
      enddo

      write(6,101)
      do i=1,nfiles
         read (1,*) filename(i)
         write(6,*) filename(i)
      enddo

*--- Array np(i,j) stores numberss of pair interaction files 
*--- for species i and j

      nf = 0
      do i=1,numat
         do j=i,numat
            nf = nf + 1
            np(i,j) = nf
            np(j,i) = nf
         enddo
      enddo

********************************************************
* Read in data from the potential files and calculate  *
* the spline coefficients                              *
********************************************************

      do i=1, nfiles
         nunit=10+i
         open (nunit,file=filename(i),status='old')    
         read(nunit,*) n(i), x0(i), xn(i), step(i), idumm1, idumm2
         read(nunit,*) (y(j), j=1,n(i))
         step(i)=1.0D0/step(i)
         do j=1,n(i)
            x(j)=x0(i)+step(i)*j
         enddo
         call SPLINE(n(i),x,y,w1,w2,w3)
         do j=1,n(i)
            a(i,j)=y(j)
            b(i,j)=w1(j)
            c(i,j)=w2(j)
            d(i,j)=w3(j)
         enddo
         close(nunit)
      enddo

*****************************************************
*Calculate the cut-off radius of atomic interaction *
*****************************************************
 
      cutoff=0.D0  
      do i=1, npr
         if (xn(i).gt.cutoff) cutoff=xn(i)
      enddo
      cutoff2=cutoff**2
      write(6, 102) cutoff

 100  format(//' WRONG NUMBER OF CHEMICAL SPECIES !!!'//)
 101  format(/' THE FOLLOWING POTENTIAL FILES ARE USED:'/)
 102  format(/' CUT-OFF RADIUS OF ATOMIC INTERACTION:  ',f8.4)
 103  format(/' CHEMICAL ELEMENTS:','   TYPE   ELEMENT'/)
 104  format(20x,i5,6x,a4)

      RETURN
      END


      SUBROUTINE SPLINE (n, x, y, b, c, d)
****************************************************************
*               Spline interpolation subroutine from 
*          "Computer Methods for Mathematical Computations"
*            by G.E. Forsythe, M.A. Malcolm and C.B. Moler
*             (Englewood Cliffs, NJ: Prentice-Hall, 1977)
****************************************************************

      implicit real*8 (a-h,o-z)
      dimension x(*), y(*), b(*), c(*), d(*)

* The coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
* for a cubic interpolating spline

* s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

* for  x(i) .le. x .le. x(i+1)

* input..

*   n = the number of data points or knots (n.ge.2)
*   x = the abscissas of the knots in strictly increasing order
*   y = the ordinates of the knots

* output..

*   b, c, d  = arrays of spline coefficients as defined above.

* using  p  to denote differentiation,

*   y(i) = s(x(i))
*   b(i) = sp(x(i))
*   c(i) = spp(x(i))/2
*   d(i) = sppp(x(i))/6  (derivative from the right)

* the accompanying function subprogram  seval  can be used
* to evaluate the spline.

      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50

* set up tridiagonal system

* b = diagonal, d = offdiagonal, c = right hand side.

      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue

* end conditions.  third derivatives at  x(1)  and  x(n)
* obtained from divided differences

      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.D0
      c(n) = 0.D0
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))

* forward elimination

   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue

* back substitution

      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue

* c(i) is now the sigma(i) of the text

* compute polynomial coefficients

      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3*c(i)
   40 continue
      c(n) = 3*c(n)
      d(n) = d(n-1)
      return

   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.D0
      d(1) = 0.D0
      b(2) = b(1)
      c(2) = 0.D0
      d(2) = 0.D0
      RETURN
      END

      FUNCTION SPL(k,u,m)
*-----------------------------------------
*  Yuri Mishin, Feb. 1996; revised in June 1998
*===================================================================
* Calculates the cubic spline function 

*   spl = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3

* and its 1st, 2nd and 3rd derivatives for an arbitrary argument u,
* where  x(i) .lt. u .lt. x(i+1).

* If  u .lt. x(1) then  i = 1  is used.
* If  u .ge. x(n) then  i = n  is used for the embedding function
*                       while spl = 0 is set for pair interaction
*                       and electronic density

* x and y are arrays of abscissas and ordinates, b, c and d are 
* arrays of spline coefficients computed by subroutine SPLINE. 
* Input parameter m indicates that the m-th derivative (1<= m <=3) 
* of spl(u) or only spl(u) (m = 0) must be returned.
* Up to maxf data sets, enumerated by index k, can be used in this
* function.
* =================================================================

      implicit real*8 (a-h,o-z)
      parameter (nspmax   = 4,    ! max number of chem. species
     &           maxpnt   = 3000, ! max number of tabulated values
     &           maxf     = nspmax*(nspmax+1)/2 + 2*nspmax  
     &                            ! max number of potential files
     &) 
      character*4 element(nspmax)
      common /chemical/ cutoff2,xmass(nspmax),np(nspmax,nspmax),
     &                  numat,npairs,npr,element 
      common/pot/a(maxf,maxpnt),b(maxf,maxpnt),c(maxf,maxpnt),
     &           d(maxf,maxpnt),x0(maxf),xn(maxf),
     &           step(maxf),n(maxf)
 
      if(u.gt.xn(k)) then
         if(k.le.npr) then 
            spl = 0.D0
            return
         endif
      endif

*------- Find i such that x(i) <= u < x(i+1) 

      ux=u-x0(k)
      s = step(k)
      if (ux.lt.s) then
         i=1
      else 
         i=INT(ux/s)
         if (i.gt.n(k)) i=n(k)
      endif

*------------ calculate spl(u) or its derivatives

      dx = ux - i*s

      if (m.eq.0) then
         spl = a(k,i) + dx*(b(k,i) + dx*(c(k,i) + dx*d(k,i)))
      else 
         if(m.eq.1) then
            spl = b(k,i) + dx*(2*c(k,i) + 3*dx*d(k,i))
         else
            if (m.eq.2) then
               spl = 2*c(k,i) + 6*dx*d(k,i)
            else
               if (m.eq.3) then
                  spl = 6*d(k,i)
               else
                  write(6,100) m
                  stop
               endif
            endif
         endif
      endif
 100  format(//' WRONG INPUT INDEX m = ', i3,
     &'  IN FUNCTION SPL !!!'//)

      RETURN
      END
