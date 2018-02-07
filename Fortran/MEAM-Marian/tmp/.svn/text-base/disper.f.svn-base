	program disper
	
*** Calculates the Phonon Dispersion Curves
	
	include "moldy.h"
	
	double precision eps, reps, rmass
	parameter (rmass = 283.00*1.0365e-28)
	parameter (eps = 0.01)
	parameter (reps = eps*alatt)
	double precision dm(3,3*nm)
	double precision sdm(3,3)
	double precision diag(3), off(3)
	
	call input

C	call forces
C	write(*,*) 'pe=',pe,'petrip=',petrip
C
C	x(1,1)  =x(1,1)+eps*b0
C	x(2,13) =x(2,13)+eps*b0
C	x(3,3)  =x(3,3)+eps*b0
C	call forces
C	do kl = 1, 13
C	   write(*,*) 'f(',kl,') = (',fx(kl),fy(kl),fz(kl)
C	enddo
C
C	stop


*** read initial and potential data
	
	do j = 1, 3
	   
*** Atom 515 is the one located at ( 0.0 0.0 0.0 )
*** only one atom's position in one dimesion changes at a time
*** conversion to unit-box units is performed in routine 'forces'
	   
	   x(j,43) = xi(j,43) + eps
	   call forces 
	   fplus = f
	   
	   x(j,43) = xi(j,43) - eps
	   call forces
	   fmin = f
	   
*** Forces come out in eV/Angstrom
	   
	   x = xi
*** we reassign the relaxed values to x so that in the next
**** iteration only one atom's position in one dimesion changes
	   
	   do ii = 1, nm
	      
	      do jj = 1, 3
		 
		 kk = 3*(ii-1) + jj
		 
		 dm(j,kk) = (fmin(jj,ii)-fplus(jj,ii))
     > 		      /(2.0*reps*dsqrt( rmass*rmass ))
		 
*** dm(3,3*nm) are the three x -- y -- z rows of the dynamical matrix correpsonding to atom 515
		 
	      enddo
	   enddo
	   print *, ' atom at the origin '
	   
	enddo
	
	do n = 1, nm
	   do m = 1, 3
	      nindex = 3*(n-1) + m
	      write (45,*) (dm(j,nindex), j = 1, 3)
	      
*** We dump all the row elements as columns in fort.45
	   enddo
	enddo
	
	
*** We now create the local dynamical matrix as a function of k
*	k-branches in bcc metals are:
*      	   k = [100], 0 < k < 2*pi/a
*      	   k = [111], 0 < k < 2*pi/a
*      	   k = [110], 0 < k < pi/a
*** a is the lattice parameter
	
	do k = 1, 20
	   
*** We take kmax pointos in each k-branch	
	   delta_k = 2.0* pi / 20.0 ! for [100] and [111] wavevectors
c	delta_k = pi / 10.0 ! for [110] wavevectors
	   
	   wkx = k*delta_k 
	   wky = k*delta_k 
c	wkz = k*delta_k 
c	wky = 0d0 
	   wkz = 0d0 
	   
	   do n = 1, 3
	      do m = 1, 3
		 
		 do j = 1, nm	
		    
		    prodk = wkx*xi(1,j) + wky*xi(2,j) + wkz*xi(3,j)    
		    sdm(n,m) = sdm(n,m) + ( dm(n,3*(j-1)+m) )*cos(prodk)
		    
		 enddo
		 
	      enddo
	   enddo
	   
*** We have generated a 3x3 matrix
*** Symmetrization
	   
	   do l = 1, 2
	      do ll = l+1, 3
		 sdm(l,ll) = ( sdm(l,ll) + sdm(ll,l) )/2d0
		 sdm(ll,l) = sdm(l,ll)
	      enddo
	   enddo
	   
	   print *, 'symmetrization performed, now entering TRED2...'
	   
*** Subroutine TRED2 (a, n, np, d, e)
* Householder reduction of a real symmetric, n by n matrix a, stored in a np by np physical array.
* On output, a is replaced by the orthogonal matrix q effecting the transformation. d returns the diagonal 
*  elements of the tridiagonal matrix, and e the off- (sub- or super-) diagonal elements, with e(1)=0. 
	   
	   call tred2 (sdm, 3, 3, diag, off)
	   
	   print *, 'out of TRED2'
	   print *, '            '
	   
	   
	   print *, 'now entering TQLI...'
	   
*** Subroutine TQLI (d, e, n, np, z)
* QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric, tridiagonal
* matrix, or of a real symmetric matrix previously reduced by tred2. d is a vector of length np. On input, its first n
* elements are the diagonal elements of the tridiagonal matrix. On output, it returns the eigenvalues. The vector e inputs
* the subdiagonal elements of the tridiagonal matrix with e(1) arbitrary. On output e is destroyed. If the eigenvectors of
* tridiagonal matrix are desired, the matrix z (n x n matrix, stored in np x np array) is input as the identity matrix.
* If the eigenvectors of a matrix that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
* In either case, the k-th column of z returns the normalized eigenvector corresponding to d(k).
	   
	   call tqli (diag, off, 3, 3, sdm)
	   
	   print *, 'out of TQLI'
	   print *, '            '
	   
*** Formetted writing of the eigenvalues (w**2.0)
	   write(47,'(f8.5,3(1x,1pe16.9))') 
     >		wkx/2./pi, (dsqrt(diag(kk))/2d12/pi, kk=1,3)
	   
	   do kk = 1, 3
	      write(48,*) (sdm(kk,mm), mm = 1, 3)
	   enddo
	   
*** Closure of the k-loop
	enddo
	
	stop
	
	
	end
	
