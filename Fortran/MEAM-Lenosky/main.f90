
! Main program to test silicon potential subroutines
        implicit real*8 (a-h,o-z)
	character*8 potential
        integer count1,count2,count_rate,count_max
! ncell*: number of cubic unit cell along x,y,z
        parameter(ncell1=8,ncell2=8,ncell3=8)
! nat: number of atoms (e.g. vacancy)
        parameter(nat=8*ncell1*ncell2*ncell3-1)
        parameter(nint=2*32+1)
        dimension rxyz0(3,nat),fxyz(3,nat),ener_ind(nat),alat(3)
        dimension displ(3,nat),simpson(nint)
!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface
!$omp parallel private(iam,npr)  
!$       iam=omp_get_thread_num()
!$       npr=omp_get_num_threads()
!$       if (iam.eq.0) write(6,*) 'number of threads',npr
!$omp end parallel 

! Bazant potential
!	potential='bazant'
! Lenosky potential
	potential='lenosky'

	write(6,*) 'Testing ',potential,' potential'


	if (nat.le.0 .or. nat.gt. 8*ncell1*ncell2*ncell3) stop 'nat'

        if (mod(nint,2).ne.1) stop 'nint has to be odd'
        simpson(1)=1.d0/3.d0
        simpson(2)=4.d0/3.d0
        do i=3,nint-2,2
        simpson(i)=2.d0/3.d0
        simpson(i+1)=4.d0/3.d0
        enddo
        simpson(nint)=1.d0/3.d0

! create silicon diamond structure
        call diamond(alat,ncell1,ncell2,ncell3,nat,rxyz0)
        write(6,'(a,x,i7,3(x,e12.5))') 'nat,alat',nat,alat
        write(6,*) '  '
        count=0.d0
! create random displacements (use sin() instead of rand())
        stepsize=1.d-4
        do iat=1,nat
        displ(1,iat)=stepsize*abs(sin(iat+.2d0))
        displ(2,iat)=stepsize*abs(sin(iat+.4d0))
        displ(3,iat)=stepsize*abs(sin(iat+.7d0))
        enddo

! calculate energy at equilibrium geometry (diamond structure) 
! and at two additional points along the displacements
       call cpu_time(t1)
       call system_clock(count1,count_rate,count_max)
        path=0.d0
        do 1000,idispl=1,nint
        do 10,irep=1,1
	if (potential.eq.'bazant') then
        call bazant(nat,alat,rxyz0,fxyz,ener,coord,ener_var,coord_var,count)
	else if (potential.eq.'lenosky') then
        call lenosky(nat,alat,rxyz0,fxyz,ener_ind,ener,coord,ener_var,coord_var,count)
	else 
	stop 'unknown potential'
	endif
10      continue
        if (idispl.eq.1) ener0=ener

        if (idispl.eq.1 .or. idispl.eq.nint) then
! check whether total force vanishes
        sumx=0.d0
        sumy=0.d0
        sumz=0.d0
! check whether sum of individual energy equal total energy
        sume=0.d0
        do 8483,iat=1,nat
            sume=sume+ener_ind(iat)
            sumx=sumx+fxyz(1,iat)
            sumy=sumy+fxyz(2,iat)
8483        sumz=sumz+fxyz(3,iat)
        write(6,*) 'ener=',ener,'  sume=',sume
        write(6,*) 'ener-sume=',ener-sume
        write(6,'(a,i5,2(x,e19.12))') 'idispl,ener,ener/nat',idispl,ener,ener/nat
        write(6,'(a,3(x,e12.5))') 'coord,ener_var,coord_var', &
                                   coord,ener_var,coord_var
        write(6,'(a,3(x,e10.3))') & 
          'Sum of x, y and z component of forces:',sumx,sumy,sumz
        endif


! integrate force*displacement 
        t1=0.d0
        t2=0.d0
        t3=0.d0
        do iat=1,nat
        t1=t1-fxyz(1,iat)*displ(1,iat) 
        t2=t2-fxyz(2,iat)*displ(2,iat) 
        t3=t3-fxyz(3,iat)*displ(3,iat)
        enddo
        path=path+simpson(idispl)*(t1+t2+t3)

! next positions along path
        do iat=1,nat
        rxyz0(1,iat)=rxyz0(1,iat)+displ(1,iat)
        rxyz0(2,iat)=rxyz0(2,iat)+displ(2,iat)
        rxyz0(3,iat)=rxyz0(3,iat)+displ(3,iat)
        enddo
1000        continue
       call cpu_time(t2)
       call system_clock(count2,count_rate,count_max)


! compare energy difference with  force*displacement to check correctness of forces
        dener=ener-ener0
        write(6,*) '                           '
        write(6,*) 'Check correctness of forces'
        write(6,*) 'Difference of total energies ',dener
        write(6,*) 'Integral force*displacement  ',path
        write(6,*) 'Difference ',path-dener
        write(6,*) '                           '
        write(6,*) 'number of force evaluations (count)',count
        write(6,*) '    CPU time ', t2-t1
        time=(count2-count1)/float(count_rate)
        write(6,*) 'elapsed time ', time
        write(6,*) 'time/count, time/(count*nat)',time/count, time/(count*nat)


        end


      subroutine diamond(alat,ncell1,ncell2,ncell3,nat,rxyz)
!     computes initial positions of the atoms in diamond structure
        implicit real*8 (a-h,o-z)
        parameter(acell=5.42981d0)
        dimension rxyz(3,nat),alat(3)

        alat(1)=acell*ncell1
        alat(2)=acell*ncell2
        alat(3)=acell*ncell3

        acell1=acell
        acell2=acell
        acell3=acell
!      perfect lattice position routine
      do 30 j3=0,ncell3-1
       do 30 j2=0,ncell2-1
        do 30 j1=0,ncell1-1
           jj=8*(j1+ncell1*j2+ncell2*ncell1*j3)

           if (jj+1.le.nat) then
           rxyz(1,jj+1)=acell1*j1
           rxyz(2,jj+1)=acell2*j2
           rxyz(3,jj+1)=acell3*j3 
           endif

           if (jj+2.le.nat) then
           rxyz(1,jj+2)=acell1*j1 +.5d0*acell1
           rxyz(2,jj+2)=acell2*j2 +.5d0*acell2
           rxyz(3,jj+2)=acell3*j3 
           endif

           if (jj+3.le.nat) then
           rxyz(1,jj+3)=acell1*j1 +.5d0*acell1
           rxyz(2,jj+3)=acell2*j2 
           rxyz(3,jj+3)=acell3*j3 +.5d0*acell3
           endif

           if (jj+4.le.nat) then
           rxyz(1,jj+4)=acell1*j1
           rxyz(2,jj+4)=acell2*j2 +.5d0*acell2
           rxyz(3,jj+4)=acell3*j3 +.5d0*acell3
           endif

           if (jj+5.le.nat) then
           rxyz(1,jj+5)=acell1*j1 + .25d0*acell1
           rxyz(2,jj+5)=acell2*j2 + .25d0*acell2 
           rxyz(3,jj+5)=acell3*j3 + .25d0*acell3
           endif

           if (jj+6.le.nat) then
           rxyz(1,jj+6)=acell1*j1 + .25d0*acell1 +.5d0*acell1
           rxyz(2,jj+6)=acell2*j2 + .25d0*acell2 +.5d0*acell2
           rxyz(3,jj+6)=acell3*j3 + .25d0*acell3
           endif

           if (jj+7.le.nat) then
           rxyz(1,jj+7)=acell1*j1 + .25d0*acell1 +.5d0*acell1
           rxyz(2,jj+7)=acell2*j2 + .25d0*acell2  
           rxyz(3,jj+7)=acell3*j3 + .25d0*acell3 +.5d0*acell3
           endif

           if (jj+8.le.nat) then
           rxyz(1,jj+8)=acell1*j1 + .25d0*acell1
           rxyz(2,jj+8)=acell2*j2 + .25d0*acell2 +.5d0*acell2
           rxyz(3,jj+8)=acell3*j3 + .25d0*acell3 +.5d0*acell3
           endif
30    continue

      return
      end            

!  end of the test program part
