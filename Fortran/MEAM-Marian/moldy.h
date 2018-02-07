***************************************************************
*   MOLDY.H
*   Include file
*****************************************************************
* number of link cells on each direction
* total number of link cells

        implicit double precision(a-h,o-z)

        character(100) meamfile
        common/meamdatafile/meamfile

	parameter (pi = 3.1415926)
	integer nlcx, nlcy, nlcz, nlc
C        parameter (nlcx = 4, nlcy = 4, nlcz = 4)
C	parameter (nlc = nlcx*nlcy*nlcz)

        common/linklist/ nlcx, nlcy, nlcz, nlc
        parameter (nlcmax = 100000)

C	real*8 alatt
C        real*8 b0
C	parameter (alatt = 4.640)
C	parameter (b0 = 4d0)

* PBC Box
        real*8 Lx, Ly, Lz
        common/box/ Lx, Ly, Lz
    
	parameter (ndim = 3)

* number of atoms in the lattice,
* average number of neighbors.
* total number of atoms in a perfect lattice
*

C        integer nm, nnbrs
C        parameter (nm = 256, nnbrs = 400)

C nmmax:     maximum number of particles
C nm:        actual number of particles
        integer nm, nmmax, nnbrs
        common/NP/nm
        parameter (nmmax=100000, nnbrs=400)
    
*
* link cell map
*
        integer*4 nix(27),niy(27),niz(27)
        data nix/ 0,-1,-1,-1,0,0,-1,1,-1, 0, 1,-1,0,1, 0,1,1, 1, 0,
     &            0, 0, 1, 1, 1,-1,-1,-1/
        data niy/ 0, 0,-1, 1,1,0, 0,0,-1,-1,-1, 1,1,1,-1,0,1,-1, 1,
     &            0,-1, 0,-1, 1, 0,-1, 1/
        data niz/ 0, 0, 0, 0,0,1, 1,1, 1, 1, 1, 1,1,1, 0,0,0, 0,-1,
     &           -1,-1,-1,-1,-1,-1,-1,-1/
*
* Adimensional arrays in unitary box
	double precision x0(nmmax), y0(nmmax), z0(nmmax)
	common/adim/ x0, y0, z0

* Relaxed positions array
	double precision xi(ndim,nmmax), x(ndim,nmmax)
	common/pos/ xi, x
    
*
* Auxiliary arrays for 'rhofs.f' 'krafs.f'
* embedding energies,
	double precision rhotot(nmmax), embf(nmmax), embfp(nmmax)
        common/embed/ rhotot, embf, embfp

* Force arrays
        double precision f(ndim,nmmax), fplus(ndim,nmmax),
     >                   fmin(ndim,nmmax)
	common/force/ fplus, fmin, f

* potential energy of each atom
        double precision atpe(nmmax), atpe3b(nmmax), pe, petrip
        common/ppote/ atpe, atpe3b, pe, petrip

* neighbors position
	double precision xnbr(nnbrs), ynbr(nnbrs), znbr(nnbrs)
        common/area18/ xnbr, ynbr, znbr

* link map calculation    
	integer ltop(nlcmax), jaddr(nnbrs), linkmp(nmmax)
	common/lnmap/ ltop, jaddr, linkmp

* force in krafs.F
	double precision fx(nmmax),fy(nmmax),fz(nmmax)
        common/area6/ fx,fy,fz

*****************************************************************
* Common blocks for MEAM potential
*
 	parameter (nnbrs2 = 500)
	real*8 legend
        integer*4 noscr, ncrys
	common/meamread/zsmeam,alphas,betas(4),esubs,asubs,ts(4),rozros,
     >	rcutmeam,cmin,cmax,repuls,attrac,legend,equidist
        common/meamread_int/ibarr,noscr,ncrys

	common/meamforces/tav(3,nmmax),c8a(nmmax),dang2(nmmax),dang1(nmmax),
     >	 rhsq(3,nmmax),cg8c(nmmax),a8b(3,nmmax),ag(3,nmmax),
     >	 b8c(3,3,nmmax),d8d(3,3,3,nmmax),res,rcutmeam2
* Data	for screening function
        double precision scrnab, dscrnab
	common/screendata/scrnab(nmmax,nnbrs2)
	common/screendife/dscrnab(3,nmmax,nnbrs2)
*
*   New list for neighbors
*
 	common/neighborlist/numneigh(nmmax),neighlabel(nmmax,nnbrs2)
     >  !,neighdist(nmmax,nnbrs2)
*
***********************************************************************
	integer id(nmmax)
	common/identifier/ id
