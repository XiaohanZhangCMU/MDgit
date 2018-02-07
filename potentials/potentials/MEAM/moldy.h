***************************************************************
*   MOLDY.H
*   Include file
*****************************************************************
* number of link cells on each direction
* total number of link cells

        implicit double precision(a-h,o-z)

	parameter (pi = 3.1415926)
	integer nlcx, nlcy, nlcz, nlc
        parameter (nlcx = 4, nlcy = 4, nlcz = 4)
	parameter (nlc = nlcx*nlcy*nlcz)

	real*8 alatt, b0
	parameter (alatt = 4.640)
	parameter (b0 = 4d0)
	parameter (ndim = 3)

* number of atoms in the lattice,
* average number of neighbors.
* total number of atoms in a perfect lattice
*
	integer nm, nnbrs
        parameter (nm = 256, nnbrs = 400)
        integer*4 nix(27),niy(27),niz(27)
*
* link cell map
*
        data nix/ 0,-1,-1,-1,0,0,-1,1,-1, 0, 1,-1,0,1, 0,1,1, 1, 0,
     &            0, 0, 1, 1, 1,-1,-1,-1/
        data niy/ 0, 0,-1, 1,1,0, 0,0,-1,-1,-1, 1,1,1,-1,0,1,-1, 1,
     &            0,-1, 0,-1, 1, 0,-1, 1/
        data niz/ 0, 0, 0, 0,0,1, 1,1, 1, 1, 1, 1,1,1, 0,0,0, 0,-1,
     &           -1,-1,-1,-1,-1,-1,-1,-1/
*
* Adimensional arrays in unitary box
	double precision x0(nm), y0(nm), z0(nm)
	common/adim/ x0, y0, z0

* Relaxed positions array
	double precision xi(ndim,nm), x(ndim,nm)
	common/pos/ xi, x
*
* Auxiliary arrays for 'rhofs.f' 'krafs.f'
* embedding energies,
	double precision rhotot(nm), embf(nm), embfp(nm)
        common/embed/ rhotot, embf, embfp

* Force arrays
        double precision f(ndim,nm), fplus(ndim,nm), fmin(ndim,nm)
	common/force/ fplus, fmin, f

* potential energy of each atom
        double precision atpe(nm), atpe3b(nm), pe, petrip
        common/ppote/ atpe, atpe3b, pe, petrip

* neighbors position
	double precision xnbr(nnbrs), ynbr(nnbrs), znbr(nnbrs)
        common/area18/ xnbr, ynbr, znbr

* link map calculation    
	integer ltop(nlc), jaddr(nnbrs), linkmp(nm)
	common/lnmap/ ltop, jaddr, linkmp

* force in krafs.F
	double precision fx(nm),fy(nm),fz(nm)
        common/area6/ fx,fy,fz

*****************************************************************
* Common blocks for MEAM potential
*
 	parameter (nnbrs2 = 500)
	real*8 legend
	integer*4 noscr
	common/meamread/zsmeam,alphas,betas(4),esubs,asubs,ts(4),rozros,
     >	rcutmeam,cmin,cmax,repuls,attrac,legend
	common/meamread_int/ibarr,noscr

	common/meamforces/tav(3,nm),c8a(nm),dang2(nm),dang1(nm),
     >	 rhsq(3,nm),cg8c(nm),a8b(3,nm),ag(3,nm),
     >	 b8c(3,3,nm),d8d(3,3,3,nm),res,rcutmeam2
* Data	for screening function
        double precision scrnab, dscrnab
	common/screendata/scrnab(nm,nnbrs2)
	common/screendife/dscrnab(3,nm,nnbrs2)
*
*   New list for neighbors
*
 	common/neighborlist/numneigh(nm),neighlabel(nm,nnbrs2)	 !,neighdist(nm,nnbrs2)
*
***********************************************************************
	integer id(nm)
	common/identifier/ id
