c meam.inc
c$Author: baskes $
c$Date: 2004/02/24 16:11:17 $
c$Source: /home/baskes/programs/MEAM/v88/RCS/meam.inc,v $
c$Revision: 1.7 $
c$Log: meam.inc,v $
cRevision 1.7  2004/02/24 16:11:17  baskes
cadded ialloy
c
cRevision 1.6  2003/12/06 04:19:54  baskes
cupdating
c
cRevision 1.5  2001/09/09 14:28:19  baskes
cadded SiO2
c
cRevision 1.4  2001/01/08 23:31:35  baskes
cfixed header
c
cRevision 1.3  2001/01/08 23:03:28  baskes
cfixed header
c
       parameter (sconst=1.3, scnres=1.0e-6 )
	real legend
      character*4 lattces
      character*4 enames
	common /pl/plocal(natmax)
      common /cmeam/ omegas(nelmax),res(nelmax,nelmax),
     .               esubs(nelmax,nelmax),alats(nelmax),zs(nelmax),
     .               alphas(nelmax,nelmax),deltas(nelmax,nelmax),
     .               betas(4,nelmax),ts(4,nelmax),asubs(nelmax),
     .               rozros(nelmax),z0s(nelmax),
     .	             repuls(nelmax,nelmax),attrac(nelmax,nelmax),
     .               legend,ibarr(nelmax),irhos(nelmax),
     .               enames(nelmax),lattces(nelmax)
        common/scr/ rscrmax,cmin(nelmax,nelmax,nelmax),
     1              cmax(nelmax,nelmax,nelmax),nscree,ipscree
        common /wrkspc1/ c(natmax),cg(natmax),a(3,natmax),
     .  b(3,3,natmax),ag(3,natmax),d(3,3,3,natmax),
     .  scrnab(natmax,neimax),dscrnab(3,natmax,neimax),tav(3,natmax),
     .  rhsq(3,natmax),tsqav(3,natmax),ialloy
	common /rhon/ irho,nrho
        common /myoutput/ fid
        integer fid
