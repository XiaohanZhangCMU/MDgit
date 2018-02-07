c utils.inc
c$Author: baskes $
c$Date: 2003/12/05 22:03:25 $
c$Source: /home/baskes/programs/EAM/dyn/utils/RCS/utils.inc,v $
c$Revision: 1.2 $
c$Log: utils.inc,v $
cRevision 1.2  2003/12/05 22:03:25  baskes
cupdating
c
cRevision 1.1  2003/08/20 17:17:42  baskes
cInitial revision
c
cRevision 1.10  2002/04/26 20:56:00  baskes
cadd vcm averaging
c
cRevision 1.9  2002/01/04 23:15:33  baskes
cmodified averaging
c
cRevision 1.8  2002/01/03 23:25:21  baskes
cadded sheardot
c
cRevision 1.7  2001/06/12 15:34:25  baskes
cadded acf
c
cRevision 1.6  2001/05/15 21:34:18  baskes
cadded table
c
cRevision 1.5  2001/01/17 01:03:33  baskes
caverage slocal added
c
cRevision 1.4  2001/01/08 23:31:47  baskes
cfixed header
c
cRevision 1.3  2001/01/08 23:02:36  baskes
cfixed header
c
      logical sortp,follow,icm
	character *8 latty
	common /vel/ temp0,ngen,iseed,icm
	common /interact/ rcutsq,rcut
      common /bigwrk/ dneighs(3,natmax*neimax/2),nneipss(natmax),
     . nneighs(natmax),neiind(0:natmax),jneighs(natmax*neimax/2)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      logical onlyoi,twoi,threei,fouri
      common /images/ onlyoi,twoi,threei,fouri
      common /particle/ rv(6,natmax),itype(natmax)
      common /scaled/ y(6,natmax)
      common /forces/ f(3,natmax),e(natmax),stresst(3,3),
     .                slocal(3,3,natmax),avslocal(3,3,natmax),noutput
      common /accels/ acc(3,natmax)
      common /density/ rho(natmax),fp(natmax)
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes
      common /integ/ t0,tend,dt,tol,dxmax,accbnd(2),nsteps,noutp,
     .               inte,iaccur,ndegfr,mnterm
      common /errs/ detops,gg
      common /calls/ frctmx,frctmn,gnetmx,gnetmn,ngtlst,nneimx,
     .               mxlstu,mxnnei,nforce
      common /cbound/ dpress,bndmas(3),bnddrg(3),bndvel(3),dstress(3),
     . ibdtyp,idynper(3),iper
      common /cneigh/ rctsqn,dradn,nneigh,nneips,nmeth
     .   ,nnindx(0:natmax),newlst,nlstmx,nnlst(natmax*neimax)
      dimension nlist3(neimax,natmax)
      equivalence (nnlst(1),nlist3(1,1))
      common /averag/ avtemp,avpres,avvol,avden,avpe,avke,avstrs(3,3),
     .  eqtim,avtmprg(0:6),avper(3),avper2(3),avtemp2,avpres2,
     .  avpe2,avke2,avvol2,avstrs2(3,3),avforce(0:6,3),avvcm(0:6,3),
     .  nequil,nskipd,nave
      common /fixit/ timefx,vectfx(3,0:100),dvctfx(3,0:100),
     .               workfx(natmax),posfx(3,natmax),modefx(0:100),
     .               ipntfx(natmax),nfixes,ifxal1,ifxal2
      common /tmpfix/ qreg(6),tmpreg(0:6),xcenter(6),xwidth(6),
     . e0,pzz,v0,
     . ycenter(6),ywidth(6),zcenter(6),zwidth(6),acttmp(0:6),
     . actvcm(0:6,3),actforce(0:6,3),
     . drag(0:6),timetmp,ifxtmp,nregs,ipntrg(natmax),natrg(0:6),
     . itempg(0:6,3),vcmrg(0:6,3),num1(6),num2(6),follow
      common /xcell/ mcell(natcmax,icmax,icmay,icmaz),natc(icmax,
     . icmay,icmaz),neigh(neimx)
        common/mono/ beta,betar,cbeta,sbeta,ctbeta,sheardot
      common /acfs/ nbuff, nval, nlimitacf, ncalcacf
      common /print/ typep,xminp,xmaxp,yminp,ymaxp,zminp,zmaxp,sortp,
     . ipinter,ipatoms,ipitera,irstrt,iconf,iconst,nminp,nmaxp,iforce
