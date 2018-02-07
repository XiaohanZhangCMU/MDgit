/*
  rebo.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Apr 22 16:40:43 2010

  FUNCTION  :  MD simulation package of C using REBO potential
  Reference:   

  Note:        Implemented based on LAMMPS (7/24/2009)

  Atom species: (need to check with LAMMPS manual)
               0 - C
               1 - H
  To do:
        1. translate MD++ neighbor list to LAMMPS neighbor list (DONE)
        2. convert LAMMPS tally_# functions to MD++ virial/force accumulation (DONE)
        3. construct a function to read a potential file and parameters. (DONE)
           ex.rcmin,rcmax,rcmaxsq 
        
*/

#include "rebo.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024
#define TOL 1.0e-9
#define PGDELTA 1

void REBOFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int size;
    size=_NP*allocmultiple;

    Realloc(REBO_numneigh,int,size);
    Realloc(closestdistsq,double,size);
    Realloc(nC,double,size);
    Realloc(nH,double,size);

}

void REBOFrame::initvars()
{
    MDPARALLELFrame::initvars();

    strcpy(incnfile,"../carbon.cn");
    strcpy(rebofile,"rebof");

    REBO_numneigh = NULL;
    REBO_firstneigh = NULL;
    REBO_firstneigh_mem = NULL;
    closestdistsq = NULL;
    nC = nH = NULL;

    DUMP(HIG"REBOFrame initvars"NOR);
}

void REBOFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("rebofile",rebofile,STRING);
    bindvar("rcut",&acut,DOUBLE);
}

int REBOFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readREBO",readAIREBO());

    return -1;
}

void REBOFrame::rebo()
{
  /* LAMMPS variables */
  int j,m,itype,jtype;
  double evdwl,ecoul,fpair;
  double r2ij,rrij,wij;
  double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
  double dwij;

  int vflag_atom=1;

  /* MD++ variables */
  int ipt, jpt;
  int n0, n1;
  Vector3 sij,rij,fij;

  DUMP("REBO");
    
  refreshneighborlist();
  REBO_neigh();

  /* Move these to the beginning of potential() function, 
     so that rebo() can be called multiple times in rebolj.cpp */
  /*
  _EPOT=0; _VIRIAL.clear();
  for(ipt=0;ipt<_NP;ipt++) { _F[ipt].clear(); _EPOT_IND[ipt]=0; }
  */

  n0=0;
  n1=_NP;

  evdwl = 0.0;
  ecoul = 0.0;

  // two-body interactions from REBO neighbor list, skip half of them
  for(ipt=n0;ipt<n1;ipt++)
  {
    if(fixed[ipt]==-1) continue; /* ignore atom if -1 */
    itype = species[ipt];

    /* move all neighboring atoms to the nearest image around atom ipt */
    _R[ipt] = _H*_SR[ipt];
    for(j=0;j<nn[ipt];j++)
    {
      jpt=nindex[ipt][j];
      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
      sij.subtract(_SR[ipt],_SR[jpt]);
      sij.subint();
      rij=_H*sij;
      _R[jpt]=_R[ipt]-rij;
    }

    for(j=0;j<nn[ipt];j++)
    {
      jpt=nindex[ipt][j];
      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
    
      /* skip if ipt >= jpt */
      if (ipt>=jpt) continue;

      jtype = species[jpt];

      rij=_R[ipt]-_R[jpt];  
      r2ij=rij.norm2();
      if(r2ij>acutsq) continue;
      rrij=sqrt(r2ij);

      wij = Sp(rrij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
      if (wij <= TOL) continue;
      
      Qij = Q[itype][jtype];
      Aij = A[itype][jtype];
      alphaij = alpha[itype][jtype];

      VR = wij*(1.0+(Qij/rrij)) * Aij*exp(-alphaij*rrij);
      pre = wij*Aij * exp(-alphaij*rrij);
      dVRdi = pre * ((-alphaij)-(Qij/r2ij)-(Qij*alphaij/rrij));
      dVRdi += VR/wij * dwij;

      VA = dVA = 0.0;
      for (m = 0; m < 3; m++) {
	term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rrij);
	VA += term;
	dVA += -Beta[itype][jtype][m] * term;
      }
      dVA += VA/wij * dwij;

      bij = bondorder(ipt,jpt,&rij[0],rrij,VA,_F,vflag_atom);
      dVAdi = bij*dVA;

      fpair = -(dVRdi+dVAdi) / rrij;

      fij.clear(); fij.addnv(fpair,rij);
      _F[ipt] += fij;
      _F[jpt] -= fij;

      if (eflag)  evdwl = VR + bij*VA;
      if (evflag) ev_tally(ipt,jpt,evdwl,ecoul,fpair,rij.x,rij.y,rij.z);

    }
  }

  /* zero out the forces on fixed atoms */
//  for(ipt=0;ipt<_NP;ipt++)
//  {
//     if(fixed[ipt]) _F[ipt].clear();
//  }

  /* restore _R to avoid unnecessary update of neighbor list */
  SHtoR();
}

//void REBOFrame::NbrList_reconstruct(int iatom)
//{
//    MDPARALLELFrame::NbrList_reconstruct();
//}

/* ----------------------------------------------------------------------
   create REBO neighbor list from main neighbor list
   REBO neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void REBOFrame::REBO_neigh()
{
  /* LAMMPS variables */
  int j,n,itype,jtype;
  double dS;

  /* MD++ variables */
  int i, ipt, jpt;
  int n0, n1;
  Vector3 sij, rij;
  double r2ij;

  n0=0;
  n1=_NP;

  Realloc(REBO_firstneigh_mem,char,(n1-n0)*_NNM*sizeof(int)+(n1-n0)*sizeof(int *));
  REBO_firstneigh=(int **)(REBO_firstneigh_mem+(n1-n0)*_NNM*sizeof(int));
  for(i=0;i<(n1-n0);i++) REBO_firstneigh[i]=(int *)(REBO_firstneigh_mem+i*_NNM*sizeof(int));
  
  // store all REBO neighs of owned atoms
  // scan full neighbor list of I
  // if J is ghost and within LJ cutoff:
  //   flag it via REBO_numneigh so its REBO neighbors will be stored below
  //   REBO requires neighbors of neighbors of i,j in each i,j LJ interaction

  /* translate MD++ neighbor list to LAMMPS neighbor list */
  for (ipt = n0; ipt< n1; ipt++) {

    if(fixed[ipt]==-1) continue; /* ignore atom if -1 */

    itype = species[ipt];  // species of atom ipt

    n = 0;
    nC[ipt] = nH[ipt] = 0.0;   // initialize nC and nH to be 0.
    
    // iterate all the neighbors of an atom ipt
    for (j = 0; j < nn[ipt]; j++) {

      jpt=nindex[ipt][j];   // the index of j-th neighbor

      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */

      jtype = species[jpt]; // species of atom jpt
      
      // calculate (atomic distance)^2 of atoms ipt and jpt
      sij.subtract(_SR[jpt],_SR[ipt]);
      sij.subint();
      _H.multiply(sij,rij);
      r2ij=rij.norm2();

      // calculate weighting funcs nC[ipt] and nH[ipt]
      // depending on the combinations of species.
      if (r2ij < rcmaxsq[itype][jtype]) {
        REBO_firstneigh[ipt][n] = jpt; 
        n ++;
	if (jtype == 0)
	  nC[ipt] += Sp(sqrt(r2ij),rcmin[itype][jtype],rcmax[itype][jtype],dS);
	else
	  nH[ipt] += Sp(sqrt(r2ij),rcmin[itype][jtype],rcmax[itype][jtype],dS);
      }
    }
    REBO_numneigh[ipt] = n;    
  }

#if 0
  /* debug: print out LAMMPS neighbor list */
  for (ipt = n0; ipt < n1; ipt ++)
  {
      //if (ipt > 2) continue;
      INFO_Printf("nn[%d] = %d  REBO_numneigh[%d] = %d\n",ipt,nn[ipt],ipt,REBO_numneigh[ipt]);
      INFO_Printf("REBO_firstneigh[%d] = ", ipt);
      for (j = 0; j < REBO_numneigh[ipt]; j ++ ) INFO_Printf(" %d,",REBO_firstneigh[ipt][j]);
      INFO_Printf("\n");
      //INFO_Printf("nindex[%d] = ", ipt);
      //for (j = 0; j < nn[ipt]; j ++ ) INFO_Printf(" %d,",nindex[ipt][j]); INFO_Printf("\n");
      //INFO_Printf("REBO_firstneigh[%d] = %x", ipt, REBO_firstneigh[ipt]);
  }
#endif
      
}
    
/* ----------------------------------------------------------------------
   Bij function
------------------------------------------------------------------------- */
double REBOFrame::bondorder(int i, int j, double rij[3],
			     double rijmag, double VA,
			     Vector3 *f, int vflag_atom)
{
  int atomi,atomj,k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int itype,jtype,ktype,ltype,ntype;
  double rik[3],rjl[3],rkn[3],rji[3],rki[3],rlj[3],rknmag,dNki,dwjl,bij;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ,Nki,Nlj,dS;
  double lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double Tij;
  double r32[3],r32mag,cos321,r43[3],r13[3];
  double dNlj;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  //double **x = atom->x;
  Vector3 *x = _R;

  atomi = i;
  atomj = j;

  //int *type = atom->type;
  //itype = map[type[i]];
  //jtype = map[type[j]];
  itype = species[i];
  jtype = species[j];

  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijC = nC[i]-(wij*kronecker(jtype,0));
  NijH = nH[i]-(wij*kronecker(jtype,1));
  NjiC = nC[j]-(wij*kronecker(itype,0));
  NjiH = nH[j]-(wij*kronecker(itype,1));
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  
  REBO_neighs = REBO_firstneigh[i];

  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs[k];

    if (atomk != atomj) {
      ktype = species[atomk];

      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];

      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) * 
	((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));

      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
	(wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) / 
	(rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg
      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);

      Etmp = Etmp+(wik*g*exp(lamdajik)); 
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS));
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  pij = pow(1.0+Etmp+PijS,-0.5); 
  tmp = -0.5*pow(pij,3.0);

  // pij forces

  REBO_neighs = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs[k];
    if (atomk != atomj) {
      //ktype = map[type[atomk]];
      ktype = species[atomk];

      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) * 
	((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      cosjik = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) / 
	(rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      dcosjikdri[0] = ((rij[0]+rik[0])/(rijmag*rikmag)) - 
	(cosjik*((rij[0]/(rijmag*rijmag))+(rik[0]/(rikmag*rikmag))));
      dcosjikdri[1] = ((rij[1]+rik[1])/(rijmag*rikmag)) - 
	(cosjik*((rij[1]/(rijmag*rijmag))+(rik[1]/(rikmag*rikmag))));
      dcosjikdri[2] = ((rij[2]+rik[2])/(rijmag*rikmag)) - 
	(cosjik*((rij[2]/(rijmag*rijmag))+(rik[2]/(rikmag*rikmag))));
      dcosjikdrk[0] = (-rij[0]/(rijmag*rikmag)) + 
	(cosjik*(rik[0]/(rikmag*rikmag)));
      dcosjikdrk[1] = (-rij[1]/(rijmag*rikmag)) + 
	(cosjik*(rik[1]/(rikmag*rikmag)));
      dcosjikdrk[2] = (-rij[2]/(rijmag*rikmag)) + 
	(cosjik*(rik[2]/(rikmag*rikmag)));
      dcosjikdrj[0] = (-rik[0]/(rijmag*rikmag)) + 
	(cosjik*(rij[0]/(rijmag*rijmag)));
      dcosjikdrj[1] = (-rik[1]/(rijmag*rikmag)) + 
	(cosjik*(rij[1]/(rijmag*rijmag)));
      dcosjikdrj[2] = (-rik[2]/(rijmag*rikmag)) + 
	(cosjik*(rij[2]/(rijmag*rijmag)));

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wik*dgdc*exp(lamdajik));
      fj[0] = -tmp2*dcosjikdrj[0];
      fj[1] = -tmp2*dcosjikdrj[1];
      fj[2] = -tmp2*dcosjikdrj[2];
      fi[0] = -tmp2*dcosjikdri[0];
      fi[1] = -tmp2*dcosjikdri[1];
      fi[2] = -tmp2*dcosjikdri[2];
      fk[0] = -tmp2*dcosjikdrk[0];
      fk[1] = -tmp2*dcosjikdrk[1];
      fk[2] = -tmp2*dcosjikdrk[2];

      tmp2 = VA*.5*(tmp*wik*g*exp(lamdajik)*4.0*kronecker(itype,1));
      fj[0] -= tmp2*(-rij[0]/rijmag);
      fj[1] -= tmp2*(-rij[1]/rijmag);
      fj[2] -= tmp2*(-rij[2]/rijmag);
      fi[0] -= tmp2*((-rik[0]/rikmag)+(rij[0]/rijmag));
      fi[1] -= tmp2*((-rik[1]/rikmag)+(rij[1]/rijmag));
      fi[2] -= tmp2*((-rik[2]/rikmag)+(rij[2]/rijmag));
      fk[0] -= tmp2*(rik[0]/rikmag);
      fk[1] -= tmp2*(rik[1]/rikmag);
      fk[2] -= tmp2*(rik[2]/rikmag);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwik*g*exp(lamdajik))/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ktype]*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2]; 

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atomk][0] += fk[0]; f[atomk][1] += fk[1]; f[atomk][2] += fk[2];

      if (vflag_atom) {
	rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
	rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
	v_tally3(atomi,atomj,atomk,fj,fk,rji,rki);
      }
    }
  }

  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
      //ltype = map[type[atoml]];
      ltype = species[atoml];

      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) * 
	((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0)) + 
	nH[atoml]-(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) / 
	(rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl)); 
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl)); 
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);  
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*pow(pji,3.0);

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
      //ltype = map[type[atoml]];
      ltype = species[atoml];

      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) * 
	((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2]))) / 
	(rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      dcosijldri[0] = (-rjl[0]/(rijmag*rjlmag)) - 
	(cosijl*rij[0]/(rijmag*rijmag));
      dcosijldri[1] = (-rjl[1]/(rijmag*rjlmag)) - 
	(cosijl*rij[1]/(rijmag*rijmag));
      dcosijldri[2] = (-rjl[2]/(rijmag*rjlmag)) - 
	(cosijl*rij[2]/(rijmag*rijmag));
      dcosijldrj[0] = ((-rij[0]+rjl[0])/(rijmag*rjlmag)) + 
	(cosijl*((rij[0]/pow(rijmag,2.0))-(rjl[0]/(rjlmag*rjlmag))));
      dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) + 
	(cosijl*((rij[1]/pow(rijmag,2.0))-(rjl[1]/(rjlmag*rjlmag))));
      dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) + 
	(cosijl*((rij[2]/pow(rijmag,2.0))-(rjl[2]/(rjlmag*rjlmag))));
      dcosijldrl[0] = (rij[0]/(rijmag*rjlmag))+(cosijl*rjl[0]/(rjlmag*rjlmag));
      dcosijldrl[1] = (rij[1]/(rijmag*rjlmag))+(cosijl*rjl[1]/(rjlmag*rjlmag));
      dcosijldrl[2] = (rij[2]/(rijmag*rjlmag))+(cosijl*rjl[2]/(rjlmag*rjlmag));

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wjl*dgdc*exp(lamdaijl));
      fi[0] = -tmp2*dcosijldri[0]; 
      fi[1] = -tmp2*dcosijldri[1]; 
      fi[2] = -tmp2*dcosijldri[2]; 
      fj[0] = -tmp2*dcosijldrj[0]; 
      fj[1] = -tmp2*dcosijldrj[1]; 
      fj[2] = -tmp2*dcosijldrj[2]; 
      fl[0] = -tmp2*dcosijldrl[0]; 
      fl[1] = -tmp2*dcosijldrl[1]; 
      fl[2] = -tmp2*dcosijldrl[2]; 

      tmp2 = VA*.5*(tmp*wjl*g*exp(lamdaijl)*4.0*kronecker(jtype,1));
      fi[0] -= tmp2*(rij[0]/rijmag);
      fi[1] -= tmp2*(rij[1]/rijmag);
      fi[2] -= tmp2*(rij[2]/rijmag);
      fj[0] -= tmp2*((-rjl[0]/rjlmag)-(rij[0]/rijmag));
      fj[1] -= tmp2*((-rjl[1]/rjlmag)-(rij[1]/rijmag));
      fj[2] -= tmp2*((-rjl[2]/rjlmag)-(rij[2]/rijmag));
      fl[0] -= tmp2*(rjl[0]/rjlmag);
      fl[1] -= tmp2*(rjl[1]/rjlmag);
      fl[2] -= tmp2*(rjl[2]/rjlmag);

      // coordination forces
      
      // dwik forces

      tmp2 = VA*.5*(tmp*dwjl*g*exp(lamdaijl))/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // PIJ forces 

      tmp2 = VA*.5*(tmp*dN2[ltype]*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];  

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atoml][0] += fl[0]; f[atoml][1] += fl[1]; f[atoml][2] += fl[2];

      if (vflag_atom) {
	rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
	v_tally3(atomi,atomj,atoml,fi,fl,rij,rlj);
      }
    }
  }	
  
  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);
  
  // piRC forces

  REBO_neighs_i = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs_i[k];
    if (atomk !=atomj) {
      //ktype = map[type[atomk]];
      ktype = species[atomk];

      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] - 
	(wik*kronecker(itype,1));
      SpN = Sp(Nki,Nmin,Nmax,dNki);

      tmp2 = VA*dN3[0]*dwik/rikmag;
      f[atomi][0] -= tmp2*rik[0]; 
      f[atomi][1] -= tmp2*rik[1]; 
      f[atomi][2] -= tmp2*rik[2]; 
      f[atomk][0] += tmp2*rik[0]; 
      f[atomk][1] += tmp2*rik[1]; 
      f[atomk][2] += tmp2*rik[2]; 

      if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

      tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;
      f[atomi][0] -= tmp2*rik[0]; 
      f[atomi][1] -= tmp2*rik[1]; 
      f[atomi][2] -= tmp2*rik[2]; 
      f[atomk][0] += tmp2*rik[0]; 
      f[atomk][1] += tmp2*rik[1]; 
      f[atomk][2] += tmp2*rik[2];

      if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

      if (fabs(dNki) > TOL) {
	REBO_neighs_k = REBO_firstneigh[atomk];
	for (n = 0; n < REBO_numneigh[atomk]; n++) {
	  atomn = REBO_neighs_k[n];
	  if (atomn != atomi) {
	    //ntype = map[type[atomn]];
            ntype = species[atomn];

	    rkn[0] = x[atomk][0]-x[atomn][0];
	    rkn[1] = x[atomk][1]-x[atomn][1];
	    rkn[2] = x[atomk][2]-x[atomn][2];
	    rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
	    Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

	    tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)/rknmag;
	    f[atomk][0] -= tmp2*rkn[0]; 
	    f[atomk][1] -= tmp2*rkn[1]; 
	    f[atomk][2] -= tmp2*rkn[2]; 
	    f[atomn][0] += tmp2*rkn[0]; 
	    f[atomn][1] += tmp2*rkn[1]; 
	    f[atomn][2] += tmp2*rkn[2];

	    if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
	  }
	}
      } 
    }
  }  
  
  // piRC forces

  REBO_neighs = REBO_firstneigh[atomj];
  for (l = 0; l < REBO_numneigh[atomj]; l++) {
    atoml = REBO_neighs[l];
    if (atoml !=atomi) {
      //ltype = map[type[atoml]];
      ltype = species[atoml];

      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	(wjl*kronecker(jtype,1));
      SpN = Sp(Nlj,Nmin,Nmax,dNlj);

      tmp2 = VA*dN3[1]*dwjl/rjlmag;
      f[atomj][0] -= tmp2*rjl[0]; 
      f[atomj][1] -= tmp2*rjl[1]; 
      f[atomj][2] -= tmp2*rjl[2];
      f[atoml][0] += tmp2*rjl[0]; 
      f[atoml][1] += tmp2*rjl[1]; 
      f[atoml][2] += tmp2*rjl[2];

      if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;
      f[atomj][0] -= tmp2*rjl[0]; 
      f[atomj][1] -= tmp2*rjl[1]; 
      f[atomj][2] -= tmp2*rjl[2];
      f[atoml][0] += tmp2*rjl[0]; 
      f[atoml][1] += tmp2*rjl[1]; 
      f[atoml][2] += tmp2*rjl[2];

      if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      if (fabs(dNlj) > TOL) {
	REBO_neighs_l = REBO_firstneigh[atoml];
	for (n = 0; n < REBO_numneigh[atoml]; n++) {
	  atomn = REBO_neighs_l[n];
	  if (atomn != atomj) {
	    //ntype = map[type[atomn]];
            ntype = species[atomn];

	    rln[0] = x[atoml][0]-x[atomn][0];
	    rln[1] = x[atoml][1]-x[atomn][1];
	    rln[2] = x[atoml][2]-x[atomn][2];
	    rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
	    Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

	    tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)/rlnmag;
	    f[atoml][0] -= tmp2*rln[0]; 
	    f[atoml][1] -= tmp2*rln[1]; 
	    f[atoml][2] -= tmp2*rln[2];
	    f[atomn][0] += tmp2*rln[0]; 
	    f[atomn][1] += tmp2*rln[1]; 
	    f[atomn][2] += tmp2*rln[2];
	    
	    if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
	  }
	}
      } 
    }
  }   
  
  Tij = 0.0;
  dN3[0] = 0.0;
  dN3[1] = 0.0;
  dN3[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3);
  Etmp = 0.0;

  if (fabs(Tij) > TOL) {
    atom2 = atomi;
    atom3 = atomj;
    r32[0] = x[atom3][0]-x[atom2][0];
    r32[1] = x[atom3][1]-x[atom2][1];
    r32[2] = x[atom3][2]-x[atom2][2];
    r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
    r23[0] = -r32[0];
    r23[1] = -r32[1];
    r23[2] = -r32[2];
    r23mag = r32mag;
    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      atom1 = atomk;
      //ktype = map[type[atomk]];
      ktype = species[atomk];

      if (atomk != atomj) {
	r21[0] = x[atom2][0]-x[atom1][0];
	r21[1] = x[atom2][1]-x[atom1][1];
	r21[2] = x[atom2][2]-x[atom1][2];
	r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
	cos321 = -1.0*((r21[0]*r32[0])+(r21[1]*r32[1])+(r21[2]*r32[2])) / 
	  (r21mag*r32mag);
	cos321 = MIN(cos321,1.0);
	cos321 = MAX(cos321,-1.0);
	Sp2(cos321,thmin,thmax,dcut321);
	sin321 = sqrt(1.0 - cos321*cos321);
	sink2i = 1.0/(sin321*sin321);
	rik2i = 1.0/(r21mag*r21mag);
	if (sin321 != 0.0) {
	  rr = (r23mag*r23mag)-(r21mag*r21mag);
	  rjk[0] = r21[0]-r23[0];
	  rjk[1] = r21[1]-r23[1];
	  rjk[2] = r21[2]-r23[2];
	  rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
	  rijrik = 2.0*r23mag*r21mag;
	  rik2 = r21mag*r21mag;
	  dctik = (-rr+rjk2)/(rijrik*rik2);
	  dctij = (rr+rjk2)/(rijrik*r23mag*r23mag);
	  dctjk = -2.0/rijrik;
	  w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
	  rijmag = r32mag;
	  rikmag = r21mag;
	  rij2 = r32mag*r32mag;
	  rik2 = r21mag*r21mag;
	  costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
	  tspjik = Sp2(costmp,thmin,thmax,dtsjik);
	  dtsjik = -dtsjik;

	  REBO_neighs_j = REBO_firstneigh[j];
	  for (l = 0; l < REBO_numneigh[j]; l++) {
	    atoml = REBO_neighs_j[l];
	    atom4 = atoml;
	    //ltype = map[type[atoml]];
            ltype = species[atoml];

	    if (!(atoml == atomi || atoml == atomk)) {
	      r34[0] = x[atom3][0]-x[atom4][0];
	      r34[1] = x[atom3][1]-x[atom4][1];
	      r34[2] = x[atom3][2]-x[atom4][2];
	      r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
	      cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) / 
		(r32mag*r34mag);
	      cos234 = MIN(cos234,1.0);
	      cos234 = MAX(cos234,-1.0);
	      sin234 = sqrt(1.0 - cos234*cos234);
	      sinl2i = 1.0/(sin234*sin234);
	      rjl2i = 1.0/(r34mag*r34mag);

	      if (sin234 != 0.0) {
		w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
		rr = (r23mag*r23mag)-(r34mag*r34mag);
		ril[0] = r23[0]+r34[0];
		ril[1] = r23[1]+r34[1];
		ril[2] = r23[2]+r34[2];
		ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
		rijrjl = 2.0*r23mag*r34mag;
		rjl2 = r34mag*r34mag;
		dctjl = (-rr+ril2)/(rijrjl*rjl2);
		dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
		dctil = -2.0/rijrjl;
		rjlmag = r34mag;
		rjl2 = r34mag*r34mag;
		costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
		tspijl = Sp2(costmp,thmin,thmax,dtsijl);
		dtsijl = -dtsijl;
		prefactor = VA*Tij;

		cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
		cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
		cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
		cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
		cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
		cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

		cwnum = (cross321[0]*cross234[0]) + 
		  (cross321[1]*cross234[1]) + (cross321[2]*cross234[2]);
		cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
		om1234 = cwnum/cwnom;
		cw = om1234;
		Etmp += ((1.0-pow(om1234,2.0))*w21*w34) * 
		  (1.0-tspjik)*(1.0-tspijl);
		
		dt1dik = (rik2i)-(dctik*sink2i*cos321);
		dt1djk = (-dctjk*sink2i*cos321);
		dt1djl = (rjl2i)-(dctjl*sinl2i*cos234);
		dt1dil = (-dctil*sinl2i*cos234);
		dt1dij = (2.0/(r23mag*r23mag))-(dctij*sink2i*cos321) - 
		  (dctji*sinl2i*cos234);
		
		dt2dik[0] = (-r23[2]*cross234[1])+(r23[1]*cross234[2]);
		dt2dik[1] = (-r23[0]*cross234[2])+(r23[2]*cross234[0]);
		dt2dik[2] = (-r23[1]*cross234[0])+(r23[0]*cross234[1]);
		
		dt2djl[0] = (-r23[1]*cross321[2])+(r23[2]*cross321[1]);
		dt2djl[1] = (-r23[2]*cross321[0])+(r23[0]*cross321[2]);
		dt2djl[2] = (-r23[0]*cross321[1])+(r23[1]*cross321[0]);
		
		dt2dij[0] = (r21[2]*cross234[1])-(r34[2]*cross321[1]) - 
		  (r21[1]*cross234[2])+(r34[1]*cross321[2]);
		dt2dij[1] = (r21[0]*cross234[2])-(r34[0]*cross321[2]) - 
		  (r21[2]*cross234[0])+(r34[2]*cross321[0]);
		dt2dij[2] = (r21[1]*cross234[0])-(r34[1]*cross321[0]) - 
		  (r21[0]*cross234[1])+(r34[0]*cross321[1]);
		
		aa = (prefactor*2.0*cw/cwnom)*w21*w34 * 
		  (1.0-tspjik)*(1.0-tspijl);
		aaa1 = -prefactor*(1.0-pow(om1234,2.0)) * 
		  (1.0-tspjik)*(1.0-tspijl);
		aaa2 = aaa1*w21*w34;
		at2 = aa*cwnum;
		
		fcijpc = (-dt1dij*at2)+(aaa2*dtsjik*dctij*(1.0-tspijl)) + 
		  (aaa2*dtsijl*dctji*(1.0-tspjik));
		fcikpc = (-dt1dik*at2)+(aaa2*dtsjik*dctik*(1.0-tspijl));
		fcjlpc = (-dt1djl*at2)+(aaa2*dtsijl*dctjl*(1.0-tspjik));
		fcjkpc = (-dt1djk*at2)+(aaa2*dtsjik*dctjk*(1.0-tspijl));
		fcilpc = (-dt1dil*at2)+(aaa2*dtsijl*dctil*(1.0-tspjik));
		
		F23[0] = (fcijpc*r23[0])+(aa*dt2dij[0]);
		F23[1] = (fcijpc*r23[1])+(aa*dt2dij[1]);
		F23[2] = (fcijpc*r23[2])+(aa*dt2dij[2]);
		
		F12[0] = (fcikpc*r21[0])+(aa*dt2dik[0]);
		F12[1] = (fcikpc*r21[1])+(aa*dt2dik[1]);
		F12[2] = (fcikpc*r21[2])+(aa*dt2dik[2]);
		
		F34[0] = (fcjlpc*r34[0])+(aa*dt2djl[0]);
		F34[1] = (fcjlpc*r34[1])+(aa*dt2djl[1]);
		F34[2] = (fcjlpc*r34[2])+(aa*dt2djl[2]);
		
		F31[0] = (fcjkpc*rjk[0]);
		F31[1] = (fcjkpc*rjk[1]);
		F31[2] = (fcjkpc*rjk[2]);
		
		F24[0] = (fcilpc*ril[0]);
		F24[1] = (fcilpc*ril[1]);
		F24[2] = (fcilpc*ril[2]);
		
		f1[0] = -F12[0]-F31[0];
		f1[1] = -F12[1]-F31[1];
		f1[2] = -F12[2]-F31[2];
		f2[0] = F23[0]+F12[0]+F24[0];
		f2[1] = F23[1]+F12[1]+F24[1];
		f2[2] = F23[2]+F12[2]+F24[2];
		f3[0] = -F23[0]+F34[0]+F31[0];
		f3[1] = -F23[1]+F34[1]+F31[1];
		f3[2] = -F23[2]+F34[2]+F31[2];
		f4[0] = -F34[0]-F24[0];
		f4[1] = -F34[1]-F24[1];
		f4[2] = -F34[2]-F24[2];
		
		// coordination forces

		tmp2 = VA*Tij*((1.0-(om1234*om1234))) * 
		  (1.0-tspjik)*(1.0-tspijl)*dw21*w34/r21mag;
		f2[0] -= tmp2*r21[0];
		f2[1] -= tmp2*r21[1];
		f2[2] -= tmp2*r21[2];
		f1[0] += tmp2*r21[0];
		f1[1] += tmp2*r21[1];
		f1[2] += tmp2*r21[2];
		
		tmp2 = VA*Tij*((1.0-(om1234*om1234))) * 
		  (1.0-tspjik)*(1.0-tspijl)*w21*dw34/r34mag;
		f3[0] -= tmp2*r34[0];
		f3[1] -= tmp2*r34[1];
		f3[2] -= tmp2*r34[2];
		f4[0] += tmp2*r34[0];
		f4[1] += tmp2*r34[1];
		f4[2] += tmp2*r34[2];  

		f[atom1][0] += f1[0]; f[atom1][1] += f1[1];
		f[atom1][2] += f1[2];
		f[atom2][0] += f2[0]; f[atom2][1] += f2[1];
		f[atom2][2] += f2[2];
		f[atom3][0] += f3[0]; f[atom3][1] += f3[1];
		f[atom3][2] += f3[2];
		f[atom4][0] += f4[0]; f[atom4][1] += f4[1];
		f[atom4][2] += f4[2];
		
		if (vflag_atom) {
		  r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
		  r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
		  v_tally4(atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43);
		}
	      }
	    }
	  }
	}
      }
    }
    
    // Tij forces now that we have Etmp

    REBO_neighs = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs[k];
      if (atomk != atomj) {
	//ktype = map[type[atomk]];
        ktype = species[atomk];

	rik[0] = x[atomi][0]-x[atomk][0];
	rik[1] = x[atomi][1]-x[atomk][1];
	rik[2] = x[atomi][2]-x[atomk][2];
	rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
	wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
	Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] - 
	  (wik*kronecker(itype,1));
	SpN = Sp(Nki,Nmin,Nmax,dNki);

	tmp2 = VA*dN3[0]*dwik*Etmp/rikmag;
	f[atomi][0] -= tmp2*rik[0]; 
	f[atomi][1] -= tmp2*rik[1]; 
	f[atomi][2] -= tmp2*rik[2]; 
	f[atomk][0] += tmp2*rik[0]; 
	f[atomk][1] += tmp2*rik[1]; 
	f[atomk][2] += tmp2*rik[2]; 

	if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;
	f[atomi][0] -= tmp2*rik[0]; 
	f[atomi][1] -= tmp2*rik[1]; 
	f[atomi][2] -= tmp2*rik[2]; 
	f[atomk][0] += tmp2*rik[0]; 
	f[atomk][1] += tmp2*rik[1]; 
	f[atomk][2] += tmp2*rik[2]; 

	if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	if (fabs(dNki) > TOL) {
	  REBO_neighs_k = REBO_firstneigh[atomk];
	  for (n = 0; n < REBO_numneigh[atomk]; n++) {
	    atomn = REBO_neighs_k[n];
	    //ntype = map[type[atomn]];
            ntype = species[atomn];

	    if (atomn != atomi) {
	      rkn[0] = x[atomk][0]-x[atomn][0];
	      rkn[1] = x[atomk][1]-x[atomn][1];
	      rkn[2] = x[atomk][2]-x[atomn][2];
	      rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
	      Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

	      tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)*Etmp/rknmag;
	      f[atomk][0] -= tmp2*rkn[0]; 
	      f[atomk][1] -= tmp2*rkn[1]; 
	      f[atomk][2] -= tmp2*rkn[2]; 
	      f[atomn][0] += tmp2*rkn[0]; 
	      f[atomn][1] += tmp2*rkn[1]; 
	      f[atomn][2] += tmp2*rkn[2]; 

	      if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
	    }
	  }
	}
      }
    }

    // Tij forces

    REBO_neighs = REBO_firstneigh[j];
    for (l = 0; l < REBO_numneigh[j]; l++) {
      atoml = REBO_neighs[l];
      if (atoml != atomi) {
	//ltype = map[type[atoml]];
        ltype = species[atoml];

	rjl[0] = x[atomj][0]-x[atoml][0];
	rjl[1] = x[atomj][1]-x[atoml][1];
	rjl[2] = x[atomj][2]-x[atoml][2];
	rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
	wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
	Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	  (wjl*kronecker(jtype,1));
	SpN = Sp(Nlj,Nmin,Nmax,dNlj);

	tmp2 = VA*dN3[1]*dwjl*Etmp/rjlmag;
	f[atomj][0] -= tmp2*rjl[0]; 
	f[atomj][1] -= tmp2*rjl[1]; 
	f[atomj][2] -= tmp2*rjl[2]; 
	f[atoml][0] += tmp2*rjl[0]; 
	f[atoml][1] += tmp2*rjl[1]; 
	f[atoml][2] += tmp2*rjl[2]; 

	if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;
	f[atomj][0] -= tmp2*rjl[0]; 
	f[atomj][1] -= tmp2*rjl[1]; 
	f[atomj][2] -= tmp2*rjl[2]; 
	f[atoml][0] += tmp2*rjl[0]; 
	f[atoml][1] += tmp2*rjl[1]; 
	f[atoml][2] += tmp2*rjl[2]; 

	if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	if (fabs(dNlj) > TOL) {
	  REBO_neighs_l = REBO_firstneigh[atoml];
	  for (n = 0; n < REBO_numneigh[atoml]; n++) {
	    atomn = REBO_neighs_l[n];
	    //ntype = map[type[atomn]];
	    ntype = species[atomn];

	    if (atomn !=atomj) {
	      rln[0] = x[atoml][0]-x[atomn][0];
	      rln[1] = x[atoml][1]-x[atomn][1];
	      rln[2] = x[atoml][2]-x[atomn][2];
	      rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
	      Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

	      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)*Etmp/rlnmag;
	      f[atoml][0] -= tmp2*rln[0];
	      f[atoml][1] -= tmp2*rln[1];
	      f[atoml][2] -= tmp2*rln[2];
	      f[atomn][0] += tmp2*rln[0];
	      f[atomn][1] += tmp2*rln[1];
	      f[atomn][2] += tmp2*rln[2];

	      if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
	    }
	  }
	}
      }
    }
  }

  bij = (0.5*(pij+pji))+piRC+(Tij*Etmp);

  return bij;
}

/* ----------------------------------------------------------------------
   G spline
------------------------------------------------------------------------- */

double REBOFrame::gSpline(double costh, double Nij, int typei,
			   double *dgdc, double *dgdN)
{
  double coeffs[6],dS,g1,g2,dg1,dg2,cut,g;
  int i,j;

  i = 0;
  j = 0;
  g = 0.0;
  cut = 0.0;
  dS = 0.0;
  dg1 = 0.0;
  dg2 = 0.0;
  *dgdc = 0.0;
  *dgdN = 0.0;

  // central atom is Carbon

  if (typei == 0) {
    if (costh < gCdom[0]) costh = gCdom[0];
    if (costh > gCdom[4]) costh = gCdom[4];
    if (Nij >= NCmax) {
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j];
	}
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      g = g2;
      *dgdc = dg2;
      *dgdN = 0.0;
    }
    if (Nij <= NCmin) {
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
	}
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      g = g1;
      *dgdc = dg1;
      *dgdN = 0.0;
    }
    if (Nij > NCmin && Nij < NCmax) {
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
	}
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j];
	}
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      cut = Sp(Nij,NCmin,NCmax,dS);
      g = g2+cut*(g1-g2);
      *dgdc = dg2+(cut*(dg1-dg2));
      *dgdN = dS*(g1-g2);
    }
  }
  
  // central atom is Hydrogen

  if (typei == 1) {
    if (costh < gHdom[0]) costh = gHdom[0];
    if (costh > gHdom[3]) costh = gHdom[3];
    for (i = 0; i < 3; i++) {
      if (costh >= gHdom[i] && costh <= gHdom[i+1]) {
	for (j = 0; j < 6; j++) coeffs[j] = gH[i][j];
      }
    }
    g = Sp5th(costh,coeffs,&dg1);
    *dgdN = 0.0;
    *dgdc = dg1;
  }

  return g;
}

/* ----------------------------------------------------------------------
   Pij spline
------------------------------------------------------------------------- */

double REBOFrame::PijSpline(double NijC, double NijH, int typei, int typej,
			     double dN2[2])
{
  int x,y,i,done;
  double Pij,coeffs[16];
  
  for (i = 0; i < 16; i++) coeffs[i]=0.0;

  x = 0;
  y = 0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  done = 0;

  // if inputs are out of bounds set them back to a point in bounds

  if (typei == 0 && typej == 0) {
    if (NijC < pCCdom[0][0]) NijC=pCCdom[0][0];
    if (NijC > pCCdom[0][1]) NijC=pCCdom[0][1];
    if (NijH < pCCdom[1][0]) NijH=pCCdom[1][0];
    if (NijH > pCCdom[1][1]) NijH=pCCdom[1][1];

    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCCf[(int) NijC][(int) NijH];
      dN2[0] = PCCdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCCdfdy[(int) NijC][(int) NijH];
      done = 1;
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCC[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
    }
  }

  // if inputs are out of bounds set them back to a point in bounds
  
   if (typei == 0 && typej == 1){
     if (NijC < pCHdom[0][0]) NijC=pCHdom[0][0];
     if (NijC > pCHdom[0][1]) NijC=pCHdom[0][1];
      if (NijH < pCHdom[1][0]) NijH=pCHdom[1][0];
      if (NijH > pCHdom[1][1]) NijH=pCHdom[1][1];

    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCHf[(int) NijC][(int) NijH];
      dN2[0] = PCHdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCHdfdy[(int) NijC][(int) NijH];
      done = 1;
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCH[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
    }
  }
  
  if (typei == 1 && typej == 0) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }


  if (typei == 1 && typej == 1) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }
  return Pij;
}

/* ----------------------------------------------------------------------
   PiRC spline
------------------------------------------------------------------------- */

double REBOFrame::piRCSpline(double Nij, double Nji, double Nijconj,
			      int typei, int typej, double dN3[3])
{
  int x,y,z,i,done;
  double piRC,coeffs[64];
  x=0;
  y=0;
  z=0;
  i=0;

  done=0;

  for (i=0; i<64; i++) coeffs[i]=0.0;

  if (typei==0 && typej==0) {
    //if the inputs are out of bounds set them back to a point in bounds
    if (Nij<piCCdom[0][0]) Nij=piCCdom[0][0];
    if (Nij>piCCdom[0][1]) Nij=piCCdom[0][1];
    if (Nji<piCCdom[1][0]) Nji=piCCdom[1][0];
    if (Nji>piCCdom[1][1]) Nji=piCCdom[1][1];
    if (Nijconj<piCCdom[2][0]) Nijconj=piCCdom[2][0];
    if (Nijconj>piCCdom[2][1]) Nijconj=piCCdom[2][1];

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
	fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCCf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCCdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCCdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCCdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCCdom[0][1]; i++)
	if (Nij>=(double) i && Nij<=(double) i+1 || Nij==(double) i) x=i;
      for (i=0; i<piCCdom[1][1]; i++)
	if (Nji>=(double) i && Nji<=(double) i+1 || Nji==(double) i) y=i;
      for (i=0; i<piCCdom[2][1]; i++)
	if (Nijconj>=(double) i && Nijconj<=(double) i+1 || 
	    Nijconj==(double) i) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCC[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }
  
  
  // CH interaction

  if (typei==0 && typej==1 || typei==1 && typej==0) {
    // if the inputs are out of bounds set them back to a point in bounds

    if (Nij<piCHdom[0][0] || Nij>piCHdom[0][1] || 
	Nji<piCHdom[1][0] || Nji>piCHdom[1][1] || 
	Nijconj<piCHdom[2][0] || Nijconj>piCHdom[2][1]) {
      if (Nij<piCHdom[0][0]) Nij=piCHdom[0][0];
      if (Nij>piCHdom[0][1]) Nij=piCHdom[0][1];
      if (Nji<piCHdom[1][0]) Nji=piCHdom[1][0];
      if (Nji>piCHdom[1][1]) Nji=piCHdom[1][1];
      if (Nijconj<piCHdom[2][0]) Nijconj=piCHdom[2][0];
      if (Nijconj>piCHdom[2][1]) Nijconj=piCHdom[2][1];
    }

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
	fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCHdom[0][1]; i++)
	if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piCHdom[1][1]; i++)
	if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piCHdom[2][1]; i++)
	if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }
  
  if (typei==1 && typej==1) {
    if (Nij<piHHdom[0][0] || Nij>piHHdom[0][1] || 
	Nji<piHHdom[1][0] || Nji>piHHdom[1][1] || 
	Nijconj<piHHdom[2][0] || Nijconj>piHHdom[2][1]) {
      Nij=0.0;
      Nji=0.0;
      Nijconj=0.0;
    }
    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
	fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piHHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piHHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piHHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piHHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    } 
    if (done==0) {
      for (i=0; i<piHHdom[0][1]; i++)
	if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piHHdom[1][1]; i++)
	if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piHHdom[2][1]; i++)
	if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piHH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }

  return piRC;
}

/* ----------------------------------------------------------------------
   Tij spline
------------------------------------------------------------------------- */

double REBOFrame::TijSpline(double Nij, double Nji,
			     double Nijconj, double dN3[3])
{
  int x,y,z,i,done;
  double Tijf,coeffs[64];

  x=0;
  y=0;
  z=0;
  i=0;
  Tijf=0.0;
  done=0;
  for (i=0; i<64; i++) coeffs[i]=0.0;

  //if the inputs are out of bounds set them back to a point in bounds

  if (Nij<Tijdom[0][0]) Nij=Tijdom[0][0];
  if (Nij>Tijdom[0][1]) Nij=Tijdom[0][1];
  if (Nji<Tijdom[1][0]) Nji=Tijdom[1][0];
  if (Nji>Tijdom[1][1]) Nji=Tijdom[1][1];
  if (Nijconj<Tijdom[2][0]) Nijconj=Tijdom[2][0];
  if (Nijconj>Tijdom[2][1]) Nijconj=Tijdom[2][1];

  if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
      fabs(Nijconj-floor(Nijconj))<TOL) {
    Tijf=Tf[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[0]=Tdfdx[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[1]=Tdfdy[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[2]=Tdfdz[(int) Nij][(int) Nji][(int) Nijconj];
    done=1;
  }

  if (done==0) {
    for (i=0; i<Tijdom[0][1]; i++) 
      if (Nij>=i && Nij<=i+1) x=i;
    for (i=0; i<Tijdom[1][1]; i++)
      if (Nji>=i && Nji<=i+1) y=i;
    for (i=0; i<Tijdom[2][1]; i++) 
      if (Nijconj>=i && Nijconj<=i+1) z=i;

    for (i=0; i<64; i++) coeffs[i]=Tijc[x][y][z][i];
    Tijf=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
  }

  return Tijf;
}

/* ----------------------------------------------------------------------
   Kronecker delta function
------------------------------------------------------------------------- */

double REBOFrame::kronecker(int a, int b)
{
  double kd;
  if (a == b) kd = 1.0;
  else kd = 0.0;
  return kd;
}

// ----------------------------------------------------------------------
// generic Spline functions
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   fifth order spline evaluation
------------------------------------------------------------------------- */

double REBOFrame::Sp5th(double x, double coeffs[6], double *df)
{
  double f, d;
  const double x2 = x*x;
  const double x3 = x2*x;

  f  = coeffs[0];
  f += coeffs[1]*x;
  d  = coeffs[1];
  f += coeffs[2]*x2;
  d += 2.0*coeffs[2]*x;
  f += coeffs[3]*x3;
  d += 3.0*coeffs[3]*x2;
  f += coeffs[4]*x2*x2;
  d += 4.0*coeffs[4]*x3;
  f += coeffs[5]*x2*x3;
  d += 5.0*coeffs[5]*x2*x2;

  *df = d;
  return f;
}

/* ----------------------------------------------------------------------
   bicubic spline evaluation
------------------------------------------------------------------------- */

double REBOFrame::Spbicubic(double x, double y,
			     double coeffs[16], double df[2])
{
  double f,xn,yn,xn1,yn1,c;
  int i,j;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;

  xn = 1.0;
  for (i = 0; i < 4; i++) {
    yn = 1.0;
    for (j = 0; j < 4; j++) {
      c = coeffs[i*4+j];

      f += c*xn*yn;
      if (i > 0) df[0] += c * ((double) i) * xn1 * yn;
      if (j > 0) df[1] += c * ((double) j) * xn * yn1;

      yn1 = yn;
      yn *= y;
    }
    xn1 = xn;
    xn *= x;
  }

  return f;
}

/* ----------------------------------------------------------------------
   tricubic spline evaluation
------------------------------------------------------------------------- */

double REBOFrame::Sptricubic(double x, double y, double z,
			      double coeffs[64], double df[3])
{
  double f,ir,jr,kr,xn,yn,zn,xn1,yn1,zn1,c;
  int i,j,k;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;
  df[2] = 0.0;

  xn = 1.0;
  for (i = 0; i < 4; i++) {
    ir = (double) i;
    yn = 1.0;
    for (j = 0; j < 4; j++) {
      jr = (double) j;
      zn = 1.0;
      for (k = 0; k < 4; k++) {
        kr = (double) k;
        c = coeffs[16*i+4*j+k];
        f += c*xn*yn*zn;
        if (i > 0) df[0] += c * ir * xn1 * yn * zn;
        if (j > 0) df[1] += c * jr * xn * yn1 * zn;
        if (k > 0) df[2] += c * kr * xn * yn * zn1;
        zn1 = zn;
        zn *= z;
      }
      yn1 = yn;
      yn *= y;
    }
    xn1 = xn;
    xn *= x;
  }

  return f;
}

/* ----------------------------------------------------------------------
   initialize spline knot values
------------------------------------------------------------------------- */

void REBOFrame::spline_init()
{
  int i,j,k;
  INFO("AIREBOFrame::spline_init()");
  
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      PCCf[i][j] = 0.0;
      PCCdfdx[i][j] = 0.0;
      PCCdfdy[i][j] = 0.0;
      PCHf[i][j] = 0.0;
      PCHdfdx[i][j] = 0.0;
      PCHdfdy[i][j] = 0.0;
    }
  }
  
  PCCf[0][2] = -0.00050;
  PCCf[0][3] = 0.0161253646;
  PCCf[1][1] = -0.010960;
  PCCf[1][2] = 0.00632624824;
  PCCf[2][0] = -0.0276030;
  PCCf[2][1] = 0.00317953083;
  
  PCHf[0][1] = 0.209336733;
  PCHf[0][2] = -0.0644496154;
  PCHf[0][3] = -0.303927546;
  PCHf[1][0] = 0.010;
  PCHf[1][1] = -0.125123401;
  PCHf[1][2] = -0.298905246;
  PCHf[2][0] = -0.122042146;
  PCHf[2][1] = -0.300529172;
  PCHf[3][0] = -0.307584705;
  
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      for (k = 0; k < 10; k++) {
	piCCf[i][j][k] = 0.0;
	piCCdfdx[i][j][k] = 0.0;
	piCCdfdy[i][j][k] = 0.0;
	piCCdfdz[i][j][k] = 0.0;
	piCHf[i][j][k] = 0.0;
	piCHdfdx[i][j][k] = 0.0;
	piCHdfdy[i][j][k] = 0.0;
	piCHdfdz[i][j][k] = 0.0;
	piHHf[i][j][k] = 0.0;
	piHHdfdx[i][j][k] = 0.0;
	piHHdfdy[i][j][k] = 0.0;
	piHHdfdz[i][j][k] = 0.0;
	Tf[i][j][k] = 0.0;
	Tdfdx[i][j][k] = 0.0;
	Tdfdy[i][j][k] = 0.0;
	Tdfdz[i][j][k] = 0.0;
      }
    }
  }
  
  for (i = 3; i < 10; i++) piCCf[0][0][i] = 0.0049586079;
  piCCf[1][0][1] = 0.021693495;
  piCCf[0][1][1] = 0.021693495;
  for (i = 2; i < 10; i++) piCCf[1][0][i] = 0.0049586079;
  for (i = 2; i < 10; i++) piCCf[0][1][i] = 0.0049586079;
  piCCf[1][1][1] = 0.05250;
  piCCf[1][1][2] = -0.002088750;
  for (i = 3; i < 10; i++) piCCf[1][1][i] = -0.00804280;
  piCCf[2][0][1] = 0.024698831850;
  piCCf[0][2][1] = 0.024698831850;
  piCCf[2][0][2] = -0.00597133450;
  piCCf[0][2][2] = -0.00597133450;
  for (i = 3; i < 10; i++) piCCf[2][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][2][i] = 0.0049586079;
  piCCf[2][1][1] = 0.00482478490;
  piCCf[1][2][1] = 0.00482478490;
  piCCf[2][1][2] = 0.0150;
  piCCf[1][2][2] = 0.0150;
  piCCf[2][1][3] = -0.010;
  piCCf[1][2][3] = -0.010;
  piCCf[2][1][4] = -0.01168893870;
  piCCf[1][2][4] = -0.01168893870;
  piCCf[2][1][5] = -0.013377877400;
  piCCf[1][2][5] = -0.013377877400;
  piCCf[2][1][6] = -0.015066816000;
  piCCf[1][2][6] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[2][1][i] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[1][2][i] = -0.015066816000;
  piCCf[2][2][1] = 0.0472247850;
  piCCf[2][2][2] = 0.0110;
  piCCf[2][2][3] = 0.0198529350;
  piCCf[2][2][4] = 0.01654411250;
  piCCf[2][2][5] = 0.013235290;
  piCCf[2][2][6] = 0.00992646749999 ;
  piCCf[2][2][7] = 0.006617644999;
  piCCf[2][2][8] = 0.00330882250;
  piCCf[3][0][1] = -0.05989946750;
  piCCf[0][3][1] = -0.05989946750;
  piCCf[3][0][2] = -0.05989946750;
  piCCf[0][3][2] = -0.05989946750;
  for (i = 3; i < 10; i++) piCCf[3][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][3][i] = 0.0049586079;
  piCCf[3][1][2] = -0.0624183760;
  piCCf[1][3][2] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[3][1][i] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[1][3][i] = -0.0624183760;
  piCCf[3][2][1] = -0.02235469150;
  piCCf[2][3][1] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[3][2][i] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[2][3][i] = -0.02235469150;

  piCCdfdx[2][1][1] = -0.026250;
  piCCdfdx[2][1][5] = -0.0271880;
  piCCdfdx[2][1][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdx[2][1][i] = -0.0271880;
  piCCdfdx[1][3][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdx[2][3][i] = 0.031209;

  piCCdfdy[1][2][1] = -0.026250;
  piCCdfdy[1][2][5] = -0.0271880;
  piCCdfdy[1][2][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdy[1][2][i] = -0.0271880;
  piCCdfdy[3][1][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdy[3][2][i] = 0.031209;

  piCCdfdz[1][1][2] = -0.0302715;
  piCCdfdz[2][1][4] = -0.0100220;
  piCCdfdz[1][2][4] = -0.0100220;
  piCCdfdz[2][1][5] = -0.0100220;
  piCCdfdz[1][2][5] = -0.0100220;
  for (i = 4; i < 9; i++) piCCdfdz[2][2][i] = -0.0033090;

  //  make top end of piCC flat instead of zero
  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCCf[i][j][k] = piCCf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){ // also enforces some symmetry
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCCf[i][j][k] = piCCf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCCf[4][4][k] = piCCf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCCf[i][j][k] = piCCf[i][j][k-1];
      }
  }

  piCHf[1][1][1] = -0.050;
  piCHf[1][1][2] = -0.050;
  piCHf[1][1][3] = -0.30;
  for (i = 4; i < 10; i++) piCHf[1][1][i] = -0.050;
  for (i = 5; i < 10; i++) piCHf[2][0][i] = -0.004523893758064;
  for (i = 5; i < 10; i++) piCHf[0][2][i] = -0.004523893758064;
  piCHf[2][1][2] = -0.250;
  piCHf[1][2][2] = -0.250;
  piCHf[2][1][3] = -0.250;
  piCHf[1][2][3] = -0.250;
  piCHf[3][1][1] = -0.10;
  piCHf[1][3][1] = -0.10;
  piCHf[3][1][2] = -0.125;
  piCHf[1][3][2] = -0.125;
  piCHf[3][1][3] = -0.125;
  piCHf[1][3][3] = -0.125;
  for (i = 4; i < 10; i++) piCHf[3][1][i] = -0.10;
  for (i = 4; i < 10; i++) piCHf[1][3][i] = -0.10;
  
  // make top end of piCH flat instead of zero
 // also enforces some symmetry

  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCHf[i][j][k] = piCHf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCHf[i][j][k] = piCHf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCHf[4][4][k] = piCHf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCHf[i][j][k] = piCHf[i][j][k-1];
      }
  }

  piHHf[1][1][1] = 0.124915958;

  Tf[2][2][1] = -0.035140;
  for (i = 2; i < 10; i++) Tf[2][2][i] = -0.0040480;
  
}


// ----------------------------------------------------------------------
// S'(t) and S(t) cutoff functions
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   cutoff function Sprime
   return cutoff and dX = derivative
------------------------------------------------------------------------- */

double REBOFrame::Sp(double Xij, double Xmin, double Xmax, double &dX)
{
  double cutoff;

  double t = (Xij-Xmin) / (Xmax-Xmin);
  if (t <= 0.0) {
    cutoff = 1.0;
    dX = 0.0;
  } 
  else if (t >= 1.0) {
    cutoff = 0.0;
    dX = 0.0;
  } 
  else {
    cutoff = 0.5 * (1.0+cos(M_PI*t));
    dX = (-0.5*M_PI*sin(M_PI*t)) / (Xmax-Xmin);
  }
  return cutoff;
}

/* ----------------------------------------------------------------------
   LJ cutoff function Sp2
   return cutoff and dX = derivative
------------------------------------------------------------------------- */

double REBOFrame::Sp2(double Xij, double Xmin, double Xmax, double &dX)
{
  double cutoff;

  double t = (Xij-Xmin) / (Xmax-Xmin);
  if (t <= 0.0) {
    cutoff = 1.0;
    dX = 0.0;
  } 
  if (t >= 1.0) {
    cutoff = 0.0;
    dX = 0.0;
  } 
  if (t>0.0 && t<1.0) {
    cutoff = (1.0-(t*t*(3.0-2.0*t)));
    dX = 6.0*(t*t-t) / (Xmax-Xmin);
  }
  return cutoff;
}

/* ----------------------------------------------------------------------
   read AIREBO potential file
------------------------------------------------------------------------- */

void REBOFrame::read_file(char *filename)
{
  int i,j,k,l,limit;
  char s[MAXLINE];
  
  // REBO Parameters (AIREBO)
  
  double rcmin_CC,rcmin_CH,rcmin_HH,rcmax_CC,rcmax_CH,
    rcmax_HH,rcmaxp_CC,rcmaxp_CH,rcmaxp_HH;
  double Q_CC,Q_CH,Q_HH,alpha_CC,alpha_CH,alpha_HH,A_CC,A_CH,A_HH;
  double BIJc_CC1,BIJc_CC2,BIJc_CC3,BIJc_CH1,BIJc_CH2,BIJc_CH3,
    BIJc_HH1,BIJc_HH2,BIJc_HH3;
  double Beta_CC1,Beta_CC2,Beta_CC3,Beta_CH1,Beta_CH2,Beta_CH3,
    Beta_HH1,Beta_HH2,Beta_HH3;
  double rho_CC,rho_CH,rho_HH;
  
  // LJ Parameters (AIREBO)
  
  double rcLJmin_CC,rcLJmin_CH,rcLJmin_HH,rcLJmax_CC,rcLJmax_CH,
    rcLJmax_HH,bLJmin_CC;
  double bLJmin_CH,bLJmin_HH,bLJmax_CC,bLJmax_CH,bLJmax_HH,
    epsilon_CC,epsilon_CH,epsilon_HH;
  double sigma_CC,sigma_CH,sigma_HH,epsilonT_CCCC,epsilonT_CCCH,epsilonT_HCCH;
  
    FILE *fp = fopen(filename,"r");
    if (fp == NULL)
        ERROR("Cannot open AIREBO potential file"<<filename);

    // skip initial comment lines

    while (1) {
      fgets(s,MAXLINE,fp);
      if (s[0] != '#') break;
    }
    
    // read parameters
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&smin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_HH);
    
    // LJ parameters
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_HCCH);
    
    // gC spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    // number-1 = # of domains for the spline 	
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gCdom[i]);
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&gC1[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&gC2[i][j]);
      }
    }
    
    // gH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);	
    
    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gHdom[i]);
    }
    
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&gH[i][j]);
      }
    }
    
    // pCC spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);	
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&pCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) pCCdom[0][1]; i++) {
      for (j = 0; j < (int) pCCdom[1][1]; j++) {
	for (k = 0; k < 16; k++) {
	  fgets(s,MAXLINE,fp);
	  sscanf(s,"%lg",&pCC[i][j][k]);
	}
      }
    }
    
    // pCH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);	
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&pCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) pCHdom[0][1]; i++) {
      for (j = 0; j < (int) pCHdom[1][1]; j++) {
	for (k = 0; k < 16; k++) {
	  fgets(s,MAXLINE,fp);
	  sscanf(s,"%lg",&pCH[i][j][k]);
	}
      }
    }
    
    // piCC cpline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&piCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) piCCdom[0][1]; i++) {
      for (j = 0; j < (int) piCCdom[1][1]; j++) {
	for (k = 0; k < (int) piCCdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&piCC[i][j][k][l]);
	  }
	}
      }
    }
    
    // piCH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&piCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) piCHdom[0][1]; i++) {
      for (j = 0; j < (int) piCHdom[1][1]; j++) {
	for (k = 0; k < (int) piCHdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&piCH[i][j][k][l]);
	  }
	}
      }
    }
    
    // piHH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&piHHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) piHHdom[0][1]; i++) {
      for (j = 0; j < (int) piHHdom[1][1]; j++) {
	for (k = 0; k < (int) piHHdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&piHH[i][j][k][l]);
	  }
	}
      }
    }
    
    // Tij spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&Tijdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) Tijdom[0][1]; i++) {
      for (j = 0; j < (int) Tijdom[1][1]; j++) {
	for (k = 0; k < (int) Tijdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&Tijc[i][j][k][l]);
	  }
	}
      }
    }
    
    fclose(fp);

  
    // store read-in values in arrays
    
    // REBO
   
    rcmin[0][0] = rcmin_CC;
    rcmin[0][1] = rcmin_CH;
    rcmin[1][0] = rcmin[0][1];
    rcmin[1][1] = rcmin_HH;
    
    rcmax[0][0] = rcmax_CC;
    rcmax[0][1] = rcmax_CH;
    rcmax[1][0] = rcmax[0][1];
    rcmax[1][1] = rcmax_HH;
    
    rcmaxsq[0][0] = rcmax[0][0]*rcmax[0][0];
    rcmaxsq[1][0] = rcmax[1][0]*rcmax[1][0];
    rcmaxsq[0][1] = rcmax[0][1]*rcmax[0][1];
    rcmaxsq[1][1] = rcmax[1][1]*rcmax[1][1];
    
    rcmaxp[0][0] = rcmaxp_CC;
    rcmaxp[0][1] = rcmaxp_CH;
    rcmaxp[1][0] = rcmaxp[0][1];
    rcmaxp[1][1] = rcmaxp_HH;
    
    Q[0][0] = Q_CC;
    Q[0][1] = Q_CH;
    Q[1][0] = Q[0][1];
    Q[1][1] = Q_HH;
    
    alpha[0][0] = alpha_CC;
    alpha[0][1] = alpha_CH;
    alpha[1][0] = alpha[0][1];
    alpha[1][1] = alpha_HH;
    
    A[0][0] = A_CC;
    A[0][1] = A_CH;
    A[1][0] = A[0][1];
    A[1][1] = A_HH;
    
    rho[0][0] = rho_CC;
    rho[0][1] = rho_CH;
    rho[1][0] = rho[0][1];
    rho[1][1] = rho_HH;
    
    BIJc[0][0][0] = BIJc_CC1;
    BIJc[0][0][1] = BIJc_CC2;
    BIJc[0][0][2] = BIJc_CC3;
    BIJc[0][1][0] = BIJc_CH1;
    BIJc[0][1][1] = BIJc_CH2;
    BIJc[0][1][2] = BIJc_CH3;
    BIJc[1][0][0] = BIJc_CH1;
    BIJc[1][0][1] = BIJc_CH2;
    BIJc[1][0][2] = BIJc_CH3;
    BIJc[1][1][0] = BIJc_HH1;
    BIJc[1][1][1] = BIJc_HH2;
    BIJc[1][1][2] = BIJc_HH3;
    
    Beta[0][0][0] = Beta_CC1;
    Beta[0][0][1] = Beta_CC2;
    Beta[0][0][2] = Beta_CC3;
    Beta[0][1][0] = Beta_CH1;
    Beta[0][1][1] = Beta_CH2;
    Beta[0][1][2] = Beta_CH3;
    Beta[1][0][0] = Beta_CH1;
    Beta[1][0][1] = Beta_CH2;
    Beta[1][0][2] = Beta_CH3;
    Beta[1][1][0] = Beta_HH1;
    Beta[1][1][1] = Beta_HH2;
    Beta[1][1][2] = Beta_HH3;
    
    // LJ
    
    rcLJmin[0][0] = rcLJmin_CC;
    rcLJmin[0][1] = rcLJmin_CH;
    rcLJmin[1][0] = rcLJmin[0][1];
    rcLJmin[1][1] = rcLJmin_HH;
    
    rcLJmax[0][0] = rcLJmax_CC;
    rcLJmax[0][1] = rcLJmax_CH;
    rcLJmax[1][0] = rcLJmax[0][1];
    rcLJmax[1][1] = rcLJmax_HH;
    
    rcLJmaxsq[0][0] = rcLJmax[0][0]*rcLJmax[0][0];
    rcLJmaxsq[1][0] = rcLJmax[1][0]*rcLJmax[1][0];
    rcLJmaxsq[0][1] = rcLJmax[0][1]*rcLJmax[0][1];
    rcLJmaxsq[1][1] = rcLJmax[1][1]*rcLJmax[1][1];

    bLJmin[0][0] = bLJmin_CC;
    bLJmin[0][1] = bLJmin_CH;
    bLJmin[1][0] = bLJmin[0][1];
    bLJmin[1][1] = bLJmin_HH;
    
    bLJmax[0][0] = bLJmax_CC;
    bLJmax[0][1] = bLJmax_CH;
    bLJmax[1][0] = bLJmax[0][1];
    bLJmax[1][1] = bLJmax_HH;
    
    epsilon[0][0] = epsilon_CC;
    epsilon[0][1] = epsilon_CH;
    epsilon[1][0] = epsilon[0][1];
    epsilon[1][1] = epsilon_HH;
    
    sigma[0][0] = sigma_CC;
    sigma[0][1] = sigma_CH;
    sigma[1][0] = sigma[0][1];
    sigma[1][1] = sigma_HH;
    
    // torsional

    thmin = -1.0;
    thmax = -0.995; 
    epsilonT[0][0] = epsilonT_CCCC;
    epsilonT[0][1] = epsilonT_CCCH;
    epsilonT[1][0] = epsilonT[0][1];
    epsilonT[1][1] = epsilonT_HCCH; 
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void REBOFrame::ev_tally(int i, int j, 
		    double evdwl, double ecoul, double fpair,
		    double delx, double dely, double delz)
{
  double evdwlhalf,ecoulhalf,epairhalf;
  Matrix33 vv;
  Vector3  dr; 

  eng_vdwl += evdwl;
  eng_coul += ecoul;

  evdwlhalf = 0.5*evdwl;
  ecoulhalf = 0.5*ecoul;
  epairhalf = 0.5 * (evdwl + ecoul);

  _EPOT_IND[i] += epairhalf;
  _EPOT_IND[j] += epairhalf;
  _EPOT += epairhalf * 2.0;

  dr.set(delx, dely, delz);
  vv.clear();
  vv.addnvv(fpair*0.5,dr,dr);

  _VIRIAL_IND[i] += vv;
  _VIRIAL_IND[j] += vv;
  _VIRIAL += vv * 2.0;
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
   fpair is magnitude of force on atom I
------------------------------------------------------------------------- */

void REBOFrame::v_tally2(int i, int j, double fpair, double *drij)
{
  Matrix33 vv;
  Vector3 dr;

  dr.set(drij[0],drij[1],drij[2]);
  vv.clear();
  vv.addnvv(fpair*0.5, dr, dr);  

  _VIRIAL_IND[i] += vv;
  _VIRIAL_IND[j] += vv;
  _VIRIAL += vv * 2.0;  
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO and Tersoff potential, newton_pair is always on
------------------------------------------------------------------------- */

void REBOFrame::v_tally3(int i, int j, int k,
		    double *fi, double *fj, double *drik, double *drjk)
{
  Matrix33 vv;
  Vector3 dr, df;

  vv.clear();
  dr.set(drik); df.set(fi); vv.addnvv(1.0/3.0, dr, df);
  dr.set(drjk); df.set(fj); vv.addnvv(1.0/3.0, dr, df);

  _VIRIAL_IND[i] += vv;
  _VIRIAL_IND[j] += vv;
  _VIRIAL_IND[k] += vv;
  _VIRIAL += vv * 3.0;
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
------------------------------------------------------------------------- */

void REBOFrame::v_tally4(int i, int j, int k, int m,
		    double *fi, double *fj, double *fk,
		    double *drim, double *drjm, double *drkm)
{
  Matrix33 vv;
  Vector3 dr, df;

  dr.set(drim); df.set(fi); vv.addnvv(0.25, dr, df);
  dr.set(drjm); df.set(fj); vv.addnvv(0.25, dr, df);
  dr.set(drkm); df.set(fk); vv.addnvv(0.25, dr, df);

  _VIRIAL_IND[i] += vv;
  _VIRIAL_IND[j] += vv;
  _VIRIAL_IND[k] += vv;
  _VIRIAL_IND[m] += vv;
  _VIRIAL += vv * 4.0;
}

void REBOFrame::potential()
{
  int ipt;

  _EPOT=0; _VIRIAL.clear();
  for(ipt=0;ipt<_NP;ipt++) { _F[ipt].clear(); _EPOT_IND[ipt]=0; }

  rebo();
}

int REBOFrame::readAIREBO()
{

    LFile::SubHomeDir(rebofile,rebofile);
    read_file(rebofile);

    spline_init();

    /* Fitting parameters for C */  
    //acut=1.0;
    _RLIST=acut*1.1;
    _SKIN=_RLIST-acut;
    acutsq=acut*acut;

    return 0;
}

#ifdef _TEST

/* Main Program Begins */
class REBOFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

