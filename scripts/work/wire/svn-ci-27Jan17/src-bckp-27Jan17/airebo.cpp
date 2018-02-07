/*
  rebo.cpp
  by Keonwook Kang  kwkang75@yonsei.ac.kr
  Last Modified : Thurs Feb 04 11:28:54 2014

  FUNCTION  :  MD simulation package of C using REBO potential
  Reference:   
  (1) D. W. Brenner et al, J. Phys.: Condens. Matter 14 (2002) 783

  Note:        Implemented based on LAMMPS (02/04/14)

  Atom species: (need to check with LAMMPS manual)
               0 - C
               1 - H
  To do:
        1. translate MD++ neighbor list to LAMMPS neighbor list (DONE)
        2. convert LAMMPS tally_# functions to MD++ virial/force accumulation (DONE)
        3. construct a function to read a potential file and parameters. (DONE)
           ex.rcmin,rcmax,rcmaxsq 
        4. Free box relaxation does not fully converge.
               potential      gradient
           (1) pure REBO       1e-5
           (2) REBO+LJ         1e-1
           (3) REBO+Torsion    1e-5
           (4) REBO+LJ+Torsion 1e0
           LAMMPS has a similar issue. 
        5. In AIREBOFrame::flj(), atomk and atomm are defined 
           only if certain conditions are satisfied. Is it OK?
           
*/

#include "airebo.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024
#define TOL 1.0e-9
#define PGDELTA 1

void AIREBOFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int size;
    size=_NP*allocmultiple;

    Realloc(REBO_R,Vector3,size);
    Realloc(REBO_numneigh,int,size);
    Realloc(closestdistsq,double,size);
    Realloc(nC,double,size);
    Realloc(nH,double,size);

}

void AIREBOFrame::initvars()
{   
    MDPARALLELFrame::initvars();

    strcpy(incnfile,"../carbon.cn");
    strcpy(rebofile,"rebof");

    REBO_numneigh = NULL;
    REBO_firstneigh = NULL;
    REBO_firstneigh_mem = NULL;
    closestdistsq = NULL;
    nC = nH = NULL;
    
    DUMP(HIG"AIREBOFrame initvars"NOR);
}

void AIREBOFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("rebofile",rebofile,STRING);
    //bindvar("rcut",&acut,DOUBLE);
    bindvar("ljcut",&cutlj,DOUBLE);
    bindvar("ljflag",&ljflag,INT);
    bindvar("torsionflag",&torsionflag,INT);
}

int AIREBOFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readREBO",readAIREBO());
    bindcommand(name,"readAIREBO",readAIREBO());

    return -1;
}

void AIREBOFrame::airebo()
{ 
    DUMP("AIREBO");
    
    refreshneighborlist();
    REBO_neigh();

    rebo();
    if(ljflag) flj();
    if(torsionflag) torsion();

    //if (!conj_fixbox) virial_fdotr_compute();
    //virial_fdotr_compute();
    
    /* zero out the forces on fixed atoms */
    //  for(ipt=0;ipt<_NP;ipt++)
    //  {
    //     if(fixed[ipt]) _F[ipt].clear();
    //  }
}

void AIREBOFrame::rebo()
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
  Vector3 si,sj,sij,rij,fij;

  n0=0;
  n1=_NP;

  evdwl = 0.0;
  ecoul = 0.0;

  SHtoR();
  for(ipt=n0;ipt<n1;ipt++) REBO_R[ipt] = _R[ipt];  

  // two-body interactions from REBO neighbor list, skip half of them
  for(ipt=n0;ipt<n1;ipt++)
  {
    if(fixed[ipt]==-1) continue; /* ignore atom if -1 */
    itype = species[ipt];
    /* move all neighboring atoms to the nearest image around atom ipt */
    REBO_R[ipt] = _H*_SR[ipt];
    for(j=0;j<REBO_numneigh[ipt];j++)
    {
      jpt=REBO_firstneigh[ipt][j];
      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
      sij.subtract(_SR[ipt],_SR[jpt]);
      sij.subint();
      rij=_H*sij;
      REBO_R[jpt]=REBO_R[ipt]-rij;
    }
    
    for(j=0;j<REBO_numneigh[ipt];j++)
    {
      jpt=REBO_firstneigh[ipt][j];
      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
    
      /* skip if ipt >= jpt */
      if (ipt>=jpt) continue;

      jtype = species[jpt];
      rij=REBO_R[ipt]-REBO_R[jpt];  
      r2ij=rij.norm2();
      //if(r2ij>acutsq) continue;
      rrij=sqrt(r2ij);

      /* cut-off fcn wij for pair interaction */
      wij = Sp(rrij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
      if (wij <= TOL) continue;
      
      Qij = Q[itype][jtype];
      Aij = A[itype][jtype];
      alphaij = alpha[itype][jtype];

      /* Repulsive pair potential VR, eqn (5) in Brenner et al (2002) */
      VR = wij*(1.0+(Qij/rrij)) * Aij*exp(-alphaij*rrij);
      pre = wij*Aij * exp(-alphaij*rrij);
      dVRdi = pre * ((-alphaij)-(Qij/r2ij)-(Qij*alphaij/rrij));
      dVRdi += VR/wij * dwij;

      /* Attractive pair potential VA, eqn (6) in Brenner et al (2002) */
      VA = dVA = 0.0;
      for (m = 0; m < 3; m++) {
        term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rrij);
        VA += term;
        dVA += -Beta[itype][jtype][m] * term;
      }
      dVA += VA/wij * dwij;

      /* Bond order term bij */
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
}

/* ----------------------------------------------------------------------
   compute LJ forces and energy
   find 3- and 4-step paths between atoms I,J via REBO neighbor lists
------------------------------------------------------------------------- */
void AIREBOFrame::flj()
{
  int j,k,m,ipt,jpt,kpt,mpt,itype,jtype,ktype,mtype;
  int atomk,atomm,testpath,npath,done;
  double evdwl,fpair;
  double best,wik,wkm,cij,dwij,dwik,dwkj,dwkm,dwmj;
  double rrij,r2ij,rrik,rrjk,rsq;
  double wkj,dC,VLJ,dVLJ,VA,Str,dStr,Stb;
  double vdw,slw,dvdw,dslw,drij,swidth,tee,tee2;
  double rljmin,rljmax,sigcut,sigmin,sigwid;
  double rrkm,rrmj,wmj,r2inv,r6inv,scale;
  double rikS,rkjS,rkmS,rmjS,wikS,dwikS;
  double wkjS,dwkjS,wkmS,dwkmS,wmjS,dwmjS;
  double fpair1,fpair2,fpair3;
  Vector3 fi,fj,fk,fm,delscale,fij;
  Vector3 sij,rij,rik,rjk,rkm,rjm,delikS,deljkS,delkmS,deljmS,delimS;
  int vflag_atom=1;
  
  /* MD++ variables */
  int n0=0,n1=_NP;
  
  // I-J interaction from full neighbor list
  // skip 1/2 of interactions since only consider each pair once

  evdwl = 0.0;
  rljmin = 0.0;
  rljmax = 0.0;
  sigcut = 0.0;
  sigmin = 0.0;
  sigwid = 0.0;

  SHtoR();
  for(ipt=n0;ipt<n1;ipt++) REBO_R[ipt] = _R[ipt];
  
  for (ipt = n0; ipt < n1; ipt++) {
    if(fixed[ipt]==-1) continue;
    itype = species[ipt];
    
    /* move all neighboring atoms to the nearest image around atom ipt */
    REBO_R[ipt] = _H*_SR[ipt];
    for(j=0;j<nn[ipt];j++)
    {
      jpt=nindex[ipt][j];
      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
      sij.subtract(_SR[ipt],_SR[jpt]);
      sij.subint();
      rij=_H*sij;
      REBO_R[jpt]=REBO_R[ipt]-rij;
    }    

    for (j = 0; j < nn[ipt]; j++) {
      jpt=nindex[ipt][j];
      if(fixed[jpt]==-1) continue;

      /* skip if ipt >= jpt */
      if (ipt>=jpt) continue;
      
      jtype = species[jpt];
      rij=REBO_R[ipt]-REBO_R[jpt];
      r2ij=rij.norm2();

      // if outside of LJ cutoff, skip
      // if outside of 4-path cutoff, best = 0.0, no need to test paths
      // if outside of 2-path cutoff but inside 4-path cutoff,
      //   best = 0.0, test 3-,4-paths
      // if inside 2-path cutoff, best = wij, only test 3-,4-paths if best < 1

      if (r2ij >= cutljsq[itype][jtype]) continue;
      rrij = sqrt(r2ij);
      if (rrij >= cut3rebo) {
        best = 0.0;
        testpath = 0;
      } else if (rrij >= rcmax[itype][jtype]) {
        best = 0.0;
        testpath = 1;
      } else {
        best = Sp(rrij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
        npath = 2;
        if (best < 1.0) testpath = 1;
        else testpath = 0;
      }

      done = 0;
      if (testpath) {

        // test all 3-body paths = I-K-J
        // I-K interactions come from atom I's REBO neighbors
        // if wik > current best, compute wkj
        // if best = 1.0, done

        for (k = 0; k < REBO_numneigh[ipt] && done==0; k++) {
          kpt = REBO_firstneigh[ipt][k];
          if (kpt == jpt) continue;
          ktype = species[kpt];

          rik=REBO_R[ipt]-REBO_R[kpt];
          rsq=rik.norm2();

          if (rsq < rcmaxsq[itype][ktype]) {
            rrik = sqrt(rsq);
            wik = Sp(rrik,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
          } else wik = 0.0;

          if (wik > best) {
            rjk=REBO_R[jpt]-REBO_R[kpt];
            rsq=rjk.norm2();

            if (rsq < rcmaxsq[ktype][jtype]) {
              rrjk = sqrt(rsq);
              wkj = Sp(rrjk,rcmin[ktype][jtype],rcmax[ktype][jtype],dwkj);
              if (wik*wkj > best) {
                best = wik*wkj;
                npath = 3;
                atomk = kpt;
                delikS = rik;
                rikS = rrik;
                wikS = wik;
                dwikS = dwik;
                deljkS = rjk;
                rkjS = rrjk;
                wkjS = wkj;
                dwkjS = dwkj;
                if (best == 1.0) {
                  done = 1;
                  break;
                }
              }
            }

            // test all 4-body paths = I-K-M-J
            // K-M interactions come from atom K's REBO neighbors
            // if wik*wkm > current best, compute wmj
            // if best = 1.0, done

            for (m = 0; m < REBO_numneigh[kpt] && done==0; m++) {
              mpt = REBO_firstneigh[kpt][m];
              if (mpt == ipt || mpt == jpt) continue;
              mtype = species[mpt];
              rkm = REBO_R[kpt]-REBO_R[mpt];
              rsq=rkm.norm2(); 

              if (rsq < rcmaxsq[ktype][mtype]) {
                rrkm = sqrt(rsq);
                wkm = Sp(rrkm,rcmin[ktype][mtype],rcmax[ktype][mtype],dwkm);
              } else wkm = 0.0;

              if (wik*wkm > best) {
                rjm = REBO_R[jpt]-REBO_R[mpt];
                rsq = rjm.norm2();

                if (rsq < rcmaxsq[mtype][jtype]) {
                  rrmj = sqrt(rsq);
                  wmj = Sp(rrmj,rcmin[mtype][jtype],rcmax[mtype][jtype],dwmj);
                  if (wik*wkm*wmj > best) {
                    best = wik*wkm*wmj;
                    npath = 4;
                    atomk = kpt;
                    delikS = rik;
                    rikS = rrik;
                    wikS = wik;
                    dwikS = dwik;
                    atomm = mpt;
                    delkmS = rkm;
                    rkmS = rrkm;
                    wkmS = wkm;
                    dwkmS = dwkm;
                    deljmS = rjm;
                    rmjS = rrmj;
                    wmjS = wmj;
                    dwmjS = dwmj;
                    if (best == 1.0) {
                      done = 1;
                      break;
                    }
                  }
                }
              }
            } // for (m = 0; ...
          } //  if (wik > best) 
        } // for (k = 0; ...)
      } // if (testpath) 

      cij = 1.0 - best;
      if (cij == 0.0) continue;

      // compute LJ forces and energy

      sigwid = 0.84;
      sigcut = 3.0;
      sigmin = sigcut - sigwid;

      rljmin = sigma[itype][jtype];
      rljmax = sigcut * rljmin;
      rljmin = sigmin * rljmin;

      if (rrij > rljmax) {
        slw = 0.0;
        dslw = 0.0;
      } else if (rrij > rljmin) {
        drij = rrij - rljmin;
        swidth = rljmax - rljmin;
        tee = drij / swidth;
        tee2 = tee*tee;
        slw = 1.0 - tee2 * (3.0 - 2.0 * tee);
        dslw = 6.0 * tee * (1.0 - tee) / rrij / swidth;
      } else {
        slw = 1.0;
        dslw = 0.0;
      }

      r2inv = 1.0/r2ij;
      r6inv = r2inv*r2inv*r2inv;

      vdw = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
      dvdw = -r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]) / rrij;

      // VLJ now becomes vdw * slw, derivaties, etc.

      VLJ = vdw * slw;
      dVLJ = dvdw * slw + vdw * dslw;

      Str = Sp2(rrij,rcLJmin[itype][jtype],rcLJmax[itype][jtype],dStr);
      VA = Str*cij*VLJ;
      if (Str > 0.0) {
        scale = rcmin[itype][jtype] / rrij;
        delscale = rij*scale;
        Stb = bondorderLJ(ipt,jpt,&delscale[0],rcmin[itype][jtype],VA,
                          &rij[0],rrij,_F,vflag_atom);
      } else Stb = 0.0;

      fpair = -(dStr * (Stb*cij*VLJ - cij*VLJ) +
                dVLJ * (Str*Stb*cij + cij - Str*cij)) / rrij;

      fij.clear(); fij.addnv(fpair,rij);
      _F[ipt] += fij;
      _F[jpt] -= fij;

      if (eflag) evdwl = VA*Stb + (1.0-Str)*cij*VLJ;
      if (evflag) ev_tally(ipt,jpt,evdwl,0.0,fpair,rij.x,rij.y,rij.z);
      
      if (cij < 1.0) {
        dC = Str*Stb*VLJ + (1.0-Str)*VLJ;
        if (npath == 2) {
          fpair = dC*dwij / rrij;
          fij.clear(); fij.addnv(fpair,rij);
          _F[ipt] += fij;
          _F[jpt] -= fij;

          if (vflag_atom) v_tally2(ipt,jpt,fpair,&rij[0]);

        } else if (npath == 3) {
          fpair1 = dC*dwikS*wkjS / rikS;
          fi = delikS*fpair1;
          fpair2 = dC*wikS*dwkjS / rkjS;
          fj = deljkS*fpair2;

          _F[ipt] += fi;
          _F[jpt] += fj;
          //_F[kpt] -= (fi+fj);
          _F[atomk] -= (fi+fj);

          //if (vflag_atom)
          //  v_tally3(ipt,jpt,kpt,&fi[0],&fj[0],&delikS[0],&deljkS[0]);
          if (vflag_atom)
            v_tally3(ipt,jpt,atomk,&fi[0],&fj[0],&delikS[0],&deljkS[0]);

        } else {
          fpair1 = dC*dwikS*wkmS*wmjS / rikS;
          fi = delikS*fpair1;

          fpair2 = dC*wikS*dwkmS*wmjS / rkmS;
          fk = delkmS*fpair2; fk -= fi;

          fpair3 = dC*wikS*wkmS*dwmjS / rmjS;
          fj = deljmS*fpair3;

          fm = delkmS*((-1.0)*fpair2); fm -= fj;

          _F[ipt] += fi;
          _F[jpt] += fj;
          //_F[kpt] += fk;
          //_F[mpt] += fm;
          _F[atomk] += fk;
          _F[atomm] += fm;

          if (vflag_atom) {
            delimS = (delikS + delkmS);
            //v_tally4(ipt,jpt,kpt,mpt,&fi[0],&fj[0],&fk[0],&delimS[0],&deljmS[0],&delkmS[0]);
            v_tally4(ipt,jpt,atomk,atomm,&fi[0],&fj[0],&fk[0],&delimS[0],&deljmS[0],&delkmS[0]);
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   torsional forces and energy
------------------------------------------------------------------------- */

void AIREBOFrame::torsion()
{
  int j,k,l;
  double evdwl,fpair;
  double cos321;
  double w21,dw21,cos234,w34,dw34;
  double cross321[3],cross321mag,cross234[3],cross234mag;
  double w23,dw23,cw2,ekijl,Ec;
  double cw,cwnum,cwnom;
  double rij,rij2,rik,rjl,tspjik,dtsjik,tspijl,dtsijl,costmp,fcpc;
  double sin321,sin234,rjk2,rik2,ril2,rjl2;
  double rjk,ril;
  double Vtors;
  double dndij[3],tmpvec[3],dndik[3],dndjl[3];
  double dcidij,dcidik,dcidjk,dcjdji,dcjdjl,dcjdil;
  double dsidij,dsidik,dsidjk,dsjdji,dsjdjl,dsjdil;
  double dxidij,dxidik,dxidjk,dxjdji,dxjdjl,dxjdil;
  double ddndij,ddndik,ddndjk,ddndjl,ddndil,dcwddn,dcwdn,dvpdcw,Ftmp[3];
  double del32[3],rsq,r32,del23[3],del21[3],r21;
  double deljk[3],del34[3],delil[3],delkl[3],r23,r34;
  double fi[3],fj[3],fk[3],fl[3];
  int itype,jtype,ktype,ltype;
  int *REBO_neighs_i,*REBO_neighs_j;
  class Vector3 *x;

  /* MD++ variables */
  int n0=0,n1=_NP;
  int ipt,jpt,kpt,lpt;
  Vector3 sij,Rij;
  
  //double **f = atom->f;
  //int *type = atom->type;
  //int *tag = atom->tag;

  SHtoR();
  for(ipt=n0;ipt<n1;ipt++) REBO_R[ipt] = _R[ipt];  

  for (ipt = n0; ipt < n1; ipt++) {
    if(fixed[ipt]==-1) continue;
    
    /* move all neighboring atoms to the nearest image around atom ipt */
    REBO_R[ipt] = _H*_SR[ipt];
    for(j=0;j<nn[ipt];j++)
    {
      jpt=nindex[ipt][j];
      if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
      sij.subtract(_SR[ipt],_SR[jpt]);
      sij.subint();
      Rij=_H*sij;
      REBO_R[jpt]=REBO_R[ipt]-Rij;
    }
    
    itype = species[ipt];
    if (itype != 0) continue; /* skip H atom */

    x = REBO_R;

    for (j = 0; j < REBO_numneigh[ipt]; j++) {
      jpt = REBO_firstneigh[ipt][j];

      if(ipt>=jpt) continue;
      
      jtype = species[jpt];
      if (jtype != 0) continue; /* skip H atom */

      del32[0] = x[jpt][0]-x[ipt][0];
      del32[1] = x[jpt][1]-x[ipt][1];
      del32[2] = x[jpt][2]-x[ipt][2];
      rsq = del32[0]*del32[0] + del32[1]*del32[1] + del32[2]*del32[2];
      r32 = sqrt(rsq);
      del23[0] = -del32[0];
      del23[1] = -del32[1];
      del23[2] = -del32[2];
      r23 = r32;
      w23 = Sp(r23,rcmin[itype][jtype],rcmax[itype][jtype],dw23);

      for (k = 0; k < REBO_numneigh[ipt]; k++) {
        kpt = REBO_firstneigh[ipt][k];
        ktype = species[kpt];
        if (kpt == jpt) continue;
        del21[0] = x[ipt][0]-x[kpt][0];
        del21[1] = x[ipt][1]-x[kpt][1];
        del21[2] = x[ipt][2]-x[kpt][2];
        rsq = del21[0]*del21[0] + del21[1]*del21[1] + del21[2]*del21[2];
        r21 = sqrt(rsq);
        cos321 = - ((del21[0]*del32[0]) + (del21[1]*del32[1]) +
                    (del21[2]*del32[2])) / (r21*r32);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        sin321 = sqrt(1.0 - cos321*cos321);
        if (sin321 < TOL) continue;

        deljk[0] = del21[0]-del23[0];
        deljk[1] = del21[1]-del23[1];
        deljk[2] = del21[2]-del23[2];
        rjk2 = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
        rjk=sqrt(rjk2);
        rik2 = r21*r21;
        w21 = Sp(r21,rcmin[itype][ktype],rcmax[itype][ktype],dw21);

        rij = r32;
        rik = r21;
        rij2 = r32*r32;
        rik2 = r21*r21;
        costmp = 0.5*(rij2+rik2-rjk2)/rij/rik;
        tspjik = Sp2(costmp,thmin,thmax,dtsjik);
        dtsjik = -dtsjik;

        REBO_neighs_j = REBO_firstneigh[jpt];
        for (l = 0; l < REBO_numneigh[jpt]; l++) {
          lpt = REBO_neighs_j[l];
          ltype = species[lpt];
          if (lpt == ipt || lpt == kpt) continue;
          del34[0] = x[jpt][0]-x[lpt][0];
          del34[1] = x[jpt][1]-x[lpt][1];
          del34[2] = x[jpt][2]-x[lpt][2];
          rsq = del34[0]*del34[0] + del34[1]*del34[1] + del34[2]*del34[2];
          r34 = sqrt(rsq);
          cos234 = (del32[0]*del34[0] + del32[1]*del34[1] +
                    del32[2]*del34[2]) / (r32*r34);
          cos234 = MIN(cos234,1.0);
          cos234 = MAX(cos234,-1.0);
          sin234 = sqrt(1.0 - cos234*cos234);
          if (sin234 < TOL) continue;
          w34 = Sp(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],dw34);
          delil[0] = del23[0] + del34[0];
          delil[1] = del23[1] + del34[1];
          delil[2] = del23[2] + del34[2];
          ril2 = delil[0]*delil[0] + delil[1]*delil[1] + delil[2]*delil[2];
          ril=sqrt(ril2);
          rjl2 = r34*r34;

          rjl = r34;
          rjl2 = r34*r34;
          costmp = 0.5*(rij2+rjl2-ril2)/rij/rjl;
          tspijl = Sp2(costmp,thmin,thmax,dtsijl);
          dtsijl = -dtsijl; //need minus sign
          cross321[0] = (del32[1]*del21[2])-(del32[2]*del21[1]);
          cross321[1] = (del32[2]*del21[0])-(del32[0]*del21[2]);
          cross321[2] = (del32[0]*del21[1])-(del32[1]*del21[0]);
          cross321mag = sqrt(cross321[0]*cross321[0]+
                             cross321[1]*cross321[1]+
                             cross321[2]*cross321[2]);
          cross234[0] = (del23[1]*del34[2])-(del23[2]*del34[1]);
          cross234[1] = (del23[2]*del34[0])-(del23[0]*del34[2]);
          cross234[2] = (del23[0]*del34[1])-(del23[1]*del34[0]);
          cross234mag = sqrt(cross234[0]*cross234[0]+
                             cross234[1]*cross234[1]+
                             cross234[2]*cross234[2]);
          cwnum = (cross321[0]*cross234[0]) +
            (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
          cwnom = r21*r34*r32*r32*sin321*sin234;
          cw = cwnum/cwnom;

          cw2 = (.5*(1.0-cw));
          ekijl = epsilonT[ktype][ltype];
          Ec = 256.0*ekijl/405.0;
          Vtors = (Ec*(powint(cw2,5)))-(ekijl/10.0);

          if (eflag) evdwl = Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);

          dndij[0] = (cross234[1]*del21[2])-(cross234[2]*del21[1]);
          dndij[1] = (cross234[2]*del21[0])-(cross234[0]*del21[2]);
          dndij[2] = (cross234[0]*del21[1])-(cross234[1]*del21[0]);

          tmpvec[0] = (del34[1]*cross321[2])-(del34[2]*cross321[1]);
          tmpvec[1] = (del34[2]*cross321[0])-(del34[0]*cross321[2]);
          tmpvec[2] = (del34[0]*cross321[1])-(del34[1]*cross321[0]);

          dndij[0] = dndij[0]+tmpvec[0];
          dndij[1] = dndij[1]+tmpvec[1];
          dndij[2] = dndij[2]+tmpvec[2];

          dndik[0] = (del23[1]*cross234[2])-(del23[2]*cross234[1]);
          dndik[1] = (del23[2]*cross234[0])-(del23[0]*cross234[2]);
          dndik[2] = (del23[0]*cross234[1])-(del23[1]*cross234[0]);

          dndjl[0] = (cross321[1]*del23[2])-(cross321[2]*del23[1]);
          dndjl[1] = (cross321[2]*del23[0])-(cross321[0]*del23[2]);
          dndjl[2] = (cross321[0]*del23[1])-(cross321[1]*del23[0]);

          dcidij = ((r23*r23)-(r21*r21)+(rjk*rjk))/(2.0*r23*r23*r21);
          dcidik = ((r21*r21)-(r23*r23)+(rjk*rjk))/(2.0*r23*r21*r21);
          dcidjk = (-rjk)/(r23*r21);
          dcjdji = ((r23*r23)-(r34*r34)+(ril*ril))/(2.0*r23*r23*r34);
          dcjdjl = ((r34*r34)-(r23*r23)+(ril*ril))/(2.0*r23*r34*r34);
          dcjdil = (-ril)/(r23*r34);

          dsidij = (-cos321/sin321)*dcidij;
          dsidik = (-cos321/sin321)*dcidik;
          dsidjk = (-cos321/sin321)*dcidjk;

          dsjdji = (-cos234/sin234)*dcjdji;
          dsjdjl = (-cos234/sin234)*dcjdjl;
          dsjdil = (-cos234/sin234)*dcjdil;

          dxidij = (r21*sin321)+(r23*r21*dsidij);
          dxidik = (r23*sin321)+(r23*r21*dsidik);
          dxidjk = (r23*r21*dsidjk);

          dxjdji = (r34*sin234)+(r23*r34*dsjdji);
          dxjdjl = (r23*sin234)+(r23*r34*dsjdjl);
          dxjdil = (r23*r34*dsjdil);

          ddndij = (dxidij*cross234mag)+(cross321mag*dxjdji);
          ddndik = dxidik*cross234mag;
          ddndjk = dxidjk*cross234mag;
          ddndjl = cross321mag*dxjdjl;
          ddndil = cross321mag*dxjdil;
          dcwddn = -cwnum/(cwnom*cwnom);
          dcwdn = 1.0/cwnom;
          dvpdcw = (-1.0)*Ec*(-.5)*5.0*powint(cw2,4) *
            w23*w21*w34*(1.0-tspjik)*(1.0-tspijl);

          Ftmp[0] = dvpdcw*((dcwdn*dndij[0])+(dcwddn*ddndij*del23[0]/r23));
          Ftmp[1] = dvpdcw*((dcwdn*dndij[1])+(dcwddn*ddndij*del23[1]/r23));
          Ftmp[2] = dvpdcw*((dcwdn*dndij[2])+(dcwddn*ddndij*del23[2]/r23));
          fi[0] = Ftmp[0];
          fi[1] = Ftmp[1];
          fi[2] = Ftmp[2];
          fj[0] = -Ftmp[0];
          fj[1] = -Ftmp[1];
          fj[2] = -Ftmp[2];

          Ftmp[0] = dvpdcw*((dcwdn*dndik[0])+(dcwddn*ddndik*del21[0]/r21));
          Ftmp[1] = dvpdcw*((dcwdn*dndik[1])+(dcwddn*ddndik*del21[1]/r21));
          Ftmp[2] = dvpdcw*((dcwdn*dndik[2])+(dcwddn*ddndik*del21[2]/r21));
          fi[0] += Ftmp[0];
          fi[1] += Ftmp[1];
          fi[2] += Ftmp[2];
          fk[0] = -Ftmp[0];
          fk[1] = -Ftmp[1];
          fk[2] = -Ftmp[2];

          Ftmp[0] = (dvpdcw*dcwddn*ddndjk*deljk[0])/rjk;
          Ftmp[1] = (dvpdcw*dcwddn*ddndjk*deljk[1])/rjk;
          Ftmp[2] = (dvpdcw*dcwddn*ddndjk*deljk[2])/rjk;
          fj[0] += Ftmp[0];
          fj[1] += Ftmp[1];
          fj[2] += Ftmp[2];
          fk[0] -= Ftmp[0];
          fk[1] -= Ftmp[1];
          fk[2] -= Ftmp[2];

          Ftmp[0] = dvpdcw*((dcwdn*dndjl[0])+(dcwddn*ddndjl*del34[0]/r34));
          Ftmp[1] = dvpdcw*((dcwdn*dndjl[1])+(dcwddn*ddndjl*del34[1]/r34));
          Ftmp[2] = dvpdcw*((dcwdn*dndjl[2])+(dcwddn*ddndjl*del34[2]/r34));
          fj[0] += Ftmp[0];
          fj[1] += Ftmp[1];
          fj[2] += Ftmp[2];
          fl[0] = -Ftmp[0];
          fl[1] = -Ftmp[1];
          fl[2] = -Ftmp[2];

          Ftmp[0] = (dvpdcw*dcwddn*ddndil*delil[0])/ril;
          Ftmp[1] = (dvpdcw*dcwddn*ddndil*delil[1])/ril;
          Ftmp[2] = (dvpdcw*dcwddn*ddndil*delil[2])/ril;
          fi[0] += Ftmp[0];
          fi[1] += Ftmp[1];
          fi[2] += Ftmp[2];
          fl[0] -= Ftmp[0];
          fl[1] -= Ftmp[1];
          fl[2] -= Ftmp[2];

          // coordination forces

          fpair = Vtors*dw21*w23*w34*(1.0-tspjik)*(1.0-tspijl) / r21;
          fi[0] -= del21[0]*fpair;
          fi[1] -= del21[1]*fpair;
          fi[2] -= del21[2]*fpair;
          fk[0] += del21[0]*fpair;
          fk[1] += del21[1]*fpair;
          fk[2] += del21[2]*fpair;

          fpair = Vtors*w21*dw23*w34*(1.0-tspjik)*(1.0-tspijl) / r23;
          fi[0] -= del23[0]*fpair;
          fi[1] -= del23[1]*fpair;
          fi[2] -= del23[2]*fpair;
          fj[0] += del23[0]*fpair;
          fj[1] += del23[1]*fpair;
          fj[2] += del23[2]*fpair;

          fpair = Vtors*w21*w23*dw34*(1.0-tspjik)*(1.0-tspijl) / r34;
          fj[0] -= del34[0]*fpair;
          fj[1] -= del34[1]*fpair;
          fj[2] -= del34[2]*fpair;
          fl[0] += del34[0]*fpair;
          fl[1] += del34[1]*fpair;
          fl[2] += del34[2]*fpair;

          // additional cut off function forces

          fcpc = -Vtors*w21*w23*w34*dtsjik*(1.0-tspijl);
          fpair = fcpc*dcidij/rij;
          fi[0] += fpair*del23[0];
          fi[1] += fpair*del23[1];
          fi[2] += fpair*del23[2];
          fj[0] -= fpair*del23[0];
          fj[1] -= fpair*del23[1];
          fj[2] -= fpair*del23[2];

          fpair = fcpc*dcidik/rik;
          fi[0] += fpair*del21[0];
          fi[1] += fpair*del21[1];
          fi[2] += fpair*del21[2];
          fk[0] -= fpair*del21[0];
          fk[1] -= fpair*del21[1];
          fk[2] -= fpair*del21[2];

          fpair = fcpc*dcidjk/rjk;
          fj[0] += fpair*deljk[0];
          fj[1] += fpair*deljk[1];
          fj[2] += fpair*deljk[2];
          fk[0] -= fpair*deljk[0];
          fk[1] -= fpair*deljk[1];
          fk[2] -= fpair*deljk[2];

          fcpc = -Vtors*w21*w23*w34*(1.0-tspjik)*dtsijl;
          fpair = fcpc*dcjdji/rij;
          fi[0] += fpair*del23[0];
          fi[1] += fpair*del23[1];
          fi[2] += fpair*del23[2];
          fj[0] -= fpair*del23[0];
          fj[1] -= fpair*del23[1];
          fj[2] -= fpair*del23[2];

          fpair = fcpc*dcjdjl/rjl;
          fj[0] += fpair*del34[0];
          fj[1] += fpair*del34[1];
          fj[2] += fpair*del34[2];
          fl[0] -= fpair*del34[0];
          fl[1] -= fpair*del34[1];
          fl[2] -= fpair*del34[2];

          fpair = fcpc*dcjdil/ril;
          fi[0] += fpair*delil[0];
          fi[1] += fpair*delil[1];
          fi[2] += fpair*delil[2];
          fl[0] -= fpair*delil[0];
          fl[1] -= fpair*delil[1];
          fl[2] -= fpair*delil[2];

          // sum per-atom forces into atom force array

          _F[ipt][0] += fi[0]; _F[ipt][1] += fi[1]; _F[ipt][2] += fi[2];
          _F[jpt][0] += fj[0]; _F[jpt][1] += fj[1]; _F[jpt][2] += fj[2];
          _F[kpt][0] += fk[0]; _F[kpt][1] += fk[1]; _F[kpt][2] += fk[2];
          _F[lpt][0] += fl[0]; _F[lpt][1] += fl[1]; _F[lpt][2] += fl[2];

          if (evflag) {
            delkl[0] = delil[0] - del21[0];
            delkl[1] = delil[1] - del21[1];
            delkl[2] = delil[2] - del21[2];
            ev_tally4(ipt,jpt,kpt,lpt,evdwl,fi,fj,fk,delil,del34,delkl);
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create REBO neighbor list from main neighbor list
   REBO neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void AIREBOFrame::REBO_neigh()
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
double AIREBOFrame::bondorder(int i, int j, double rij[3],
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
  Vector3 *x = REBO_R;
  //Vector3 *x = _R;

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

      /* 2nd term in eqn. (8) in Brenner et al (2002) */
      Etmp = Etmp+(wik*g*exp(lamdajik)); 
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS));
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  /* 3rd term in eqn. (8) in Brenner et al (2002) */
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  /* sigma-pi bond order term pij, eqn. (8) in Brenner et al (2002) */  
  pij = pow(1.0+Etmp+PijS,-0.5); 
  tmp = -0.5*pow(pij,3.0);

  // pij forces

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
    } //if (atomk != atomj)
  } //  for (k = 0; k < REBO_numneigh[i]; k++) {

  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
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
      /* 2nd term in eqn. (8) in Brenner et al (2002), i and j switched */
      Etmp = Etmp+(wjl*g*exp(lamdaijl)); 
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl)); 
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  /* 3rd term in eqn. (8) in Brenner et al (2002), i and j switched */
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);  
  /* sigma-pi bond order term pji, eqn. (8) in Brenner et al (2002) 
     indices i and j switched */  
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*pow(pji,3.0);

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
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
  /* radical part of bond order piRC, eqn. (14) in Brenner et al (2002) */
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);
  
  // piRC forces

  REBO_neighs_i = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs_i[k];
    if (atomk !=atomj) {
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
          }  // if (atomn != atomi)
        }  // for (n = 0; n < REBO_numneigh[atomk]; n++)
      }  // if (fabs(dNki) > TOL)
    }  // if (atomk !=atomj)
  }  // for (k = 0; k < REBO_numneigh[i]; k++)  
  
  // piRC forces

  REBO_neighs = REBO_firstneigh[atomj];
  for (l = 0; l < REBO_numneigh[atomj]; l++) {
    atoml = REBO_neighs[l];
    if (atoml !=atomi) {
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
   Bij* function
------------------------------------------------------------------------- */

double AIREBOFrame::bondorderLJ(int i, int j, double rij[3], double rijmag,
                               double VA, double rij0[3], double rij0mag,
                               Vector3 *f, int vflag_atom)
{
  int k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int atomi,atomj,itype,jtype,ktype,ltype,ntype;
  double rik[3], rjl[3], rkn[3],rknmag,dNki;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ;
  double Nki,Nlj,dS,lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosijldrj[3],dcosijldrl[3],dcosjikdrj[3],dwjl;
  double Tij,crosskij[3],crosskijmag;
  double crossijl[3],crossijlmag,omkijl;
  double tmppij,tmppji,dN2PIJ[2],dN2PJI[2],dN3piRC[3],dN3Tij[3];
  double bij,tmp3pij,tmp3pji,Stb,dStb;
  double r32[3],r32mag,cos321;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double dNlj;
  double PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;
  double F12[3],F23[3],F34[3],F31[3],F24[3];
  double fi[3],fj[3],fk[3],fl[3],f1[3],f2[3],f3[3],f4[4];
  double rji[3],rki[3],rlj[3],r13[3],r43[3];

  Vector3 *x = REBO_R;
  //int *type = atom->type;

  atomi = i;
  atomj = j;
  itype = species[atomi];
  jtype = species[atomj];
  wij = Sp(rij0mag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijC = nC[atomi]-(wij*kronecker(jtype,0));
  NijH = nH[atomi]-(wij*kronecker(jtype,1));
  NjiC = nC[atomj]-(wij*kronecker(itype,0));
  NjiH = nH[atomj]-(wij*kronecker(itype,1));

  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  Stb = 0.0;
  dStb = 0.0;

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
      Nki = nC[atomk]-(wik*kronecker(itype,0)) +
        nH[atomk]-(wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp += (wik*g*exp(lamdajik));
      tmp3 += (wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS));
    }
  }

  PijS = 0.0;
  dN2PIJ[0] = 0.0;
  dN2PIJ[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2PIJ);
  pij = pow(1.0+Etmp+PijS,-0.5);
  tmppij = -.5*pow(pij,3.0);
  tmp3pij = tmp3;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
      ltype = species[atoml];
      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] -
        (wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp += (wjl*g*exp(lamdaijl));
      tmp3 += (wjl*dgdN*exp(lamdaijl));
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2PJI[0] = 0.0;
  dN2PJI[1] = 0.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2PJI);
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmppji = -.5*pow(pji,3.0);
  tmp3pji = tmp3;

  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3piRC);
  Tij = 0.0;
  dN3Tij[0] = 0.0;
  dN3Tij[1] = 0.0;
  dN3Tij[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3Tij);

  Etmp = 0.0;
  if (fabs(Tij) > TOL) {
    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      ktype = species[atomk];
      if (atomk != atomj) {
        rik[0] = x[atomi][0]-x[atomk][0];
        rik[1] = x[atomi][1]-x[atomk][1];
        rik[2] = x[atomi][2]-x[atomk][2];
        rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
        cos321 = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
          (rijmag*rikmag);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);

        rjk[0] = rik[0]-rij[0];
        rjk[1] = rik[1]-rij[1];
        rjk[2] = rik[2]-rij[2];
        rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
        rij2 = rijmag*rijmag;
        rik2 = rikmag*rikmag;
        costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
        tspjik = Sp2(costmp,thmin,thmax,dtsjik);

        if (sqrt(1.0 - cos321*cos321) > sqrt(TOL)) {
          wik = Sp(rikmag,rcmin[itype][ktype],rcmaxp[itype][ktype],dwik);
          REBO_neighs_j = REBO_firstneigh[j];
          for (l = 0; l < REBO_numneigh[j]; l++) {
            atoml = REBO_neighs_j[l];
            ltype = species[atoml];
            if (!(atoml == atomi || atoml == atomk)) {
              rjl[0] = x[atomj][0]-x[atoml][0];
              rjl[1] = x[atomj][1]-x[atoml][1];
              rjl[2] = x[atomj][2]-x[atoml][2];
              rjlmag = sqrt(rjl[0]*rjl[0] + rjl[1]*rjl[1] + rjl[2]*rjl[2]);
              cos234 = -((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
                (rijmag*rjlmag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);

              ril[0] = rij[0]+rjl[0];
              ril[1] = rij[1]+rjl[1];
              ril[2] = rij[2]+rjl[2];
              ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
              rijrjl = 2.0*rijmag*rjlmag;
              rjl2 = rjlmag*rjlmag;
              costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
              tspijl = Sp2(costmp,thmin,thmax,dtsijl);

              if (sqrt(1.0 - cos234*cos234) > sqrt(TOL)) {
                wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dS);
                crosskij[0] = (rij[1]*rik[2]-rij[2]*rik[1]);
                crosskij[1] = (rij[2]*rik[0]-rij[0]*rik[2]);
                crosskij[2] = (rij[0]*rik[1]-rij[1]*rik[0]);
                crosskijmag = sqrt(crosskij[0]*crosskij[0] +
                                   crosskij[1]*crosskij[1] +
                                   crosskij[2]*crosskij[2]);
                crossijl[0] = (rij[1]*rjl[2]-rij[2]*rjl[1]);
                crossijl[1] = (rij[2]*rjl[0]-rij[0]*rjl[2]);
                crossijl[2] = (rij[0]*rjl[1]-rij[1]*rjl[0]);
                crossijlmag = sqrt(crossijl[0]*crossijl[0] +
                                   crossijl[1]*crossijl[1] +
                                   crossijl[2]*crossijl[2]);
                omkijl = -1.0*(((crosskij[0]*crossijl[0]) +
                                (crosskij[1]*crossijl[1]) +
                                (crosskij[2]*crossijl[2])) /
                               (crosskijmag*crossijlmag));
                Etmp += ((1.0-omkijl*omkijl)*wik*wjl) *
                  (1.0-tspjik)*(1.0-tspijl);
              }
            }
          }
        }
      }
    }
  }

  bij = (.5*(pij+pji))+piRC+(Tij*Etmp);
  Stb = Sp2(bij,bLJmin[itype][jtype],bLJmax[itype][jtype],dStb);
  VA = VA*dStb;

  if (dStb != 0.0) {
    tmp = tmppij;
    dN2[0] = dN2PIJ[0];
    dN2[1] = dN2PIJ[1];
    tmp3 = tmp3pij;

    // pij forces

    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      if (atomk != atomj) {
        lamdajik = 0.0;
        rik[0] = x[atomi][0]-x[atomk][0];
        rik[1] = x[atomi][1]-x[atomk][1];
        rik[2] = x[atomi][2]-x[atomk][2];
        rikmag = sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
        lamdajik = 4.0*kronecker(itype,1) *
          ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
        wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
        cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
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

    tmp = tmppji;
    tmp3 = tmp3pji;
    dN2[0] = dN2PJI[0];
    dN2[1] = dN2PJI[1];
    REBO_neighs  =  REBO_firstneigh[j];
    for (l = 0; l < REBO_numneigh[j]; l++) {
      atoml = REBO_neighs[l];
      if (atoml !=atomi) {
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
          (cosijl*((rij[0]/(rijmag*rijmag))-(rjl[0]/(rjlmag*rjlmag))));
        dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) +
          (cosijl*((rij[1]/(rijmag*rijmag))-(rjl[1]/(rjlmag*rjlmag))));
        dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) +
          (cosijl*((rij[2]/(rijmag*rijmag))-(rjl[2]/(rjlmag*rjlmag))));
        dcosijldrl[0] = (rij[0]/(rijmag*rjlmag)) +
          (cosijl*rjl[0]/(rjlmag*rjlmag));
        dcosijldrl[1] = (rij[1]/(rijmag*rjlmag)) +
          (cosijl*rjl[1]/(rjlmag*rjlmag));
        dcosijldrl[2] = (rij[2]/(rijmag*rjlmag)) +
          (cosijl*rjl[2]/(rjlmag*rjlmag));

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

        tmp2=VA*.5*(tmp*tmp3*dwjl)/rjlmag;
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

    // piRC forces

    dN3[0] = dN3piRC[0];
    dN3[1] = dN3piRC[1];
    dN3[2] = dN3piRC[2];

    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      if (atomk != atomj) {
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

    // piRC forces to J side

    REBO_neighs = REBO_firstneigh[j];
    for (l = 0; l < REBO_numneigh[j]; l++) {
      atoml = REBO_neighs[l];
      if (atoml != atomi) {
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

    if (fabs(Tij) > TOL) {
      dN3[0] = dN3Tij[0];
      dN3[1] = dN3Tij[1];
      dN3[2] = dN3Tij[2];
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
        ktype = species[atomk];
        if (atomk != atomj) {
          r21[0] = x[atom2][0]-x[atom1][0];
          r21[1] = x[atom2][1]-x[atom1][1];
          r21[2] = x[atom2][2]-x[atom1][2];
          r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
          cos321 = ((r21[0]*rij[0])+(r21[1]*rij[1])+(r21[2]*rij[2])) /
            (r21mag*rijmag);
          cos321 = MIN(cos321,1.0);
          cos321 = MAX(cos321,-1.0);
          sin321 = sqrt(1.0 - cos321*cos321);
          sink2i = 1.0/(sin321*sin321);
          rik2i = 1.0/(r21mag*r21mag);

          if (sin321 != 0.0) {
            rr = (rijmag*rijmag)-(r21mag*r21mag);
            rjk[0] = r21[0]-rij[0];
            rjk[1] = r21[1]-rij[1];
            rjk[2] = r21[2]-rij[2];
            rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
            rijrik = 2.0*rijmag*r21mag;
            rik2 = r21mag*r21mag;
            dctik = (-rr+rjk2)/(rijrik*rik2);
            dctij = (rr+rjk2)/(rijrik*rijmag*rijmag);
            dctjk = -2.0/rijrik;
            w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
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
              ltype = species[atoml];
              if (!(atoml == atomi || atoml == atomk)) {
                r34[0] = x[atom3][0]-x[atom4][0];
                r34[1] = x[atom3][1]-x[atom4][1];
                r34[2] = x[atom3][2]-x[atom4][2];
                r34mag = sqrt(r34[0]*r34[0] + r34[1]*r34[1] + r34[2]*r34[2]);
                cos234 = -1.0*((rij[0]*r34[0])+(rij[1]*r34[1]) +
                               (rij[2]*r34[2]))/(rijmag*r34mag);
                cos234 = MIN(cos234,1.0);
                cos234 = MAX(cos234,-1.0);
                sin234 = sqrt(1.0 - cos234*cos234);
                sinl2i = 1.0/(sin234*sin234);
                rjl2i = 1.0/(r34mag*r34mag);

                if (sin234 != 0.0) {
                  w34 = Sp(r34mag,rcmin[jtype][ltype],
                           rcmaxp[jtype][ltype],dw34);
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
                  dtsijl = -dtsijl; //need minus sign
                  prefactor = VA*Tij;

                  cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
                  cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
                  cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
                  cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
                  cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
                  cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

                  cwnum = (cross321[0]*cross234[0]) +
                    (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
                  cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
                  om1234 = cwnum/cwnom;
                  cw = om1234;
                  Etmp += ((1.0-(om1234*om1234))*w21*w34) *
                    (1.0-tspjik)*(1.0-tspijl);

                  dt1dik = (rik2i)-(dctik*sink2i*cos321);
                  dt1djk = (-dctjk*sink2i*cos321);
                  dt1djl = (rjl2i)-(dctjl*sinl2i*cos234);
                  dt1dil = (-dctil*sinl2i*cos234);
                  dt1dij = (2.0/(r23mag*r23mag)) -
                    (dctij*sink2i*cos321)-(dctji*sinl2i*cos234);

                  dt2dik[0] = (-r23[2]*cross234[1])+(r23[1]*cross234[2]);
                  dt2dik[1] = (-r23[0]*cross234[2])+(r23[2]*cross234[0]);
                  dt2dik[2] = (-r23[1]*cross234[0])+(r23[0]*cross234[1]);

                  dt2djl[0] = (-r23[1]*cross321[2])+(r23[2]*cross321[1]);
                  dt2djl[1] = (-r23[2]*cross321[0])+(r23[0]*cross321[2]);
                  dt2djl[2] = (-r23[0]*cross321[1])+(r23[1]*cross321[0]);

                  dt2dij[0] = (r21[2]*cross234[1]) -
                    (r34[2]*cross321[1])-(r21[1]*cross234[2]) +
                    (r34[1]*cross321[2]);
                  dt2dij[1] = (r21[0]*cross234[2]) -
                    (r34[0]*cross321[2])-(r21[2]*cross234[0]) +
                    (r34[2]*cross321[0]);
                  dt2dij[2] = (r21[1]*cross234[0]) -
                    (r34[1]*cross321[0])-(r21[0]*cross234[1]) +
                    (r34[0]*cross321[1]);

                  aa = (prefactor*2.0*cw/cwnom)*w21*w34 *
                    (1.0-tspjik)*(1.0-tspijl);
                  aaa1 = -prefactor*(1.0-(om1234*om1234)) *
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

      REBO_neighs = REBO_firstneigh[i];
      for (k = 0; k < REBO_numneigh[i]; k++) {
        atomk = REBO_neighs[k];
        if (atomk != atomj) {
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
              ntype = species[atomn];
              if (atomn !=atomi) {
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
              ntype = species[atomn];
              if (atomn != atomj) {
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
  }

  return Stb;
}

/* ----------------------------------------------------------------------
   G spline
------------------------------------------------------------------------- */

double AIREBOFrame::gSpline(double costh, double Nij, int typei,
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

double AIREBOFrame::PijSpline(double NijC, double NijH, int typei, int typej,
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

double AIREBOFrame::piRCSpline(double Nij, double Nji, double Nijconj,
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
	if ((Nij>=(double) i) && (Nij<=(double) i+1) || (Nij==(double) i)) x=i;
      for (i=0; i<piCCdom[1][1]; i++)
	if ((Nji>=(double) i) && (Nji<=(double) i+1) || (Nji==(double) i)) y=i;
      for (i=0; i<piCCdom[2][1]; i++)
	if ((Nijconj>=(double) i) && (Nijconj<=(double) i+1) || 
	    (Nijconj==(double) i)) z=i;

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

double AIREBOFrame::TijSpline(double Nij, double Nji,
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

double AIREBOFrame::kronecker(int a, int b)
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

double AIREBOFrame::Sp5th(double x, double coeffs[6], double *df)
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

double AIREBOFrame::Spbicubic(double x, double y,
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

double AIREBOFrame::Sptricubic(double x, double y, double z,
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

void AIREBOFrame::spline_init()
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
   cutoff function Sprime, eqn (20) in Brenner et al (2002)
   return cutoff and dX = derivative
------------------------------------------------------------------------- */

double AIREBOFrame::Sp(double Xij, double Xmin, double Xmax, double &dX)
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

double AIREBOFrame::Sp2(double Xij, double Xmin, double Xmax, double &dX)
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

void AIREBOFrame::read_file(char *filename)
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

void AIREBOFrame::ev_tally(int i, int j, 
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
   tally eng_vdwl and virial into global and per-atom accumulators
   called by AIREBO potential, newton_pair is always on
 ------------------------------------------------------------------------- */

void AIREBOFrame::ev_tally4(int i, int j, int k, int m, double evdwl,
                     double *fi, double *fj, double *fk,
                     double *drim, double *drjm, double *drkm)
{
    double epairfourth;
    Vector3 df, dr;
    Matrix33 vv;

    _EPOT += evdwl;
    epairfourth = 0.25 * evdwl;
    _EPOT_IND[i] += epairfourth;
    _EPOT_IND[j] += epairfourth;
    _EPOT_IND[k] += epairfourth;
    _EPOT_IND[m] += epairfourth;

    vv.clear();
    df.set(fi); dr.set(drim); vv.addnvv(0.25,dr,df);
    df.set(fj); dr.set(drjm); vv.addnvv(0.25,dr,df);
    df.set(fk); dr.set(drkm); vv.addnvv(0.25,dr,df);
    
    _VIRIAL_IND[i] += vv;
    _VIRIAL_IND[j] += vv;
    _VIRIAL_IND[k] += vv;
    _VIRIAL_IND[m] += vv;
    _VIRIAL += vv * 4.0;
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO potential, newton_pair is always on
   fpair is magnitude of force on atom I
------------------------------------------------------------------------- */

void AIREBOFrame::v_tally2(int i, int j, double fpair, double *drij)
{
  Matrix33 vv;
  Vector3 dr;

  dr.set(drij);
  vv.clear(); vv.addnvv(fpair*0.5, dr, dr);  

  _VIRIAL_IND[i] += vv;
  _VIRIAL_IND[j] += vv;
  _VIRIAL += vv * 2.0;  
}

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators
   called by AIREBO and Tersoff potential, newton_pair is always on
------------------------------------------------------------------------- */

void AIREBOFrame::v_tally3(int i, int j, int k,
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

void AIREBOFrame::v_tally4(int i, int j, int k, int m,
		    double *fi, double *fj, double *fk,
		    double *drim, double *drjm, double *drkm)
{
  Matrix33 vv;
  Vector3 dr, df;

  vv.clear();
  dr.set(drim); df.set(fi); vv.addnvv(0.25, dr, df);
  dr.set(drjm); df.set(fj); vv.addnvv(0.25, dr, df);
  dr.set(drkm); df.set(fk); vv.addnvv(0.25, dr, df);

  _VIRIAL_IND[i] += vv;
  _VIRIAL_IND[j] += vv;
  _VIRIAL_IND[k] += vv;
  _VIRIAL_IND[m] += vv;
  _VIRIAL += vv * 4.0;
}

/* ----------------------------------------------------------------------
   compute global pair virial via summing F dot r over own & ghost atoms
   at this point, only pairwise forces have been accumulated in atom->f
------------------------------------------------------------------------- */

void AIREBOFrame::virial_fdotr_compute()
{
    int n0=0,n1=_NP;
    Vector3 *x = _R;
    //Vector3 *x = REBO_R;
    Matrix33 vv;

    // sum over force on all particles including ghosts
    _VIRIAL.clear();
    for (int i = n0; i < n1; i++) 
    {
        vv.clear();
        vv.addnvv(1.0, _F[i], x[i]);
        _VIRIAL += vv;
    }
}


void AIREBOFrame::potential()
{
  int ipt;

  _EPOT=0; _VIRIAL.clear();
  for(ipt=0;ipt<_NP;ipt++) { _F[ipt].clear(); _EPOT_IND[ipt]=0; }

  airebo();
}

int AIREBOFrame::readAIREBO()
{
    int ii,jj;

    LFile::SubHomeDir(rebofile,rebofile);
    read_file(rebofile);

    spline_init();


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

    // use C-C values for these cutoffs since C atoms are biggest
    // rcmax_CC = 2.0, rcmax_CH = 1.8, rcmax_HH = 1.7

    // cut3rebo = 3 REBO distances
    cut3rebo = 3.0 * rcmax[0][0];
    INFO("cut3rebo="<<cut3rebo<<" (Angstroms)");
    
    // cutljrebosq = furthest distance from an owned atom a ghost atom can be
    //               to need its REBO neighs computed
    // interaction = M-K-I-J-L-N with I = owned and J = ghost
    //   this insures N is in the REBO neigh list of L
    //   since I-J < rcLJmax and J-L < rmax
    //double cutljrebo = rcLJmax[0][0] + rcmax[0][0];
    //cutljrebosq = cutljrebo * cutljrebo;
    
    // cutmax = furthest distance from an owned atom
    //          at which another atom will feel force, i.e. the ghost cutoff
    // for REBO term in potential:
    //   interaction = M-K-I-J-L-N with I = owned and J = ghost
    //   I to N is max distance = 3 REBO distances
    // for LJ term in potential:
    //   short interaction = M-K-I-J-L-N with I = owned, J = ghost, I-J < rcLJmax
    //   rcLJmax + 2*rcmax, since I-J < rcLJmax and J-L,L-N = REBO distances
    //   long interaction = I-J with I = owned and J = ghost
    //   cutlj*sigma, since I-J < LJ cutoff
    // cutghost = REBO cutoff used in REBO_neigh() for neighbors of ghosts

    //double cutmax = cut3rebo;
    acut = cut3rebo;
    if (ljflag) {
      //cutmax = MAX(cutmax,rcLJmax[0][0] + 2.0*rcmax[0][0]);
      //cutmax = MAX(cutmax,cutlj*sigma[0][0]);
      acut = MAX(acut,rcLJmax[0][0] + 2.0*rcmax[0][0]);
      acut = MAX(acut,cutlj*sigma[0][0]);
    }
    
    for(ii=0;ii<2;ii++)
    {
      for(jj=ii;jj<2;jj++)
      {
        cutljsq[ii][jj] = cutlj*sigma[ii][jj] * cutlj*sigma[ii][jj];
        lj1[ii][jj] = 48.0 * epsilon[ii][jj] * pow(sigma[ii][jj],12.0);
        lj2[ii][jj] = 24.0 * epsilon[ii][jj] * pow(sigma[ii][jj],6.0);
        lj3[ii][jj] = 4.0 * epsilon[ii][jj] * pow(sigma[ii][jj],12.0);
        lj4[ii][jj] = 4.0 * epsilon[ii][jj] * pow(sigma[ii][jj],6.0);
       
        //if(ii==jj) continue;
        //cutljsq[jj][ii] = cutljsq[ii][jj];
        //lj1[jj][ii] = lj1[ii][jj];
        //lj2[jj][ii] = lj2[ii][jj];
        //lj3[jj][ii] = lj3[ii][jj];
        //lj4[jj][ii] = lj4[ii][jj];
      }
    }
    INFO("cutlj[0][0]="<<cutlj*sigma[0][0]<<" (Angstroms)");

    INFO("acut="<<acut<<" (Angstroms)");
    _RLIST=acut*1.1;
    _SKIN=_RLIST-acut;
    acutsq=acut*acut;

    return 0;
}

#ifdef _TEST

/* Main Program Begins */
class AIREBOFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

