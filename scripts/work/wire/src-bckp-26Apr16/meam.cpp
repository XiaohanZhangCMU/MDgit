/*
  meam.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Jan 18 17:16:25 2007

  FUNCTION  : MD++ implementation of EAM/MEAM potential

  M.I.Baskes, PRB 46, 2727 (1992)
  "Modified embedded-atom potentials for cubic materials and impurities"

  To Do:
  
   This code is currently slower than Baskes's Fortran code dyn88.f
     by a factor of 10 (at rcut = 6 A) and
     by a factor of 5  (at rcut = 4.5 A)
     See timing data below.
     
  1. Using the intel compiler triple the speed and bring the comparison
      to Baskes's code to a factor of 3.6 (at rcut = 6 A) and
                          a factor of 1.7 (at rcut = 4.5 A)

  2. The timing of Baskes's code scales as O(rcut^2) while
      while MD++ shows a nearly O(rcut^4) behavior
      (most time is spent in screen and dscrfor)
     Need to understand how to achieve lower-order scaling
  
  3. Make a temporary array to store the relative position vectors
      and their magnitudes of neigbhoring atom pairs
      Build a temporary list as in sw.cpp
     
  4. Remove redundant calculations by merging energy and force
      calculation routines
      
  5. Test the speed of individual functions such as rscrn and dscrn
      and compare with Fortran counterparts
      
  6. Try to reduce the memory usage, improve cache hit rate

  7. Profile the code with gprof

  8. MD++ now faster than Baskes's code, but still 1.4 times slower than LAMMPS
     eliminate repeated calculation

  9. Compute energy, force and virial for random configurations using
     MD++ and LAMMPS

  Timing data for 100 evaluation of potential and forces 512 atoms
   (rcut in Angstroms,  time in seconds)

   rcut    gcc (MD++)    intel (MD++)    Baskes (g77)
  -------------------------------------------------------------------------
   6.0      66.9          21.9              5.98
   5.0      25.6          8.90              3.95
   4.5      16.3          5.62              3.27    (result essentially the same as 6.0)
   4.0      9.71          3.67              2.63    (start to be different from 6.0)

(after correcting dscrfor)
   
   rcut    gcc (MD++)    intel (MD++)    Baskes (g77)     Baskes (ifort)
  -------------------------------------------------------------------------
   6.0      --            6.08              6.02               4.44
   5.0      --            3.44              3.98               3.04
   4.5      --            2.69              3.24               2.48
   4.0      --            2.05              2.61               2.03

(after reordering cell list, tabulate neighbor list in screen)   
   rcut     intel (MD++)  MD++ (no double-count)  Baskes (w. ifort)  both with -xK -ipo without -pg
  -------------------------------------------------------------------------
   6.0      3.36             2.88                   4.00
   5.0      2.36             1.97                   2.72
   4.5      2.00             1.60                   2.32
   4.0      1.63             1.32                   1.94
                    (modified screen and kraMEAM)

cell = 10x10x60 (48000 atoms)  1 step evaluation

perfect crystal
rcut    MD++ (intel)   MD++ (no double-count) Baskes (w. ifort)  LAMMPS (g++)  MD++/LAMMPS
------------------------------------------------------------------------------------------
 6.0      2.2  (2.62)     1.83                   2.7                 1.12       0.898
 4.0      1.32 (1.40)     1.06                   --                  0.614      0.483
                    (modified screen and kraMEAM)
random position
rcut   MD++ (intel)  MD++ (no double-count)   Baskes (w. ifort)   LAMMPS (g++) MD++/LAMMPS
-------------------------------------------------------------------------------------------
 6.0    5.27             3.94                   (4.2)                  --       1.95
 4.0     --               --                     --                    --       0.89
                                         (neighborlist takes too long)
perfect crystal Stillinger Weber (SW)
MD++ (intel)    LAMMPS (g++)
-----------------------------------------------------------
   0.181           0.277

MD++ 46 neighbors
LAMMPS 23 neighbors (avoid double counting)

Link MD++ with LAMMPS lib/meam
cell = 10x10x60 -> cylinder (28620 atoms)  1 step evaluation
rcut = 6.0   MD++ (meam.cpp)       MD++(meam-lammps.cpp)
------------------------------------------------------------
  1 cpu      0.800                 0.498       
  2 cpu      0.578                 0.301

MD 1 step (measured after 11 MD steps)
rcut = 6.0   MD++ (meam.cpp)       MD++(meam-lammps.cpp)
------------------------------------------------------------
  1 cpu      1.025                 0.497
  2 cpu      0.582                 0.298

*/

#include "meam.h"

//#define _EXTRACUTOFF

int MEAMFrame::readMEAM()
{
    FILE *fp;
    char string[500], *p, *q;
    
    LFile::SubHomeDir(potfile,potfile);
    
    fp=fopen(potfile,"r");
    if(fp==NULL)
    {
        FATAL("MEAMFrame::readmeam file ("<<potfile<<") not found!");
        return -1;
    }
    fgets(string,500,fp);
    sscanf(string,"%lf",&zsmeam);
    INFO_Printf("zsmeam\t = %g\n",zsmeam);
    
    fgets(string,500,fp);
    sscanf(string,"%lf",&alphas);
    INFO_Printf("alphas\t = %g\n",alphas);
    
    fgets(string,500,fp);
    sscanf(string,"%lf %lf %lf %lf",betas,betas+1,betas+2,betas+3);
    INFO_Printf("betas\t = %g %g %g %g\n",betas[0],betas[1],betas[2],betas[3]);
    
    fgets(string,500,fp);
    sscanf(string,"%lf %lf",&esubs,&asubs);
    INFO_Printf("esubs\t = %g asubs = %g\n",esubs,asubs);

    fgets(string,500,fp);
    sscanf(string,"%lf %lf %lf %lf",ts,ts+1,ts+2,ts+3);
    INFO_Printf("ts\t = %g %g %g %g\n",ts[0],ts[1],ts[2],ts[3]);
    
    fgets(string,500,fp);
    sscanf(string,"%lf",&rozros);
    INFO_Printf("rozros\t = %g\n",rozros);
    
    fgets(string,500,fp);
    sscanf(string,"%d",&ibarr);
    INFO_Printf("ibarr\t = %d\n",ibarr);
    
    fgets(string,500,fp);
    sscanf(string,"%lf",&rcutmeam);
    INFO_Printf("rcut\t = %g\n",rcutmeam);

    fgets(string,500,fp);
    sscanf(string,"%lf %lf",&cmin,&cmax);
    INFO_Printf("cmin\t = %g  cmax = %g\n",cmin,cmax);
    
    fgets(string,500,fp);
    sscanf(string,"%lf %lf",&repuls,&attrac);
    INFO_Printf("repuls\t = %g  attrac = %g\n",repuls,attrac);
    
    fgets(string,500,fp);
    sscanf(string,"%lf",&legend);
    INFO_Printf("legend\t = %g (Legendre correction)\n",legend);
    
    fgets(string,500,fp);
    sscanf(string,"%d",&noscr);
    INFO_Printf("noscr\t = %d\n",noscr);

    fgets(string,500,fp);
    sscanf(string,"%lf",&alat);
    INFO_Printf("alat\t = %g\n",alat);

    fgets(string,500,fp);
    sscanf(string,"%lf",&ielement);
    INFO_Printf("ielement\t = %g\n",ielement);

    fgets(string,500,fp);
    p = strchr(string,'\'');  /* first quotation mark */
    if(p!=NULL)
    {
        q = strchr(p+1,'\''); /* second quotation mark */
        if(q!=NULL)
        {
            q[0]=0;
            strcpy(lat,p+1);
        }
    }
    INFO_Printf("lat\t = %s\n",lat);

    fgets(string,500,fp);
    p = strchr(string,'\'');  /* first quotation mark */
    if(p!=NULL)
    {
        q = strchr(p+1,'\''); /* second quotation mark */
        if(q!=NULL)
        {
            q[0]=0;
            strcpy(elt,p+1);
        }
    }
    INFO_Printf("elt\t = %s\n",elt);

    _RLIST = 1.1*rcutmeam;
    _SKIN = _RLIST - rcutmeam;
    rcutmeam2 = SQR(rcutmeam);

    if(strcmp(lat,"dia")==0)
    { /* diamond cubic structure */
        res = alat*sqrt(3.0)/4.0;
        INFO_Printf("diamond-cubic: res = %f\n",res);
    }
    else if (strcmp(lat,"fcc")==0)
    { /* FCC structure */
        res = alat/sqrt(2.0);
        INFO_Printf("FCC: res = %f\n",res);
    }
    else if (strcmp(lat,"bcc")==0)
    { /* BCC structure */
        res = alat*sqrt(3.0)/2.0;
        INFO_Printf("BCC: res = %f\n",res);
    }
    else
    {
        res = 0;
        FATAL("unknown lat ("<<lat<<"), res not set!");
    }
    
    return 0;
}

#ifdef _TORSION_OR_BENDING
#include "../cookies/src/meam_torsion_bending.cpp"
#else /* rhoMEAM() dscrfor() kraMEAM() affected */
void MEAMFrame::rhoMEAM()
{ /* requires screen() to compute scrnab[i][j] and phiid() */
    
    //INFO("rhoMEAM");

    int i, j, jpt, idx, jdx, itmp, jtmp, ktmp, ibar;
    Vector3 sij, rij, rcosv;
    double r2ij, rmagg, screenij, phifs, rej, rzj, drh;
    double betaj0, betaj1, betaj2, betaj3, rhj;
    double rogj1, rogj2, rogj3, rinve, rhocon;
    double t1, t2, t3, asub, esub, rozeroi, zmeam, z0meam;
    double rho0, rho1, rho2, rho3, rhobr, rhobar;

    /* lowest allowable rho value */
    rhocon = 1e-10;

    for(i=0;i<_NP;i++)
    {
        atpe2b[i]=0;     /* two-body potential energy per atom */
        atpe3b[i]=0;     /* three-body potential energy per atom */
        _EPOT_IND[i]=0;  /* potential energy per atom */
        rhotot[i]=0;     /* total electron density on each atom  */
        embf  [i]=0;     /* F(rho) for each atom */
        c8a   [i]=0;     /* a term in Eq.(8a) */
        cg8c  [i]=0;      /* terms inside the bracket of the second term in Eq.(8c) */

        /* Vector3 */
        tav[i].clear();
        ag[i].clear();
        a8b[i].clear();  /* a term in Eq.(8b) */

        /* Matrix33 */
        b8c[i].clear(); /* a term in Eq.(8c) */

        /* 3-dimensional matrix */
        d8d[i].clear(); /* a term in Eq.(8a) */

        _F[i].clear();  /* force on every atom */
        _VIRIAL_IND[i].clear();
    }
    _EPOT=0;
    _VIRIAL.clear();
    
    /*	start calculation
     *  calculate pair contribution to total energy and forces
     */

    /* do on i-particles */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue; /* -1 means not existing */                  
            jdx = species[jpt]; /* type of atom (j) */

            if(i==jpt) continue; /* there is double-counting in pair potential */
            
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            rmagg=sqrt(r2ij);

#ifdef _EXTRACUTOFF
            if(rmagg>rcutmeam) continue;
#endif            
            /* pre-computed screening value S_ij */
            screenij = scrnab[i][j];  /* computed by function screen() */

            if(screenij<scnres) continue; /* this line is in Baskes's code, but not in MDCASK */
            
            /* compute pair potential phi(r), Eq.(4), depending on atom species */
            phifs = phiid(rmagg,idx);

//            if(i==0)
//                INFO_Printf("pair_pot(%d,%d)=%f = %f * %f\n",i,jpt,phifs*screenij,phifs,screenij);
            
            atpe2b[i] += 0.5*phifs*screenij; /* add to the two-body potential for atom i */

            /* copy input variables to local variables (modify here for multi-component systems) */
            rej = res;        /* res    is set in readMEAM, computed from lat and alat */
            rzj = rozros;     /* rozros and betas are read from meamdata file */
            betaj0 = betas[0];
            betaj1 = betas[1];
            betaj2 = betas[2];
            betaj3 = betas[3];
            
            /* compute atomic electron density, rho_i^{a(l)} */
            rhj   = rhof(rmagg,betaj0,rej,rzj);  /* Eqs.(16a) and (16b), atomic electron density */
            rogj1 = rhof(rmagg,betaj1,rej,rzj);
            rogj2 = rhof(rmagg,betaj2,rej,rzj);
            rogj3 = rhof(rmagg,betaj3,rej,rzj);

            rinve = 1./rmagg;

            /* rcosv is x_{ij}^{alpha} = R_{ij}^{alpha} / R_{ij}, beneath Eq.(8) */
            rcosv = rij * rinve;

            c8a[i] += rhj*screenij;          /* Eq.(8a) */
            for(itmp=0;itmp<3;itmp++)     /* itmp is alpha */
            {
                t1 = screenij*rcosv[itmp];
                a8b[i][itmp] += rogj1*t1; /* term inside bracket of Eq.(8b) for each alpha */
                ag [i][itmp] += rogj3*t1; /* for the legendre correction, see below */
                
                for(jtmp=0;jtmp<3;jtmp++) /* jtmp is beta */
                {
                    t2 = t1*rcosv[jtmp];
                    b8c[i][itmp][jtmp] += rogj2*t2; /* terms inside the first term in Eq.(8c)
                                                       for each alpha,beta */
                    
                    for(ktmp=0;ktmp<3;ktmp++) /* ktmp is gamma */
                    {                               
                        /* terms inside the bracket of Eq.(8d) for each alpha,beta,gamma */
                        d8d[i].element[itmp][jtmp][ktmp] += rogj3*t2*rcosv[ktmp];
                    }
                }
            }

            tav[i][0] += ts[1]*rhj*screenij; /* where is this used ? */
            tav[i][1] += ts[2]*rhj*screenij;
            tav[i][2] += ts[3]*rhj*screenij;
            
            cg8c[i] += rogj2*screenij; /* terms inside the bracket of the second term in Eq.(8c) */

        }
    }

    /* compute F(rho) for each particle */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        rho0 = c8a[i]; /* Eq.(8a), rho_i^(0) */

        if(rho0<=0)
        {
            embf[i] = 0;
            rhotot[i] = rhocon;
            rhsq[i].clear();
        }
        else
        {
            t1 = tav[i][0]/rho0;
            t2 = tav[i][1]/rho0;
            t3 = tav[i][2]/rho0;            
            ibar = ibarr;       /* set local variables from input file */
            asub = asubs;       /* scaling parameter for embedding function */
            esub = esubs;       /* cohesive energy */
            rozeroi = rozros;  /* atomic density scaling parameter for multicomponent system */
            zmeam = zsmeam;     /* coordination number */

            /* z0meam = zbar(ibar,zmeam,ncrys,t1,t2,t3) */
            z0meam = zmeam;

            /* Start computing Eqs. (8b), (8c), (8d) */

            /* Second term in Eq.(8c) */
            rho2 = -(cg8c[i]*cg8c[i])/3.0;

            rho1 = rho3 = 0.0;

            /* perform legendre correction, M.I.Baskes, Mater.Sci.Eng. A 261, 165 (1999). */
            for(itmp=0;itmp<3;itmp++)
                rho3 -= legend*SQR(ag[i][itmp]); /* legend = 0 or 3/5 */

            for(itmp=0;itmp<3;itmp++)
            {
                rho1 += SQR(a8b[i][itmp]); /* square and sum up all contributions to Eq.(8b) */
                for(jtmp=0;jtmp<3;jtmp++)  
                {
                    rho2 += SQR(b8c[i][itmp][jtmp]); /* add contributions in the first term of Eq.(8c) */
                    for(ktmp=0;ktmp<3;ktmp++)
                    {
                        rho3 += SQR(d8d[i].element[itmp][jtmp][ktmp]); /* add contributions to Eq.(8d) */
                    }
                }
            }

            /* Save values for later calculations of forces, Eq.(8b,8c,8d), Eq(8a) is stored in c8a */
            rhsq[i].set(rho1,rho2,rho3);

            /* Eq. 9(a), sum over three terms to compute rhobar_i */
            drh   = t1*rho1 + t2*rho2 + t3*rho3; /* drh is Gamma in MSMSE 13, 1309 (2005), Eq.(6) */
            rhobr = bar(rho0,drh,ibar,zmeam*rozeroi,dang1+i,dang2+i);
            rhobar= rhobr / (z0meam*rozeroi);

            if(rhobr<0)
            {
                ERROR("rhobr("<<i<<") = "<<rhobr);
            }
            rhotot[i] = rhobr; /* save the un-scaled density to memory */
            
            embf[i] = frhoi(rhobar,asub,esub); /* F(rho) */
            
            atpe3b[i] = embf[i]; /* three-body potential energy for atom i */
        }
    }

    for(i=0;i<_NP;i++)
    {
        _EPOT_IND[i] = atpe2b[i] + atpe3b[i];
        _EPOT += _EPOT_IND[i];
    }
}

void MEAMFrame::dscrfor()
{
    int i, j, k, jpt, kpt, idx, jdx, kdx, ibar, ib, itmp, jtmp, ktmp;
    Vector3 sij, rij, sik, rik, rjk, rkm, r2ij_p, r2ij_m, r2ik_p, r2ik_m, r2jk_p, r2jk_m;
    double r2ij, r2jk, r2ik, screenij;
    double hhmeam, h2meam;
    Vector3 rcosv, scrpk, scrmk, dscrnabk, df;
    double t1, t2, t3, s1, s2, s3, asub, esub, rozeroi, zmeam, z0meam, re, roz;
    double rho0, dfrho, drbd0, drbdAA, tmp0, rmagg, rinve;
    double rozero, beta0, beta1, beta2, beta3, rhs1, rhs2, rhs3, rhs4;
    double factor, ddd, phiij, scr1, drho0s, drho1s, drho1sg, drho2s, drho3s;
    double tmp1, tmp2, tmp3, tmp1g, term2, rml, rmln;
//    int icount;
    

    /* compute the forces on atoms due to the screen functions */
    //INFO("dscrfor");

//    icount = 0;
    
    if(noscr)
    {
        INFO("no need to call dscrfor");
        return;
    }

    hhmeam = SQR(hmeam);  /* hmeam is the step size for numerical differentiation */
    h2meam = 2.0 * hmeam;    

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */

        if(c8a[i]<=0) /* c8a is rho_i^(0), Eq.(8a) */
            continue;

        /* define weighting factors for atomic densities */
        t1 = tav[i][0]/c8a[i];
        t2 = tav[i][1]/c8a[i];
        t3 = tav[i][2]/c8a[i];
        ibar = ibarr;
        ib = abs(ibar);
        asub = asubs;
        esub = esubs;
        rozeroi = rozros;
        zmeam = zsmeam;
        z0meam = zmeam;
        s1 = 0.0;
        s2 = 0.0;
        s3 = 0.0;

        /* divide rho by z0(rozeroi) here */
        roz = rhotot[i] / (z0meam*rozeroi);
        rho0 = c8a[i];
        if(roz<0)
        {
            ERROR("roz["<<i<<"] = "<<roz);
        }
        dfrho = dfrhoi(roz,asub,esub)/(z0meam*rozeroi);

        /* d rhobar / d rho0 */
        drbd0 = rhotot[i]/rho0 + rho0*dang1[i];

        /* drhobar / d AA */
        drbdAA = rho0*dang2[i];

        tmp0 = (2.0/3.0) * cg8c[i];
        
        
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue; /* -1 means not existing */        
            jdx = species[jpt]; /* type of atom (j) */

            if(i==jpt) continue; /* there is double-counting in pair potential */
            
            /* scrnab has been calculated by screen() before */
            screenij=scrnab[i][j];
            if((screenij<=scnres)||(screenij>=1.0-scnres)) continue;
            
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            rmagg=sqrt(r2ij);

#ifdef _EXTRACUTOFF
            if(rmagg>rcutmeam) continue;
#endif            
            rinve = 1./rmagg;

            /* rcosv is x_{ij}^{alpha} = R_{ij}^{alpha} / R_{ij}, beneath Eq.(8) */
            rcosv = rij * rinve;

            re = res;
            rozero = rozros;
            beta0 = betas[0];
            beta1 = betas[1];
            beta2 = betas[2];
            beta3 = betas[3];

            rhs1 = rhof(rmagg,beta0,re,rozero);
            rhs2 = rhof(rmagg,beta1,re,rozero);
            rhs3 = rhof(rmagg,beta2,re,rozero);
            rhs4 = rhof(rmagg,beta3,re,rozero);
            
            if(ibar<=0)
            {
                factor = 0.0;
            }
            else
            {
                factor = SQR(rho0/zmeam);
            }

            s1 = s2 = s3 = 0.0;

            ddd = (ts[1]-t1)*(rhsq[i][0]-s1*factor) +
                  (ts[2]-t2)*(rhsq[i][1]-s2*factor) +
                  (ts[3]-t3)*(rhsq[i][2]-s3*factor);
            
            phiij = phiid(rmagg,idx);

//            icount++;            
            /* three-body loop */
            for(k=0;k<nn[i];k++)
            {                
                kpt = nindex[i][k];
                if(fixed[kpt]==-1) continue; /* -1 means not existing */        
                kdx = species[kpt];

                if(kpt==i)   continue;
                if(kpt==jpt) continue;

                sik=_SR[kpt]-_SR[i];
                sik.subint();
                rik=_H*sik;
                r2ik=rik.norm2();

#ifdef _EXTRACUTOFF
                if(r2ik>rcutmeam2) continue;
#endif                
                rjk = rik - rij;
                r2jk=rjk.norm2();

#ifdef _EXTRACUTOFF
                if(r2jk>rcutmeam2) continue;
#endif                
                /* find proper coordinates for compilation of VIRIAL */
                rkm = rik*2.0 - rij;

                /* find screening */
                scr1 = dscrn(r2ij,r2ik,r2jk,cmin,cmax);

                /* See if derivative of screening term is non-zero */
                if((scr1<scnres)||(scr1>1.0-scnres)) continue;

                /* perform numerical derivatives with respect to atom k */
                r2ik_p = rik*( h2meam); r2ik_p += r2ik + hhmeam;
                r2ik_m = rik*(-h2meam); r2ik_m += r2ik + hhmeam;
                
                r2jk_p = rjk*( h2meam); r2jk_p += r2jk + hhmeam;
                r2jk_m = rjk*(-h2meam); r2jk_m += r2jk + hhmeam;

                scrpk.x = dscrn(r2ij,r2ik_p.x,r2jk_p.x,cmin,cmax);
                scrpk.y = dscrn(r2ij,r2ik_p.y,r2jk_p.y,cmin,cmax);
                scrpk.z = dscrn(r2ij,r2ik_p.z,r2jk_p.z,cmin,cmax);

                scrmk.x = dscrn(r2ij,r2ik_m.x,r2jk_m.x,cmin,cmax);
                scrmk.y = dscrn(r2ij,r2ik_m.y,r2jk_m.y,cmin,cmax);
                scrmk.z = dscrn(r2ij,r2ik_m.z,r2jk_m.z,cmin,cmax);
//                icount+=6;
                
                /* find derivative */
                dscrnabk = scrpk - scrmk;
                dscrnabk *= (screenij/(scr1*hmeam*2));
                
                /* Initialize accumulation variables */
                drho0s  = rhs1;
                drho1s  = 0.0;
                drho1sg = 0.0;
                drho2s  = -rhs3 * tmp0;
                drho3s  = 0.0;

                for(itmp=0;itmp<3;itmp++)
                {
                    tmp1 = 2.0 * a8b[i][itmp];
                    tmp1g = 2.0 * ag[i][itmp];
                    drho1s += tmp1*rhs2*rcosv[itmp];
                    drho1sg += tmp1g*rhs4*rcosv[itmp];

                    for(jtmp=0;jtmp<3;jtmp++)
                    {
                        rml = rcosv[itmp]*rcosv[jtmp];
                        tmp2 = 2.0 * b8c[i][itmp][jtmp];
                        drho2s += tmp2 * rhs3 * rml;

                        for(ktmp=0;ktmp<3;ktmp++)
                        {
                            rmln = rml * rcosv[ktmp];
                            tmp3 = 2.0 * d8d[i].element[itmp][jtmp][ktmp];
                            drho3s += tmp3 * rhs4 * rmln;;
                        }
                    }
                }

                term2 = t1*drho1s + t2*drho2s
                      + t3*(drho3s-legend*drho1sg) + ddd*drho0s/rho0;

                /* sum forces on the central atom */
                df = dscrnabk * (0.5*phiij + dfrho * (drbd0*drho0s + drbdAA*term2));

                _F[kpt] -= df;

                _VIRIAL.addnvv(-0.5,df,rkm);
                _VIRIAL_IND[kpt].addnvv(-0.5,df,rkm);
            }
        }
    }
//    INFO_Printf("dscrfor: icount=%d\n",icount);
}


void MEAMFrame::kraMEAM()
{
    int i, j, ia, jpt, idx, jdx, ibar, ib, kk, itmp, jtmp, ktmp;
    double t1, t2, t3, asub, esub, rozeroi, zmeam, z0meam, s1, s2, s3, roz, rho0;
    double dfrho, drbd0, drbdAA, screenij, r2ij, rmagg, rinve;
    double beta0, beta1, beta2, beta3, re, rozero, rhs1, rhs2, rhs3, rhs4, drhs1, drhs2, drhs3, drhs4;
    double factor, ddd, tmp0;
    double drho0, drho1, drho1g, drho1sg, drho2, drho3;
    double drho0s, drho1s, drho2s, drho3s, tmp1, tmp2, tmp3, tmp4, tmp1g, term1, term2;
    double t1j, t2j, t3j, rozj, rho0j, dfrhoj, drbd0j, drbdAAj;
    double drho0j, drho1j, drho1gj, drho1sgj, drho2j, drho3j, dddj, tmp0j;
    double drho0sj, drho1sj, drho2sj, drho3sj, tmp1j, tmp2j, tmp3j, tmp4j, tmp1gj, term1j, term2j;
    double rml, rmln, drmllm, df, phip, phim, dphif, phiij;
    Vector3 sij, rij, rcosv, dfi, dfj, dfip, dfjp;
    Matrix33 drcos;
    
    //INFO("kraMEAM");

    /* reset force vectors */
    for(i=0;i<_NP;i++)
        _F[i].clear();

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */

        /* define weighting factors for atomic densities */
        t1 = tav[i][0]/c8a[i];
        t2 = tav[i][1]/c8a[i];
        t3 = tav[i][2]/c8a[i];
        ibar = ibarr;
        ib = abs(ibar);
        asub = asubs;
        esub = esubs;
        rozeroi = rozros;
        zmeam = zsmeam;
        z0meam = zmeam;
        s1 = 0.0;
        s2 = 0.0;
        s3 = 0.0;

        /* divide rho by z0(rozeroi) here */
        roz = rhotot[i] / (z0meam*rozeroi);
        rho0 = c8a[i];
        if(roz<0)
        {
            ERROR("roz["<<i<<"] = "<<roz);
        }
        dfrho = dfrhoi(roz,asub,esub)/(z0meam*rozeroi);

        /* d rhobar / d rho0 */
        drbd0 = rhotot[i]/rho0 + rho0*dang1[i];

        /* drhobar / d AA */
        drbdAA = rho0*dang2[i];

        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue; /* -1 means not existing */        
            jdx = species[jpt]; /* type of atom (j) */

            //if(i==jpt) continue; /* there is double-counting in Baskes's code */
            if(i>=jpt) continue; /* avoid double-counting here (not in Baskes's code */

            screenij = scrnab[i][j];
            if(screenij<=scnres) continue;
            
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            
            r2ij=rij.norm2();
            rmagg=sqrt(r2ij);
            
#ifdef _EXTRACUTOFF
            if(rmagg>rcutmeam) continue;
#endif
                
            rinve = 1./rmagg;

            /* rcosv is x_{ij}^{alpha} = R_{ij}^{alpha} / R_{ij}, beneath Eq.(8) */
            rcosv = rij * rinve;

            /* Find direction cosines */
            drcos = dcmij(rij,rmagg);
                        
            /* repeat same procedure above for atom j */
            /* defining weighting factors for atomic densities */
            t1j = tav[jpt][0]/c8a[jpt];
            t2j = tav[jpt][1]/c8a[jpt];
            t3j = tav[jpt][2]/c8a[jpt];

            /* divide rho by z0(rozeroi) here */
            rozj = rhotot[jpt] / (z0meam*rozeroi);
            rho0j = c8a[jpt];
            if(rozj<0)
            {
                ERROR("roz["<<jpt<<"] = "<<rozj);
            }
            dfrhoj = dfrhoi(rozj,asub,esub)/(z0meam*rozeroi);

            /* d rhobar / d rho0 */
            drbd0j = rhotot[jpt]/rho0j + rho0j*dang1[jpt];
            
            /* drhobar / d AA */
            drbdAAj = rho0j*dang2[jpt];

            re = res;
            rozero = rozros;
            beta0 = betas[0];
            beta1 = betas[1];
            beta2 = betas[2];
            beta3 = betas[3];

            rhs1 = rhof(rmagg,beta0,re,rozero);
            rhs2 = rhof(rmagg,beta1,re,rozero);
            rhs3 = rhof(rmagg,beta2,re,rozero);
            rhs4 = rhof(rmagg,beta3,re,rozero);
            
            drhs1 = drhof(beta0,re,rhs1);
            drhs2 = drhof(beta1,re,rhs2);
            drhs3 = drhof(beta2,re,rhs3);
            drhs4 = drhof(beta3,re,rhs4);

            if(ibar<=0)
                factor = 0.0;
            else
                factor = SQR(rho0/zmeam);

            ddd = (ts[1]-t1)*(rhsq[i][0]-s1*factor) +
                  (ts[2]-t2)*(rhsq[i][1]-s2*factor) +
                  (ts[3]-t3)*(rhsq[i][2]-s3*factor);
            
            tmp0 = (2.0/3.0) * cg8c[i];
            
            dddj = (ts[1]-t1j)*(rhsq[jpt][0]-s1*factor) +
                   (ts[2]-t2j)*(rhsq[jpt][1]-s2*factor) +
                   (ts[3]-t3j)*(rhsq[jpt][2]-s3*factor);
            
            tmp0j = (2.0/3.0) * cg8c[jpt];

            for(kk=0;kk<3;kk++)
            {
                drho0 = drhs1*screenij*rcosv[kk];
                drho1 = 0.0;
                drho1g = 0.0;
                drho1sg = 0.0;
                drho2 = -tmp0*drhs3*screenij*rcosv[kk];
                drho3 = 0.0;
                drho0s = rhs1;
                drho1s = 0.0;
                drho2s = -rhs3*tmp0;           
                drho3s = 0.0;

                drho0j = drhs1*screenij*rcosv[kk];
                drho1j = 0.0;
                drho1gj = 0.0;
                drho1sgj = 0.0;
                drho2j = -tmp0j*drhs3*screenij*rcosv[kk];
                drho3j = 0.0;
                drho0sj = rhs1;
                drho1sj = 0.0;
                drho2sj = -rhs3*tmp0j;           
                drho3sj = 0.0;
                
                for(itmp=0;itmp<3;itmp++)
                {
                    tmp1 = 2.0*a8b[i][itmp];
                    tmp1g = 2.0*ag[i][itmp];
                    
                    drho1 += tmp1*(drhs2*rcosv[kk]*rcosv[itmp] +
                                   rhs2*drcos[kk][itmp])*screenij;
                    drho1g += tmp1g*(drhs4*rcosv[kk]*rcosv[itmp] +
                                     rhs4*drcos[kk][itmp])*screenij;
                    drho1s += tmp1*rhs2*rcosv[itmp];
                    drho1sg += tmp1g*rhs4*rcosv[itmp];
                    
                    tmp1j = 2.0*a8b[jpt][itmp];
                    tmp1gj = 2.0*ag[jpt][itmp];
                    
                    drho1j += tmp1j*(drhs2*rcosv[kk]*rcosv[itmp] +
                                   rhs2*drcos[kk][itmp])*screenij;
                    drho1gj += tmp1gj*(drhs4*rcosv[kk]*rcosv[itmp] +
                                     rhs4*drcos[kk][itmp])*screenij;
                    drho1sj += tmp1j*rhs2*rcosv[itmp];
                    drho1sgj += tmp1gj*rhs4*rcosv[itmp];
                    
                    for(jtmp=0;jtmp<3;jtmp++)
                    {
                        rml = rcosv[itmp]*rcosv[jtmp];
                        drmllm = drcos[kk][itmp]*rcosv[jtmp]+
                                 rcosv[itmp]*drcos[kk][jtmp];
                        tmp2 = 2.0*b8c[i][itmp][jtmp];
                        drho2 += tmp2*(drhs3*rcosv[kk]*rml + rhs3*drmllm)*screenij;
                        drho2s += tmp2*rhs3*rml;

                        tmp2j = 2.0*b8c[jpt][itmp][jtmp];
                        drho2j += tmp2j*(drhs3*rcosv[kk]*rml + rhs3*drmllm)*screenij;
                        drho2sj += tmp2j*rhs3*rml;
                                                
                        for(ktmp=0;ktmp<3;ktmp++)
                        {
                            rmln = rml*rcosv[ktmp];
                            tmp3 = 2.0*d8d[i].element[itmp][jtmp][ktmp];
                            drho3 += tmp3*(drhs4*rcosv[kk]*rmln
                                    + rhs4*(drmllm*rcosv[ktmp] + rml*drcos[kk][ktmp]))
                                    * screenij;
                            drho3s += tmp3*rhs4*rmln;

                            tmp3j = 2.0*d8d[jpt].element[itmp][jtmp][ktmp];
                            drho3j += tmp3j*(drhs4*rcosv[kk]*rmln
                                    + rhs4*(drmllm*rcosv[ktmp] + rml*drcos[kk][ktmp]))
                                    * screenij;
                            drho3sj += tmp3j*rhs4*rmln;                            
                        }
                    }
                }

                term1 = t1*drho1 + t2*drho2 + t3*(drho3-legend*drho1g) + ddd*drho0/rho0;

                term2 = t1*drho1s + t2*drho2s + t3*(drho3s-legend*drho1sg) + ddd*drho0s/rho0;

                term1j = t1*drho1j - t2*drho2j + t3*(drho3j-legend*drho1gj) + dddj*drho0j/rho0j;

                term2j = t1*drho1sj - t2*drho2sj + t3*(drho3sj-legend*drho1sgj) + dddj*drho0sj/rho0j;
                
                /* Compile force on atom i */
                for(ia=0;ia<nn[jpt];ia++)
                    if(nindex[jpt][ia]==i) break;
                if(ia==nn[jpt])
                    FATAL(jpt<<"in neighbor list of "<<i<<" but the reverse it not true");

                tmp4 = dscrnab[i][3*j+kk];
                df = dfrho*((drho0-drho0s*tmp4)*drbd0 + (term1-term2*tmp4)*drbdAA);
                dfi[kk]=df;

                tmp4 = dscrnab[jpt][3*ia+kk];
                df = dfrho*((drho0+drho0s*tmp4)*drbd0 + (term1+term2*tmp4)*drbdAA);
                dfj[kk]=df;
                
                tmp4j = dscrnab[jpt][3*ia+kk];
                df = dfrhoj*((-drho0j-drho0sj*tmp4j)*drbd0j + (term1j+term2j*tmp4j)*drbdAAj);
                dfip[kk]=df;

                tmp4j = dscrnab[i][3*j+kk];
                df = dfrhoj*((-drho0j+drho0sj*tmp4j)*drbd0j + (term1j-term2j*tmp4j)*drbdAAj);
                dfjp[kk]=df;                
            }
            /* Baskes's code has double counting, hence a factor of 0.5 */
            /*
            _F[i] += dfi;
            _VIRIAL.addnvv(-0.5,dfi,rij);
            _VIRIAL_IND[i].addnvv(-0.5,dfi,rij); 

            _F[jpt] -= dfj;
            _VIRIAL.addnvv(-0.5,dfj,rij);
            _VIRIAL_IND[jpt].addnvv(-0.5,dfj,rij);
            */
            /* Avoid double counting here, not in Baskes's code */
            _F[i]   += dfi;   _F[jpt] -= dfj;
            _F[i]   -= dfjp;  _F[jpt] += dfip;
            
            _VIRIAL.addnvv(-0.5,dfi,rij);
            _VIRIAL_IND[i].addnvv(-0.5,dfi,rij); 

            _VIRIAL.addnvv(-0.5,dfj,rij);
            _VIRIAL_IND[jpt].addnvv(-0.5,dfj,rij);

            _VIRIAL.addnvv( 0.5,dfjp,rij);
            _VIRIAL_IND[i].addnvv( 0.5,dfjp,rij); 

            _VIRIAL.addnvv( 0.5,dfip,rij);
            _VIRIAL_IND[jpt].addnvv( 0.5,dfip,rij);

            /* Add forces due to two body potential */
            phip = phiid(rmagg+hmeam,idx);
            phim = phiid(rmagg-hmeam,idx);
   
            dphif=(phip-phim)/(2.*hmeam);

            phiij=phiid(rmagg,idx);

            dfi = rcosv*(dphif*screenij);
            dfi[0] -= phiij*dscrnab[i][3*j+0];
            dfi[1] -= phiij*dscrnab[i][3*j+1];
            dfi[2] -= phiij*dscrnab[i][3*j+2];

            /* Baskes's code has double counting, hence factor of 0.5 */
            _F[i] += dfi;
            _VIRIAL.addnvv(-0.5,dfi,rij);
            _VIRIAL_IND[i].addnvv(-0.5,dfi,rij);

            /* The following lines are not in Baskes's code
             * avoid double counting
             */
            dfi = rcosv*(dphif*screenij*(-1.0));
            dfi[0] -= phiij*dscrnab[jpt][3*ia+0];
            dfi[1] -= phiij*dscrnab[jpt][3*ia+1];
            dfi[2] -= phiij*dscrnab[jpt][3*ia+2];
            _F[jpt] += dfi;
            _VIRIAL.addnvv( 0.5,dfi,rij);
            _VIRIAL_IND[jpt].addnvv( 0.5,dfi,rij);
        }
    }
}
#endif

void MEAMFrame::screen()
{
    int i, j, k, jpt, kpt, idx, jdx, kdx, iid;
    Vector3 sij, rij, sik, rik, rjk;
    double r2ij, r2jk, r2ik, screenij;
    Vector3 *rij_table; double *r2ij_table;
//    int icount;
    
    /* compute all screen functions S_ij */
    //INFO("screen");

    /* three-body screening function, results stored for each pair i-j
     *  Baskes, Angelo, Bisson, Modell.Simul.Mater.Sci.Eng. 2, 505 (1994)
     */

//    icount = 0;
    rij_table  = (Vector3 *) malloc(sizeof(Vector3)*_NNM);
    r2ij_table = (double  *) malloc(sizeof(double )*_NNM);
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
       
        /* pre-compute and store neighbor positions */
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            jdx = species[jpt]; /* type of atom (j) */
            if(i==jpt) continue; 
            
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            rij_table[j] = rij;
            r2ij_table[j] = r2ij;
        }

        
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue; /* -1 means not existing */        
            jdx = species[jpt]; /* type of atom (j) */

            //if(i==jpt) continue; /* there is double-counting in Baskes's code */
            if(i>=jpt) continue; /* avoid double-counting here (not in Baskes's code */

            iid=-1;
            for(k=0;k<nn[jpt];k++)
            {
                if(nindex[jpt][k]==i)
                {
                    iid = k;
                    break;
                }
            }
            if(iid==-1)
                FATAL("screen: i="<<i<<" j="<<j<<" jpt="<<jpt<<" iid="<<iid);            

            
            scrnab[i][j] = 0;  /* reset */
            scrnab[jpt][iid] = 0;
            
            rij = rij_table[j];
            r2ij = r2ij_table[j];

#ifdef _EXTRACUTOFF
            /* this line is not in Baskes's code */
            if(r2ij>rcutmeam2) continue;
#endif
            /* two-body screening */
            screenij = rscrn(r2ij);
//            icount ++;

            if((noscr==1)||(screenij<=scnres)) continue;
            for(k=0;k<nn[i];k++)
            {
                if(k==j)   continue;
                kpt = nindex[i][k];
                if(fixed[kpt]==-1) continue; /* -1 means not existing */        
                kdx = species[kpt];

                if(kpt==i)   continue;
                if(kpt==jpt) continue;

                rik = rij_table[k];
                r2ik = r2ij_table[k];
                
#ifdef _EXTRACUTOFF
                /* this line is not in Baskes's code */
                if(r2ik>rcutmeam2) continue;
#endif
                
                rjk = rik - rij;
                r2jk=rjk.norm2();

#ifdef _EXTRACUTOFF
                /* this line is not in Baskes's code */
                if(r2jk>rcutmeam2) continue;
#endif                
                /* three-body screening contribution *
                 *  S_ik = prod_k S_ijk, Eq.(A1)     */
                screenij *= dscrn(r2ij,r2ik,r2jk,cmin,cmax);
//                icount ++;
                
                /* because atoms in neighbor lists are ordered
                 * differently, this code and Baskes's code have
                 * different number of evaluations (count)
                 */
                if(screenij<=scnres) break;
            }
            scrnab[i][j] = screenij;
            scrnab[jpt][iid] = screenij;
        }        
    }
//    INFO_Printf("screen: icount=%d\n",icount);
    free(rij_table);
    free(r2ij_table);
}

void MEAMFrame::dscreen()
{
    int i, j, k, jpt, kpt, idx, jdx, kdx;
    Vector3 sij, rij, sik, rik, rjk, r2ij_p, r2ij_m, r2ik_p, r2ik_m;
    double r2ij, r2jk, r2ik, screenij;
    double hhmeam, h2meam, scrpx, scrpy, scrpz, scrmx, scrmy, scrmz;
//    int icount;
    
    /* compute the derivatives of the screen functions, dscrnab */
    //INFO("dscreen");

    if(noscr)
    {
        INFO("no need to call dscreen");
        return;
    }

    hhmeam = SQR(hmeam);  /* hmeam is the step size for numerical differentiation */
    h2meam = 2.0 * hmeam;    

//    icount = 0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue; /* -1 means not existing */        
        /* modify here for binary systems (0/1) */
        idx = species[i]; /* type of atom (i) */
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            if(fixed[jpt]==-1) continue; /* -1 means not existing */        
            jdx = species[jpt]; /* type of atom (j) */

            if(i==jpt) continue; /* there is double-counting in pair potential */
            
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();

#ifdef _EXTRACUTOFF
            /* this line is not in Baskes's code */
            if(r2ij>rcutmeam2) continue;
#endif            
            /* scrnab has been calculated by screen() before */
            screenij=scrnab[i][j];

            /* if screen is too big or too small,
             * set derivatives to zero and exit from neighbor loop */
            if((screenij>1-scnres)||(screenij<scnres))
            {
                dscrnab[i][j*3] = dscrnab[i][j*3+1] = dscrnab[i][j*3+2] = 0;
                continue;
            }

            /* perform numerical differentiation */
            r2ij_p = rij*(-h2meam);  r2ij_p += hhmeam + r2ij;
            r2ij_m = rij*( h2meam);  r2ij_m += hhmeam + r2ij;
            
            scrpx = rscrn(r2ij_p.x);
            scrpy = rscrn(r2ij_p.y);
            scrpz = rscrn(r2ij_p.z);

            scrmx = rscrn(r2ij_m.x);
            scrmy = rscrn(r2ij_m.y);
            scrmz = rscrn(r2ij_m.z);

//            icount += 6;
                        
            /* three-body loop */
            for(k=0;k<nn[i];k++)
            {
                kpt = nindex[i][k];
                if(fixed[kpt]==-1) continue; /* -1 means not existing */        
                kdx = species[kpt];

                if(kpt==i)   continue;
                if(kpt==jpt) continue;

                sik=_SR[kpt]-_SR[i];
                sik.subint();
                rik=_H*sik;
                r2ik=rik.norm2();

#ifdef _EXTRACUTOFF
                if(r2ik>rcutmeam2) continue;
#endif
                rjk = rik - rij;
                r2jk=rjk.norm2();

#ifdef _EXTRACUTOFF
                if(r2jk>rcutmeam2) continue;
#endif
                /* three-body screening contribution *
                 *  S_ik = prod_k S_ijk, Eq.(A1)     */
                
                r2ik_p = rik*(-h2meam); r2ik_p += hhmeam + r2ik;
                r2ik_m = rik*( h2meam); r2ik_m += hhmeam + r2ik;
                
                scrpx *= dscrn(r2ij_p.x,r2ik_p.x,r2jk,cmin,cmax);
                scrpy *= dscrn(r2ij_p.y,r2ik_p.y,r2jk,cmin,cmax);
                scrpz *= dscrn(r2ij_p.z,r2ik_p.z,r2jk,cmin,cmax);

                scrmx *= dscrn(r2ij_m.x,r2ik_m.x,r2jk,cmin,cmax);
                scrmy *= dscrn(r2ij_m.y,r2ik_m.y,r2jk,cmin,cmax);
                scrmz *= dscrn(r2ij_m.z,r2ik_m.z,r2jk,cmin,cmax);

//                icount += 6;
            }
            dscrnab[i][j*3+0] = (scrpx-scrmx)/(2*hmeam);
            dscrnab[i][j*3+1] = (scrpy-scrmy)/(2*hmeam);
            dscrnab[i][j*3+2] = (scrpz-scrmz)/(2*hmeam);
        }        
    }
//    INFO_Printf("dscreen: icount=%d\n",icount);
}

/* MEAM library functions */
double MEAMFrame::phiid(double r, int it)
{
    double charge, re, alpha, xmjc, phizbl, frac, fdimer;
    double rozero, esub, repul, attra;
    double asub, ztmp, t1, t2, t3, z0meam, rh0, roz, zdimer;
    double beta0, beta1, beta2, beta3, rh1, rh2, rh3, arg, dum1, dum2;
    double phi;
    int ibar;

    //INFO_Printf("MEAMFrame::phiid: r = %f  it=%d\n",r, it);
    
    charge = ielement;
    re = res;
    alpha = alphas;

    beta0=betas[0];
    beta1=betas[1];
    beta2=betas[2];
    beta3=betas[3];
    t1   =ts[1];
    t2   =ts[2];
    t3   =ts[3];
    
    rozero=rozros;
    esub  =esubs;
    repul =repuls;
    attra =attrac;
    asub  =asubs;
    ibar  =ibarr;
    ztmp  =zsmeam;
    
    z0meam=zbar(ibar,ztmp,lat,t1,t2,t3);
//    INFO_Printf("ztmp=%g z0meam=%g\n",ztmp,z0meam);
    
    xmjc=alpha*(r/re-1.0); /* Eq.(14b) */
//    INFO_Printf("xmjc = %g\n",xmjc);

    /* These lines are in MDCASK but not in Baskes's code */
//    if(xmjc<=xzbl)
//    {
//        ret=zbl(r,charge,charge);
//        //INFO_Printf("MEAMFrame::phiid = %e: xmjc = %f < xzbl = %f\n",ret,xmjc,xzbl);
//        return ret;
//    }
//    else if(xmjc>=xzbl0)
//    {
//        phizbl=0.0;
//        frac=1.0;
//    }
//    else
//    {
//        frac=xcut((xmjc-xzbl0)/(xzbl-xzbl0));
//        phizbl=zbl(r,charge,charge);
//        //INFO_Printf("MEAMFrame::phiid: phizbl = %g  frac = %g  xmjc=%f in (%f,%f)\n",phizbl,xmjc,xzbl,xzbl0);
//    }
    
    /* These lines are in Baskes's code but not in MDCASK */
    rh0=rhof(r,beta0,re,1.0); /* Eq.(8a) */
    rh1=rhof(r,beta1,re,1.0);
    rh2=rhof(r,beta2,re,1.0);
    rh3=rhof(r,beta3,re,1.0);
    if(xmjc<=xzbl0)
    {
        arg=t1*SQR(rh1)+(2.0/3.0)*t2*SQR(rh2)+t3*(1-legend)*SQR(rh3);
        zdimer=1;
        roz=bar(zdimer*rh0,arg,ibar,zdimer,&dum1,&dum2)/z0meam;
        fdimer=frhoi(roz,asub,esub);

        phizbl=zbl(r,charge,charge);

        if(xmjc<xzbl)
            frac=0.0;
        else /* xzbl <= xjmc < xzbl0, default xzbl=-3, xzbl0=-1 */
            frac=xcut((xmjc-xzbl0)/(xzbl-xzbl0));
    }
    else
    {
        frac=1;
        phizbl=0;
        fdimer=0;
    }
    
    /* Baskes's code has a statement here
     *  if(nn) then ... endif
     *  don't know its meaning and purpose */

    

    //INFO_Printf("rh0=%g rh1=%g rh2=%g rh3=%g\n",rh0,rh1,rh2,rh3);
    
    /* need to calculate rh1, rh2 */
    
    

    /* in Baskes's code, rh1, rh2, rh3 is computed with nrho=2,4,6, if irho==1 */
    
    /* compute roz, which is rho_i^0(R) / Z_i
     * rho_i^0(R) is the background electron density for the reference structure
     *  for atom i and R is the nearest-neighbor distance
     * rho_i^0(R) is computed using the same function bar() which computes
     *  partial electron density contribution
     */
    rh0=rhof(r,beta0,re,1.0); /* Eq.(8a) */
//    INFO_Printf("r=%g beta0=%g re=%g rh0=%g\n\n",r,beta0,re,rh0);
    /* Eq.(8b,8c,8d,9) */    
    if(strcmp(lat,"dia")==0)
    {
        rh3=rhof(r,beta3,re,1.0);     /* in baskes, rh3 is computed with nrho = 6 */
//        INFO_Printf("rh0=%g rh3=%g\n",rh0,rh3);
        arg=(32.0/9.0)*t3*SQR(rh3);   /* arg is Gamma = sum_{j=1}^3 tj*(rho_i^(j))^2 */
        rh0=bar(z0meam*rh0,arg,ibar,z0meam,&dum1,&dum2);
        roz=rh0/z0meam;
//        INFO_Printf("arg=%g rh0=%g ibar=%d z0meam=%g\n",arg,rh0,ibar,z0meam);
    }
    else if(strcmp(lat,"hcp")==0)
    {
        rh3=rhof(r,beta3,re,1.0);
        arg=(1.0/3.0)*ts[3]*SQR(rh3);
        rh0=bar(z0meam*rh0,arg,ibar,z0meam,&dum1,&dum2);
        roz=rh0/z0meam;
    }
    else if(strcmp(lat,"dim")==0)
    {
        rh1=rhof(r,beta1,re,1.0);
        rh2=rhof(r,beta2,re,1.0);
        rh3=rhof(r,beta3,re,1.0);
        /* need to calculate rh1, rh2 */
        arg=t1*SQR(rh1)+(2.0/3.0)*t2*SQR(rh2)+t3*(1-legend)*SQR(rh3);
        rh0=bar(z0meam*rh0,arg,ibar,z0meam,&dum1,&dum2);
        roz=rh0/z0meam;
    }
    else
    {
        FATAL("lat = "<<lat<<" not recognized");
    }

    /* Eq.(4), pair interaction */
    phi=(2.0/ztmp)*( erose(r,re,alpha,esub,repul,attra)
                     - frhoi(roz,asub,esub) );

    /* Kuo, Clancy, MSMSE, 13, 1309 (2005), Eq.(15) */
    /* This line is not in Baskes's code, but is in MDCASK */
    //if(r<attra) phi+=repul*SQR(r-attra);
    
//    INFO_Printf("r = %g frac = %g phi = %g phizbl = %g\n",r,frac,phi,phizbl);

    /* this line is the same as MDCASK but is different from Baskes's code */
    //phi=frac*phi+(1.0-frac)*phizbl;     /* frac depends on xzbl */

    if(!enable_zbl_fdimer) fdimer=0.0;   /* LAMMPS code correspond to fdimer=0.0 */
//    INFO_Printf("r=%20.14e zbl=%20.14e\n",r,zbl(r,charge,charge));

    /* this line is the same Baskes's code */
    phi=frac*phi+(1.0-frac)*(phizbl-2*fdimer);

//    INFO_Printf("phi = %g ******** (after mix) fdimer = %g\n",phi,fdimer);
//    INFO_Printf("roz = %g asub = %g esub = %g\n",roz,asub,esub);
//    INFO_Printf("erose = %g frhoi = %g\n",erose(r,re,alpha,esub,repul,attra),
//                frhoi(roz,asub,esub));
//    INFO_Printf("xmjc = %g alpha = %g re = %g\n",xmjc,alpha,re);
    
    //if(phi>10000) FATAL("phi too big!");

    return phi;
}

double MEAMFrame::zbl(double r,double z1val,double z2val)
{
    int i; double z1, z2;
    double c[4]={0.028171,0.28022,0.50986,0.18175};
    double d[4]={0.20162,0.40290,0.94229,3.1998};
    double azero = 0.4685; /* azero = (9pi^2/128)^1/3 (0.529) Angstroms */
    double cc = 14.3997;
    double x, a, ret;

    //INFO_Printf("MEAMFrame::zbl r=%g charge=%g,%g\n",r,z1val,z2val);

    z1 = round(z1val);
    z2 = round(z2val);
    
    a=azero/(pow(z1,0.23)+pow(z2,0.23));
    ret=0.0;

    x = r/a;
    for(i=0;i<4;i++)
        ret += c[i]*exp(-d[i]*x);
            
    ret *= (z1*z2/r*cc);

    return ret;    
}

double MEAMFrame::bar(double rho0,double A,int ibar,double z,double *dang1,double *dang2)
{ /* calculates rhobar, Eq.(9a) as well as the derivatives wrt rho0 and A (angular density) */
    int ib;
    double ang, ex, x;
    double x0=-0.99; int n=99;
    
    if(rho0<=0) return 1e-20;

    ib = abs(ibar);
    if((ib==0)||(ib==4))
    { /* ang = sqrt(1+A/SQR(rho0)) */   /* ibarr = 0 for metals, e.g. Au, Kuo, Clancy, MSMSE, 13, 1309 (2005). */
        ang = 1+A/SQR(rho0);            /* this is Eq.(9) in the original MEAM potential for Si */
        /* Baskes differs from MDCASK */
//        if(ang<0)
//        {
//            *dang1=0.0;
//            *dang2=0.0;
//            return 1e-20;
//        }
        if(ang<1+x0)
        {
            x=ang-1;
            ang=sqrt((1+x0)*pow(x0/x,n));
            *dang1=ang*n/rho0;
            *dang2=-ang*n/(2*A);
        }
        else
        {
            ang=sqrt(ang);
            *dang1 = -A/(rho0*rho0*rho0*ang);
            *dang2 = 1/(2*rho0*rho0*ang);
        }
    }
    else if(ib==1)
    { /* ang = exp(0.5 A / rho0^2) */
        ang = exp(0.5*A/SQR(rho0));
        *dang1 = -ang*A/(rho0*rho0*rho0);
        *dang2 = ang/(2*rho0*rho0);
    }
    else if(ib==2)
    { /* ang = exp(0.5 A/z^2) */
        ang = exp(0.5*A/SQR(z));
        *dang1 = 0;
        *dang2 = ang/(2.0*SQR(z));
    }
    else if(ib==3)
    { /* ang = 2/(1+exp(-A/rho0^2) */
        ex = exp(-A/SQR(rho0));
        ang = 2.0/(1.0 + ex);
        *dang1 = -SQR(ang)*ex*A/(rho0*rho0*rho0);
        *dang2 = SQR(ang)/2.0*ex/SQR(rho0);
    }               
    return rho0*ang;
}

inline double MEAMFrame::erose(double r,double re,double alpha,double esub,double repuls,double attrac)
{
    double x, an3, E;
//    double f;
    
    /* Eqns. (14a) */

    x=alpha*(r/re-1.0);
    if(x>=0)
        an3=attrac;
    else
        an3=repuls;

    E=-esub*(1.+x+an3*x*x*x/(r/re))*exp(-x);

    /* the following lines are commented out in Baskes's code
     * but they appear in MDCASK */
    
//    if(x>=xzbl0)
//        E=-esub*(1.+x+an3*x*x*x/(r/re))*exp(-x);
//    else if(x<=xzbl)
//        E=0.0;
//    else
//    {
//        f=xcut((x-xzbl0)/(xzbl-xzbl0));
//        E=-f*esub*(1.+x+an3*x*x*x/(r/re))*exp(-x);
//    }
    
//    INFO_Printf("erose: x=%g f=%g E=%g\n",x,f,E);
    return E;
}

//double MEAMFrame::xcut(double x)
//{
//    return SQR(1.0 - SQR(SQR(x)));
//}

Matrix33 MEAMFrame::dcmij(Vector3 rr,double rs)
{
    int i, j; double r1, r13, tmp1;
    Matrix33 drcos;

    r1 = 1.0 / rs;
    r13 = r1 * r1 * r1;
    for(i=0;i<3;i++)
    {
        tmp1 = r13 * rr[i];
        for(j=0;j<3;j++)
        {
            drcos[j][i] = -tmp1 * rr[j];
        }
        drcos[i][i] += r1;
    }
    return drcos;
}

inline double MEAMFrame::zbar(int ibar,double z,char *lattice,double t1,double t2,double t3)
{
    if(ibar<=0)
        return z;
    else
    {
        FATAL("ibar = "<<ibar<<" not implemented!");
        return 0;
    }
}

inline double MEAMFrame::rhof(double r,double abc,double re,double rozero)
{
    /* there is a generalization of this function
     * with a pre-factor r^nrho in Baskes's code
     */
    return rozero*exp(-abc*(r/re-1.0));  /* Eqs.(16a) and (16b) */
}

inline double MEAMFrame::frhoi(double rhotp,double asub,double esub)
{ /* Embedding function */
    if(rhotp>0)
        return asub*esub*rhotp*log(rhotp); /* Eq.(15) */
    else
        return asub*esub*rhotp*(-100.0);
}

inline double MEAMFrame::dfrhoi(double rho,double asub,double esub)
{ /* Embedding function */
    if(rho>0)
        return asub*esub*(log(rho)+1.0); /* Eq.(15) */
    else
        return asub*esub*(-100.0);
}

inline double MEAMFrame::rscrn(double r2ij)
{ /* r cut-off function */
    double fac, rij, frcut=0.9;

    if(enable_square_rscrn)
    {
        //INFO_Printf("modified rscrn using rij squared\n");
        fac=(r2ij-frcut*rcutmeam2)/((1.0-frcut)*rcutmeam2);
    }
    else
    {
        //INFO_Printf("original rscrn\n");
        rij = sqrt(r2ij);
        fac = 1.0-(rcutmeam-rij)/(1.0-frcut);
    }            
    if(fac <= 0.0)
        return 1.0;
    else if(fac >= 1.0)
        return 0.0;
    else
        return xcut(fac);
}
            
double MEAMFrame::dscrn(double rij, double rik, double rjk, double cmin, double cmax)
{
    double dotil, dotjl, fac;
    double argb, argt, arg;

    if((rik>sconst*rij)||(rjk>sconst*rij)) return 1.0;

    /* the following is in Baskes's code but not in MDCASK */
    dotil = rik+rij-rjk;
    if(dotil<=2.0e-3) return 1.0;
    dotjl = rij+rjk-rik;
    if(dotjl<=2.0e-3) return 1.0;
    argb=SQR(rij)-SQR(rik-rjk);
    argt=4.0*(rik*rij)-SQR(rik-rjk+rij);
    arg=argt/argb;
    //INFO_Printf("argt=%g argb=%g arg = %g cmin=%g cmax=%g\n",argt,argb,arg,cmin,cmax);
    if(arg>=cmax) return 1.0;
    else if(arg<=cmin) return 0.0;
    else
    {
        fac = (cmax-arg)/(cmax-cmin);
        //INFO_Printf("argt=%g argb=%g arg = %g fac=%g xcut=%g\n",argt,argb,arg,fac,xcut(fac));        
        return xcut(fac);
    }

#if 0 /* MDCASK code */    
    argb = 1.0 - SQR((rik-rjk)/rij);
    if(argb==0) return 1.0;

    argt = 4.0*(rik/rij)-SQR((rik-rjk+rij)/rij);
    arg=argt/argb;

    /* arg is
     * C = [2(X_ij+X_jk) - (X_ij-X_jk)^2 - 1] / [1 - (X_ij-X_jk)^2]
     * X_ij = (rij/rik)^2 X_jk = (rjk/rik)^2
     */

    /* M.I.Baskes, J.E.Angelo, C.L.Bisson, Modell.Sim.Mater.Sci.Eng. 2, 505 (1994)
     *
     * S_ijk = 0                            if C<=Cmin
     *         exp(-[(Cmax-C)/(C-Cmin)]^2)  if Cmin<C<Cmax
     *         1                            if C>=Cmax
     *
     * However, here it uses a different form for S_ijk (polynomial)
     *
     * See M.I.Baskes, Mater.Chem.Phys. 50, 152 (1997)
     *     Kuo, Clancy, MSMSE, 13, 1309 (2005), Eq.(15)
     */
    
    if((arg>cmax)||(arg<0.0))
        return 1.0;
    else if(arg<cmin)
        return 0.0;
    else
        return xcut( (cmax-arg)/(cmax-cmin) );
#endif
}





void MEAMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int i, size;
    size=_NP*allocmultiple;

    /* MEAM potential storage space */
    Realloc(atpe2b,double,size);
    Realloc(atpe3b,double,size);
    Realloc(rhotot,double,size);
    Realloc(embf,double,size);
    Realloc(c8a,double,size);
    Realloc(dang1,double,size);
    Realloc(dang2,double,size);
    Realloc(cg8c,double,size);
    Realloc(tav,Vector3,size);
    Realloc(rhsq,Vector3,size);
    Realloc(a8b,Vector3,size);
    Realloc(ag,Vector3,size);
    Realloc(b8c,Matrix33,size);
    Realloc(d8d,Matrix333,size);

    /* Allocate 2-dimensional array for S_ij */
    Realloc(scrnab_mem,char,size*_NNM*sizeof(double)+size*sizeof(double *));
    scrnab=(double **)(scrnab_mem+size*_NNM*sizeof(double));
    for(i=0;i<size;i++) scrnab[i]=(double *)(scrnab_mem+i*_NNM*sizeof(double));

    Realloc(dscrnab_mem,char,size*_NNM*3*sizeof(double)+size*sizeof(double *));
    dscrnab=(double **)(dscrnab_mem+size*_NNM*3*sizeof(double));
    for(i=0;i<size;i++) dscrnab[i]=(double *)(dscrnab_mem+i*_NNM*3*sizeof(double));
    

    bindvar("atpe2b",atpe2b,DOUBLE);
    bindvar("atpe3b",atpe3b,DOUBLE);
    bindvar("rhotot",rhotot,DOUBLE);
    bindvar("embf",embf,DOUBLE);
    bindvar("c8a",c8a,DOUBLE);
    bindvar("dang1",dang1,DOUBLE);
    bindvar("dang2",dang2,DOUBLE);
    bindvar("cg8c",cg8c,DOUBLE);
    bindvar("tav",tav,DOUBLE);
    bindvar("rhsq",rhsq,DOUBLE);
    bindvar("a8b",a8b,DOUBLE);
    bindvar("ag",ag,DOUBLE);
    bindvar("b8c",b8c,DOUBLE);
    bindvar("d8d",d8d,DOUBLE);

}    

void MEAMFrame::potential()
{
    refreshneighborlist();    
    screen();               /* compute all screen functions S_ij */
    dscreen();              /* derivatives of the screen functions */
    rhoMEAM();              /* electron density, pair potential, and total potential energy */
    kraMEAM();              /* forces on atoms */
    if(noscr==0) dscrfor(); /* force due to screening function */
}

void MEAMFrame::initvars()
{
    MDPARALLELFrame::initvars();

    betas[0]=betas[1]=betas[2]=betas[3]=0;
    ts[0]=ts[1]=ts[2]=ts[3]=0;
    strcpy(potfile,"meamdata");

    strcpy(lat,"unknown");
    strcpy(elt,"unknown");
}
void MEAMFrame::initparser()
{
    MDPARALLELFrame::initparser();

    bindvar("xzbl",&xzbl,DOUBLE);
    bindvar("xzbl0",&xzbl0,DOUBLE);
    bindvar("enable_zbl_fdimer",&enable_zbl_fdimer,INT);
    bindvar("enable_square_rscrn",&enable_square_rscrn,INT);
}

int MEAMFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readMEAM",readMEAM());
    bindcommand(name,"printpairpot",printpairpot());

#ifdef _PARALLEL
    bindcommand(name,"Broadcast_MEAM_Param",Broadcast_MEAM_Param());
#endif
    return -1;
}


void MEAMFrame::printpairpot()
{    
    int elti, eltj;
    double rmin, dr, rmax, r, phi, phip;
    char pairfile[200];
    FILE *fp;
    
    elti = (int) input[0];
    eltj = (int) input[1];
    rmin = input[2];
    dr   = input[3];
    rmax = input[4];

    strcpy(pairfile,"pairpot_baskes.dat");
    fp=fopen(pairfile,"w");
    if(fp==NULL)
    {
        FATAL("printpairpot: file "<<pairfile<<" open failure.");
    }
    if(elti!=eltj)
    {
        FATAL("elti ("<<elti<<") != eltj("<<eltj<<") not implemented");
    }
    if((rmax<rmin)||(dr<=0))
        FATAL("rmax cannot be smaller than rmin, dr must be positive");
    
    phip = 0;
    for(r=rmin;r<=rmax;r+=dr)
    {
        phi = phiid(r,elti);
        fprintf(fp,"%21.14e  %21.14e %21.14e\n",r,phi,phip);
    }
    fclose(fp);
    INFO("results written to "<<pairfile);
}


#ifdef _PARALLEL
void MEAMFrame::Broadcast_MEAM_Param()
{
    double *buf_double;
    int nparam;

    /* master asking slaves to call the same function */
    if(myDomain==0) Master_to_Slave("Broadcast_MEAM_Param");    

    /* packing parameters */
    nparam = 36;
    buf_double = (double *) malloc(nparam*sizeof(double));
    if(myDomain==0)
    {
    buf_double[0]  = zsmeam;
    buf_double[1]  = alphas;
    buf_double[2]  = betas[0];
    buf_double[3]  = betas[1];
    buf_double[4]  = betas[2];
    buf_double[5]  = betas[3];
    buf_double[6]  = esubs;
    buf_double[7]  = asubs;
    buf_double[8]  = ts[0];
    buf_double[9]  = ts[1];
    buf_double[10] = ts[2];
    buf_double[11] = ts[3];
    buf_double[12] = rozros;
    buf_double[13] = rcutmeam;
    buf_double[14] = cmin;
    buf_double[15] = cmax;
    buf_double[16] = repuls;
    buf_double[17] = attrac;
    buf_double[18] = legend;
    buf_double[19] = (double) ibarr;
    buf_double[20] = (double) noscr;
    buf_double[21] = res;
    buf_double[22] = alat;
    buf_double[23] = ielement;
    buf_double[24] = rcutmeam2;
    buf_double[25] = sconst;
    buf_double[26] = scnres;
    buf_double[27] = xzbl;
    buf_double[28] = xzbl0;
    buf_double[29] = hmeam;
    buf_double[30] = _RLIST;
    buf_double[31] = _SKIN;
    buf_double[32] = rcutmeam2;
    buf_double[33] = res;
    buf_double[34] = _NNM;
    buf_double[35] = _NIC;
    }
    /* broadcasting */
    MPI_Bcast(buf_double,nparam,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(myDomain!=0)
    {
    /* unpacking parameters
     * the following lines can be generated from the above lines by
     * replace regexp:
     *  \(buf_double\[[0-9]+\]\)[ ]+= \([^;]+\);
     *  \2 = \1;
    */
    zsmeam = buf_double[0];
    alphas = buf_double[1];
    betas[0] = buf_double[2];
    betas[1] = buf_double[3];
    betas[2] = buf_double[4];
    betas[3] = buf_double[5];
    esubs = buf_double[6];
    asubs = buf_double[7];
    ts[0] = buf_double[8];
    ts[1] = buf_double[9];
    ts[2] = buf_double[10];
    ts[3] = buf_double[11];
    rozros = buf_double[12];
    rcutmeam = buf_double[13];
    cmin = buf_double[14];
    cmax = buf_double[15];
    repuls = buf_double[16];
    attrac = buf_double[17];
    legend = buf_double[18];
    ibarr = (int) buf_double[19];
    noscr = (int) buf_double[20];
    res = buf_double[21];
    alat = buf_double[22];
    ielement = buf_double[23];
    rcutmeam2 = buf_double[24];
    sconst = buf_double[25];
    scnres = buf_double[26];
    xzbl = buf_double[27];
    xzbl0 = buf_double[28];
    hmeam = buf_double[29];
    _RLIST = buf_double[30];
    _SKIN = buf_double[31];
    rcutmeam2 = buf_double[32];
    res = buf_double[33];
    _NNM = (int) buf_double[34];
    _NIC = (int) buf_double[35];
    }
    
    free(buf_double);

    /* broadcasting string parameters */
    MPI_Bcast(elt,10,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(lat,10,MPI_CHAR,0,MPI_COMM_WORLD);
}
#endif//_PARALLEL

#ifdef _TEST

/* Main Program Begins */
class MEAMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

























































