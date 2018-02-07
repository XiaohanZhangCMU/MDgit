/*
  siedip.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Apr 22 16:40:56 2010

  FUNCTION  :  MD simulation package of Si using EDIP potential
*/

#include "siedip.h"

void SIEDIPFrame::initvars()
{
    A = 5.6714030;
    B = 2.0002804;
    rh = 1.2085196;
    a = 3.1213820;
    sig = 0.5774108;
    lam = 1.4533108;
    gam = 1.1247945;
    b = 3.1213820;
    c = 2.5609104;
    delta = 78.7590539;
    mu = 0.6966326;
    Qo = 312.1341346;
    palp = 1.4074424;
    bet = 0.0070975;
    alp = 3.1083847;
    
    u1 = -0.165799;
    u2 = 32.557;
    u3 = 0.286198;
    u4 = 0.66; 
    
    INFO_Printf("\n EDIP-Si Parameters: \n");
    INFO_Printf("%e %e %e %e %e \n",A,B,rh,a,sig);
    INFO_Printf("%e %e %e %e %e \n",lam,gam,b,c,delta);
    INFO_Printf("%e %e %e %e %e \n",mu,Qo,palp,bet,alp);
    bg=a;
    eta = delta/Qo;
    
    _RLIST=a*1.1;
    _SKIN=_RLIST-a;

    strcpy(incnfile,"../siedip.cn");
    DUMP(HIG"SIEDIPFrame initvars"NOR);
    MDFrame::initvars();
}

void SIEDIPFrame::edip()
{
    /*------------------------- VARIABLE DECLARATIONS -------------------------*/
    
    int i,j;
    double dx,dy,dz,r,rsqr,asqr;
    double rinv,rmainv,xinv,xinv3,den,Z,fZ;
    double dV2j,dV2ijx,dV2ijy,dV2ijz,pZ,dp;
    double temp0,temp1;
    double Qort,muhalf,u5;
    double rmbinv,winv,dwinv,tau,dtau,lcos,x,H,dHdx,dhdl;
    double dV3rij,dV3rijx,dV3rijy,dV3rijz;
    double dV3rik,dV3rikx,dV3riky,dV3rikz;
    double dV3l,dV3ljx,dV3ljy,dV3ljz,dV3lkx,dV3lky,dV3lkz;
    double dV2dZ,dxdZ,dV3dZ;
    double dEdrl,dEdrlx,dEdrly,dEdrlz;
    double bmc,cmbinv;
    double fjx,fjy,fjz,fkx,fky,fkz;
    
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    typedef struct{      
        double t0,t1,t2,t3;  /* various V2 functions and derivatives */
        double dx,dy,dz;     /* unit separation vector */
        double r;		 /* bond length (only needed for virial) */
    } store2;
    store2 s2[NNM]; /* two-body interaction quantities, r<a */  
    int n2;                /* size of s2[] */
    int num2[NNM];  /* atom ID numbers for s2[] */
    
    typedef struct{
        double g,dg;         /* 3-body radial function and its derivative */
        double rinv;         /* 1/r */
        double dx,dy,dz;     /* unit separation vector */
        double r;		 /* bond length (only needed for virial) */
    } store3;
    store3 s3[NNM]; /* three-body interaction quantities, r<b */
    int n3;                /* size of s3[] */
    int num3[NNM];  /* atom ID numbers for s3[] */
    
    typedef struct{
        double df;           /* derivative of neighbor function f'(r) */
        double sum;		 /* array to accumulate coordination force prefactors */
        double dx,dy,dz;     /* unit separation vector */
        double r;		 /* bond length (only needed for virial) */
    } storez;
    storez sz[NNM]; /* coordination number stuff, c<r<b */
    int nz;                /* size of sz[] */
    int numz[NNM];  /* atom ID numbers for sz[] */
    
    int nj,nk,nl;         /* indices for the store arrays */

    int n0, n1;
    int ipt, jpt, kpt, lpt;
    double V2, V3, virial; /* virial is the diagonal part of _VIRIAL */
    double eterm2, eterm3;
    Vector3 sij, rij, fgj, rgj, fgk, rgk;
    
    DUMP("EDIP");
    
    refreshneighborlist();

//    n0=(int)rint(MP_cpu()*_NP/MP_ncpu());
//    n1=(int)rint((MP_cpu()+1)*_NP/MP_ncpu());
    n0=0;
    n1=_NP;

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    for(i=n0;i<n1;i++)
    { _F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    _EPOT=V2=V3=0; /* _EPOT = V2 + V3 */

    /* COMBINE COEFFICIENTS */
    
    asqr = a*a;
    Qort = sqrt(Qo);
    muhalf = mu*0.5;
    u5 = u2*u4;    
    bmc = b-c;
    cmbinv = 1.0/(c-b);

    /* strange, why needed? */
    dx=dy=dz=0;
    
    /*--- LEVEL 1: OUTER LOOP OVER ATOMS ---*/
    for(ipt=n0; ipt<n1; ipt++)
    {   /* not parallelized yet */
        /* RESET COORDINATION AND NEIGHBOR NUMBERS */

        Z = 0.0;
        n2 = 0;  
        n3 = 0;
        nz = 0;

        /*--- LEVEL 2: LOOP PREPASS OVER PAIRS ---*/

        for(j=0; j<nn[ipt]; j++)
        {
            jpt=nindex[ipt][j];
            sij=_SR[jpt]-_SR[ipt];
            sij.subint();
            rij=_H*sij;
            dx=rij.x; dy=rij.y; dz=rij.z;

            /* TEST IF WITHIN OUTER CUTOFF */
            rsqr=rij.norm2();
            if(rsqr>=asqr) continue;
            r=sqrt(rsqr);
            
            /* PARTS OF TWO-BODY INTERACTION r<a */
 
            num2[n2] = jpt;
            rinv = 1.0/r;
            dx *= rinv;
            dy *= rinv;
            dz *= rinv;
            rmainv = 1.0/(r-a);
            s2[n2].t0 = A*exp(sig*rmainv);
            s2[n2].t1 = pow(B*rinv,rh);
            s2[n2].t2 = rh*rinv;
            s2[n2].t3 = sig*rmainv*rmainv;
            s2[n2].dx = dx;
            s2[n2].dy = dy;
            s2[n2].dz = dz;
            s2[n2].r = r;
            n2++;

            /* RADIAL PARTS OF THREE-BODY INTERACTION r<b */
            
            if(r < bg) {
                
                num3[n3] = jpt;
                rmbinv = 1.0/(r-bg);
                temp1 = gam*rmbinv;
                temp0 = exp(temp1);
#if V3g_on
                s3[n3].g = temp0;
                s3[n3].dg = -rmbinv*temp1*temp0;
#else
                s3[n3].g = 1;
                s3[n3].dg = 0;
#endif
                s3[n3].dx = dx;
                s3[n3].dy = dy;
                s3[n3].dz = dz;
                s3[n3].rinv = rinv;
                s3[n3].r = r;
                n3++;
                
                /* COORDINATION AND NEIGHBOR FUNCTION c<r<b */

                if(r<b) {
                    if(r < c) 
                        Z += 1.0;
                    else {
                        xinv = bmc/(r-c);
                        xinv3 = xinv*xinv*xinv;
                        den = 1.0/(1 - xinv3);
                        temp1 = alp*den;
                        fZ = exp(temp1);
                        Z += fZ;
                        numz[nz] = jpt;
                        sz[nz].df = fZ*temp1*den*3.0*xinv3*xinv*cmbinv;   /* df/dr */
                        sz[nz].dx = dx;
                        sz[nz].dy = dy;
                        sz[nz].dz = dz;
                        sz[nz].r = r;
                        nz++;
                    }      
                }
            }
        }
        //local[ipt].z=Z;

        //if(measure) coord_total += Z;

        /* ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES */

        for(nl=0; nl<nz; nl++) sz[nl].sum=0.0;

        /* ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION */

#if V2_on
        
#if V2Z_on
        temp0 = bet*Z;
        pZ = palp*exp(-temp0*Z);         /* bond order */
        dp = -2.0*temp0*pZ;         /* derivative of bond order */
#else
        pZ = palp*exp(-bet*16);
        dp=0.0;
#endif

        /*--- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---*/

        for(nj=0; nj<n2; nj++)
        {

            temp0 = s2[nj].t1 - pZ;

            /* two-body energy V2(rij,Z) */
            eterm2 = temp0*s2[nj].t0; 
            V2 += eterm2;
            jpt = num2[nj];

            _EPOT_IND[ipt]+=eterm2*0.5;
            _EPOT_IND[jpt]+=eterm2*0.5;
            
            /* two-body forces */

            dV2j = - (s2[nj].t0) * ((s2[nj].t1)*(s2[nj].t2) 
                                    + temp0 * (s2[nj].t3));   /* dV2/dr */
            dV2ijx = dV2j * s2[nj].dx;
            dV2ijy = dV2j * s2[nj].dy;
            dV2ijz = dV2j * s2[nj].dz;
            fgj.set(dV2ijx, dV2ijy, dV2ijz);
            _F[ipt] += fgj;
            _F[jpt] -= fgj;
            
            /* dV2/dr contribution to virial */
            virial -= s2[nj].r * (dV2ijx*s2[nj].dx 
                              + dV2ijy*s2[nj].dy + dV2ijz*s2[nj].dz);

            rgj.set(s2[nj].dx, s2[nj].dy, s2[nj].dz);
            _VIRIAL.addnvv(-s2[nj].r,fgj,rgj);
            
            /*--- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---*/
            
            dV2dZ = - dp * s2[nj].t0;
            for(nl=0; nl<nz; nl++) sz[nl].sum += dV2dZ;
            
        }
#endif

        /* COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION */
        
#if V3_on
        
#if V3Z_on
        winv = Qort*exp(-muhalf*Z); /* inverse width of angular function */    
        dwinv = -muhalf*winv;       /* its derivative */
        temp0 = exp(-u4*Z);
        tau = u1+u2*temp0*(u3-temp0); /* -cosine of angular minimum */
        dtau = u5*temp0*(2*temp0-u3); /* its derivative */
#else
        winv =  Qort*exp(-muhalf*4);
        dwinv = 0.0;
        tau = 1.0/3.0;
        dtau=0.0;
#endif

        /*--- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---*/
        
        for(nj=0; nj<(n3-1); nj++)
        {
            jpt=num3[nj];
            
            /*--- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---*/
            for(nk=nj+1; nk<n3; nk++)
            {
                kpt=num3[nk];

                /* angular function h(l,Z) */
                
                lcos = s3[nj].dx * s3[nk].dx + s3[nj].dy * s3[nk].dy
                    + s3[nj].dz * s3[nk].dz;
                x = (lcos + tau)*winv;
                temp0 = exp(-x*x);
#if V3h_on
                H = lam*(1 - temp0 + eta*x*x);
                dHdx = 2*lam*x*(temp0 + eta);
                dhdl = dHdx*winv;
#else
                H=1.0;
                dhdl=0.0;
#endif

                /* three-body energy */
                
                temp1 = s3[nj].g * s3[nk].g;
                eterm3 = temp1*H;
                
                V3 += eterm3;
                _EPOT_IND[ipt] += eterm3/3;
                _EPOT_IND[jpt] += eterm3/3;
                _EPOT_IND[kpt] += eterm3/3;
                
                /* (-) radial force on atom j */
                dV3rij = s3[nj].dg * s3[nk].g * H;
                dV3rijx = dV3rij * s3[nj].dx;
                dV3rijy = dV3rij * s3[nj].dy;
                dV3rijz = dV3rij * s3[nj].dz;
                fjx = dV3rijx;
                fjy = dV3rijy;
                fjz = dV3rijz;
                
                /* (-) radial force on atom k */
                
                dV3rik = s3[nj].g * s3[nk].dg * H;
                dV3rikx = dV3rik * s3[nk].dx;
                dV3riky = dV3rik * s3[nk].dy;
                dV3rikz = dV3rik * s3[nk].dz;
                fkx = dV3rikx;
                fky = dV3riky;
                fkz = dV3rikz;
                
                /* (-) angular force on j */
                
                dV3l = temp1*dhdl;
                dV3ljx = dV3l * (s3[nk].dx - lcos * s3[nj].dx) * s3[nj].rinv;
                dV3ljy = dV3l * (s3[nk].dy - lcos * s3[nj].dy) * s3[nj].rinv;
                dV3ljz = dV3l * (s3[nk].dz - lcos * s3[nj].dz) * s3[nj].rinv;
                fjx += dV3ljx;
                fjy += dV3ljy;
                fjz += dV3ljz;
                
                /* (-) angular force on k */
                
                dV3lkx = dV3l * (s3[nj].dx - lcos * s3[nk].dx) * s3[nk].rinv;
                dV3lky = dV3l * (s3[nj].dy - lcos * s3[nk].dy) * s3[nk].rinv;
                dV3lkz = dV3l * (s3[nj].dz - lcos * s3[nk].dz) * s3[nk].rinv;
                fkx += dV3lkx;
                fky += dV3lky;
                fkz += dV3lkz;
                
                /* apply radial + angular forces to i, j, k */

                fgj.set(fjx, fjy, fjz);
                fgk.set(fkx, fky, fkz);

                _F[jpt] -= fgj;
                _F[kpt] -= fgk;
                _F[ipt] += fgj; _F[ipt] += fgk;

                /* dV3/dR contributions to virial */
                
                virial -= s3[nj].r * (fjx*s3[nj].dx 
                                      + fjy*s3[nj].dy + fjz*s3[nj].dz); 
                virial -= s3[nk].r * (fkx*s3[nk].dx 
                                      + fky*s3[nk].dy + fkz*s3[nk].dz);
                
                rgj.set(s3[nj].dx,s3[nj].dy,s3[nj].dz); 
                rgk.set(s3[nk].dx,s3[nk].dy,s3[nk].dz);
                _VIRIAL.addnvv(-s3[nj].r,fgj,rgj);
                _VIRIAL.addnvv(-s3[nk].r,fgk,rgk);
                
                /* prefactor for 4-body forces from coordination */
#if V3Z_on
                dxdZ = dwinv*(lcos + tau) + winv*dtau;
                dV3dZ = temp1*dHdx*dxdZ;
                
                /*--- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---*/
                
                for(nl=0; nl<nz; nl++) sz[nl].sum += dV3dZ; 
#endif
            }
        }
#endif

        /*--- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---*/
        
        for(nl=0; nl<nz; nl++)
        {
            dEdrl = sz[nl].sum * sz[nl].df;
            dEdrlx = dEdrl * sz[nl].dx;
            dEdrly = dEdrl * sz[nl].dy;
            dEdrlz = dEdrl * sz[nl].dz;
            _F[ipt].x += dEdrlx;
            _F[ipt].y += dEdrly;
            _F[ipt].z += dEdrlz;
            lpt = numz[nl];
            _F[lpt].x -= dEdrlx;
            _F[lpt].y -= dEdrly;
            _F[lpt].z -= dEdrlz;

            /* dE/dZ*dZ/dr contribution to virial */
            
            virial -= sz[nl].r * (dEdrlx*sz[nl].dx 
                    + dEdrly*sz[nl].dy + dEdrlz*sz[nl].dz); 
        }
    }

    _EPOT = V2 + V3;
    virial /= 3.0;

    /* zero out the forces on fixed atoms */
//    for(i=0;i<_NP;i++)
//    {
//        if(fixed[i])
//            _F[i].clear();
//    }

#if 0
    DUMP("EDIP Collect");
    /* _EPOT, _F, _EPOT_IND, _VIRIAL will be summed */
    MP_AccumulateTo(_shmEPOTptr,&(_EPOT),1);
    MP_AccumulateTo((double *)_shmF,(double *)_F,3*_NP);
    MP_AccumulateTo(_shmEPOT_IND,_EPOT_IND,_NP);
    MP_AccumulateTo((double *)((*_shmVIRIALptr)[0]),
                        (double *)(_VIRIAL[0]),9);
    
    StrMemCpy(&(_EPOT),_shmEPOTptr,sizeof(double));
    StrMemCpy((double *)_F,(double *)_shmF,3*_NP*sizeof(double));
    StrMemCpy(_EPOT_IND,_shmEPOT_IND,_NP*sizeof(double));
    StrMemCpy((double *)(_VIRIAL[0]),
              (double *)((*_shmVIRIALptr)[0]),9*sizeof(double));
    
    DUMP("EDIP complete");
#endif
}
    
void SIEDIPFrame::potential()
{
    edip();
}




#ifdef _TEST

/* Main Program Begins */
class SIEDIPFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

