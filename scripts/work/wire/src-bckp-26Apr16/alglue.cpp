/*
  alglue.cpp
  by Wei Cai  caiwei@mit.edu, Maurice de Koning maurice@mit.edu
  Last Modified : Mon Jan  1 21:25:34 2007

  FUNCTION  : MD++ implementation of Al Glue potential
              (Ercolessi-Adams Glue potential)
*/

#include "alglue.h"


int ALGLUEFrame::initglue()
{
    /*   Initialize linear grid points reduntantly on each node to save
     *   transfer time.
     */
    int i;
    double c, phi, dphi, d2phi;
    double rho, drho, d2rho, rsq, r;
    double glue, dglue, d2glue;


//    INFO("ALGLUEFrame::initglue");
    
    rcutoff=5.5580544182;
    rmin=rcutoff-5.0;
    rsqmin=rmin*rmin;
    deltarsq=(rcutoff*rcutoff-rsqmin)/(NGRID-1);
    invdeltarsq=1./deltarsq;
    coordmax=1.2;
    deltacoord=coordmax/(NGRID-1);
    invdeltacoord=1./deltacoord;

    for(i=0;i<NGRID;i++)
    {
        rsq=rsqmin+i*deltarsq;
        r=sqrt(rsq);
        v2(r,&phi,&dphi,&d2phi);
        rh(r,&rho,&drho,&d2rho);
        phitab[i]=phi;
        dphitab[i]=-dphi/r;
        rhotab[i]=rho;
        drhotab[i]=-drho/r;
    }
    
    for(i=0;i<NGRID;i++)
    {
        c=i*deltacoord;
        uu(c,&glue,&dglue,&d2glue);
        utab[i]=glue;
        dutab[i]=dglue;
    }
    _RLIST = 1.1 * rcutoff;
    _SKIN = _RLIST - rcutoff;

    return 0;
}

int ALGLUEFrame::al_glue()
{
    int i, j, k, jpt;
    Vector3 sij,rij,fij;
    double r2ij;
    double rcutoff_sq, rk, weight, rho, phi, dphi, drho;

    refreshneighborlist();
    
    rcutoff_sq = rcutoff*rcutoff;
    
    for(i=0;i<_NP;i++)
    {
        _F[i].clear(); _EPOT_IND[i]=0;
        coord[i]=0;
    }
    _EPOT=0;
    _VIRIAL.clear();
    
    for(j=0;j<_NP*_NNM;j++)
    {
        //INFO("j="<<j);
        //rijstore[j].clear();
    }

    /* compute rho's on each atom */
//    INFO_Printf("nn[0]=%d\n",nn[0]);
    /* do on i-particles */
    for(i=0;i<_NP;i++)
    {
        /* do on j-particles */
        for(j=0;j<nn[i];j++)
        {
            /* modify here for binary systems (0/1) */
            jpt=nindex[i][j];
            if(i>=jpt) continue;
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();

            //rijstore[i*_NNM+j]=rij;

            if(r2ij>=rcutoff_sq) continue;

            rk=(r2ij-rsqmin)*invdeltarsq; 
            k=(int)floor(rk); 
            if(k<0) k=0;
            weight=rk-k;
            rho=weight*rhotab[k+1]+(1-weight)*rhotab[k];
            coord[i]+=rho;
            coord[jpt]+=rho;
//            if(i==0) INFO_Printf("coord[0]=%e\n",coord[i]);
        }
    }

    for(i=0;i<_NP;i++)
    {
        rk=coord[i]*invdeltacoord; 
        k=(int)floor(rk); 
        if(k>=NGRID-1) k=NGRID-2;
        weight=rk-k;
        _EPOT_IND[i]=weight*utab[k+1]+(1-weight)*utab[k];
        _EPOT+=_EPOT_IND[i];
        deru[i]=weight*dutab[k+1]+(1-weight)*dutab[k];
//        if(i==0)
//            INFO_Printf("i=0 k=%d weight=%e utab[k]=%e invdeltarcoord=%e\n",
//                             k,weight,utab[k],invdeltacoord);
    }
//    INFO_Printf("EPOT_IND[0]=%e\n",_EPOT_IND[0]);

    /* compute forces and virial tensor */
    for(i=0;i<_NP;i++)
    {
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(i>=jpt) continue;

            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;            
            //rij=rijstore[i*_NNM+j];
            
            r2ij=rij.norm2();
            if(r2ij>=rcutoff_sq) continue;

            rk=(r2ij-rsqmin)*invdeltarsq;
            k=(int)floor(rk);
            if(k<0) k=0;
            weight = rk-k;
            phi=weight*phitab[k+1]+(1.0-weight)*phitab[k];
            _EPOT_IND[i]+=0.5*phi;
            _EPOT_IND[jpt]+=0.5*phi;
            _EPOT+=phi;
            dphi=weight*dphitab[k+1]+(1.0-weight)*dphitab[k];
            drho=weight*drhotab[k+1]+(1.0-weight)*drhotab[k];
            dphi+=(deru[i]+deru[jpt])*drho;
            
            fij=rij*dphi;
            _F[i]-=fij;
            _F[jpt]+=fij;

            if(!(fixed[i]||fixed[jpt]))
                _VIRIAL.addnvv(dphi,rij,rij);
            else if(!(fixed[i]&&fixed[jpt]))
                _VIRIAL.addnvv(0.5*dphi,rij,rij);
        }
    }
//    INFO_Printf("EPOT_IND[0]=%e\n",_EPOT_IND[0]);
    
//    INFO("initglue completed.");
    return 0;
}


void ALGLUEFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int size;
    size=_NP*allocmultiple;
    
    /* Shared Memory */
    Realloc(deru,double,size);
    Realloc(coord,double,size);
    //Realloc(rijstore,Vector3,size*_NNM);
}    

void ALGLUEFrame::potential()
{
    al_glue();
}

void ALGLUEFrame::initvars()
{
//    INFO("initvars");
    MDPARALLELFrame::initvars();
    initglue();
}

void ALGLUEFrame::initparser()
{
    MDPARALLELFrame::initparser();
}

int ALGLUEFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    return -1;
}
   
#ifdef _TEST

/* Main Program Begins */
class ALGLUEFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
