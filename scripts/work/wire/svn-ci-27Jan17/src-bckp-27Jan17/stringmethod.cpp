// Last Modified : Fri Apr  4 11:10:35 2008

/* an alternative implementation of stringrelax, may be generalized
 * to finite temperature string method in the future
 * not fully tested.  do not use */

#include "stringmethod.h"

void MDFrame::runstringmethod()
{/*
  * The string method
  * J. of Chemical Physics, vol.126, 164103 (2007)
  * Apr. 01 2008, Keonowook Kang
  * Must constrain whole atoms.
  */
    int i, j, k, m, n, ipt, size;
    Vector3 ds, ds0, dr, dt, dR, **_RcStore;
    Matrix33 hinv;
    double s, y_csp, dRnorm, dRnorm2, Xi0, XiN, s_j, hi, w, wbar, yi, yim1;
    double *Ec, *Lc;
    void *c;
    FILE *fp;

    INFO("The String Method");

    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL)
            ERROR("config1 is not set!");
        else
            ERROR("config2 is not set!");
        return;
    }

    if(_Rc==NULL || _Fc==NULL)
    {
        AllocChain();
        caldscom12();
        initRchain();
    }
    
    n=constrainatoms[0]; 
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    Lc=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    size=sizeof(double *)*(3*n) + sizeof(double)*(3*n)*(_CHAINLENGTH+1) + 10;

    size=sizeof(Vector3 *)*(_CHAINLENGTH+1)+sizeof(Vector3)*(_CHAINLENGTH+1)*n + 10;
    c=malloc(size); memset(c,0,size);
    _RcStore=(Vector3 **) c;   _RcStore[0]=(Vector3 *)(_RcStore+(_CHAINLENGTH+1));
    for(i=1;i<=_CHAINLENGTH;i++) _RcStore[i]=_RcStore[i-1]+n;

    fp=fopen("stringeng.out","w");

    /* compute current constrain value */
    constrain_dist=0; s=0;
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        constrain_dist+=ds0.norm2();
    }
    INFO_Printf("string method: constrain_dist=%e\n",constrain_dist);

    if(nebspec[0]==1)
    {/* fix constrained atoms for relaxing surrounding atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }
    
    LaVectorDouble DD(_CHAINLENGTH+1);
    LaVectorDouble DL(_CHAINLENGTH);
    LaVectorDouble DU(_CHAINLENGTH);
    LaVectorDouble b(_CHAINLENGTH+1), X(_CHAINLENGTH+1);
    LaGenMatDouble sigma(_CHAINLENGTH+1,3*n);
    LaTridiagFactDouble Tf;
    
#define get_Xi0(y0,y1,y2,y3) protect(\
        Xi0 = DL(0)*( y0/(_nebinterp[0]-_nebinterp[1])/(_nebinterp[0]-_nebinterp[2])/(_nebinterp[0]-_nebinterp[3]) \
                    + y1/(_nebinterp[1]-_nebinterp[0])/(_nebinterp[1]-_nebinterp[2])/(_nebinterp[1]-_nebinterp[3]) \
                    + y2/(_nebinterp[2]-_nebinterp[0])/(_nebinterp[2]-_nebinterp[1])/(_nebinterp[2]-_nebinterp[3]) \
                    + y3/(_nebinterp[3]-_nebinterp[0])/(_nebinterp[3]-_nebinterp[1])/(_nebinterp[3]-_nebinterp[2])); \
        )
#define get_XiN(y_nm3,y_nm2,y_nm1,y_n) protect(\
        XiN = DU(_CHAINLENGTH-1)*( y_nm3/(_nebinterp[_CHAINLENGTH-3]-_nebinterp[_CHAINLENGTH-2])/(_nebinterp[_CHAINLENGTH-3]-_nebinterp[_CHAINLENGTH-1])/(_nebinterp[_CHAINLENGTH-3]-_nebinterp[_CHAINLENGTH]) \
                    + y_nm2/(_nebinterp[_CHAINLENGTH-2]-_nebinterp[_CHAINLENGTH-3])/(_nebinterp[_CHAINLENGTH-2]-_nebinterp[_CHAINLENGTH-1])/(_nebinterp[_CHAINLENGTH-2]-_nebinterp[_CHAINLENGTH]) \
                    + y_nm1/(_nebinterp[_CHAINLENGTH-1]-_nebinterp[_CHAINLENGTH-3])/(_nebinterp[_CHAINLENGTH-1]-_nebinterp[_CHAINLENGTH-2])/(_nebinterp[_CHAINLENGTH-1]-_nebinterp[_CHAINLENGTH]) \
                    + y_n/(_nebinterp[_CHAINLENGTH]-_nebinterp[_CHAINLENGTH-3])/(_nebinterp[_CHAINLENGTH]-_nebinterp[_CHAINLENGTH-2])/(_nebinterp[_CHAINLENGTH]-_nebinterp[_CHAINLENGTH-1])); \
        )
#define get_cubicspline(h_i,w,wbar,sigma_i,sigma_im1,y_i,y_im1) protect(\
        y_csp = w*y_i + wbar*y_im1 + SQR(h_i)*((pow(w,3)-w)*sigma_i + (pow(wbar,3)-wbar)*sigma_im1); \
        )
               
    step0 = curstep;
    for(curstep=step0;curstep<(step0 + totalsteps);curstep++)    
    {
        /* get forces */
        for(j=0;j<=_CHAINLENGTH;j++)
        {
            s=_nebinterp[j];

            interpCN(s); /* interpolate all atoms for a given s */
            /* constrained atoms at chain value */
            copyRchaintoCN(j);

            if(nebspec[0]==1)
            {
                INFO("Relaxation for the point "<<j<<" in the energy hill");
                relax(); 
            }

            call_potential();
            Ec[j]=_EPOT;

            for(i=0;i<n;i++)
            {
                ipt=constrainatoms[i+1];
                _Fc[j][i]=_F[ipt];
            }
        }

        if(curstep%printfreq==0) 
        {
            INFO_Printf("curstep = %d\n",curstep);
            for(j=0;j<=_CHAINLENGTH;j++)
                INFO_Printf("%20.12e %25.15e %25.15e\n",_nebinterp[j],Ec[j]-Ec[0]);
            fprintf(fp,"%d ",curstep);
            for(j=0;j<=_CHAINLENGTH;j++)
                fprintf(fp,"%20.12e %25.15e ",_nebinterp[j], Ec[j]-Ec[0]);
            fprintf(fp,"\n");
            fflush(fp);            
        }
        
        /* move along force - This step is more like steepest descent relaxation. */
        for(j=0;j<=_CHAINLENGTH;j++)
           for(i=0;i<n;i++)
           {
                _Rc[j][i]+=_Fc[j][i]*_TIMESTEP;
           }

        copyRchaintoCN(0); setconfig1();
        copyRchaintoCN(_CHAINLENGTH); setconfig2();
        caldscom12();
        
        /* compute current constrain value */
        constrain_dist=0;
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
            constrain_dist+=ds0.norm2();
        }
        INFO_Printf("string method: constrain_dist=%e\n",constrain_dist);
        
        /* calculate current reaction coordinates s */
        Lc[0]=0;
        for(j=1;j<=_CHAINLENGTH;j++)
        {
            dRnorm2 = 0;
            for(i=0;i<n;i++)
            {
                dR = _Rc[j][i]-_Rc[j-1][i];
                dRnorm2 += dR.norm2();
            }
            dRnorm = sqrt(dRnorm2);
            Lc[j]=Lc[j-1]+dRnorm;
        }
        _nebinterp[0]=0;
        for(j=1;j<=_CHAINLENGTH;j++)
            _nebinterp[j] = Lc[j]/Lc[_CHAINLENGTH];

        /* compute DD, DL and DU which are tri-diagonals of Matrix T */
        DD(0)=-1; DU(0)=1;
        for(j=1;j<_CHAINLENGTH;j++)
        {
            DL(j-1)=_nebinterp[j]-_nebinterp[j-1];
            DU(j)=_nebinterp[j+1]-_nebinterp[j];
            DD(j)=2*(DL(j-1)+DU(j));
        }
        DL(_CHAINLENGTH-1)=-1; DD(_CHAINLENGTH)=1;
        LaTridiagMatDouble T(DD,DL,DU);
        LaTridiagMatFactorize (T, Tf);
        
        for(k=0;k<n;k++)
            for(i=0;i<3;i++)
            {/* compute the vector b, T*sigma=b */
                get_Xi0(_Rc[0][k][i],_Rc[1][k][i],_Rc[2][k][i],_Rc[3][k][i]);
                b(0) = Xi0;
                for(j=1;j<_CHAINLENGTH;j++)
                    b(j)=(_Rc[j+1][k][i]-_Rc[j][k][i])/DU(j) - (_Rc[j][k][i]-_Rc[j-1][k][i])/DL(j-1);
                get_XiN(_Rc[_CHAINLENGTH-3][k][i],_Rc[_CHAINLENGTH-2][k][i],_Rc[_CHAINLENGTH-1][k][i],_Rc[_CHAINLENGTH][k][i]);
                b(_CHAINLENGTH)=XiN;
            /* store sigma */
                LaLinearSolve(Tf,X,b);
                sigma(LaIndex(),3*k+i).inject(X);
            }

        for(k=0;k<n;k++)
            for(j=0;j<_CHAINLENGTH;j++)
                _RcStore[j][k] = _Rc[j][k];
        
        /* redistribut string evenly using spline interpolation
           - parametrization by equal arc length */
        for(j=1;j<_CHAINLENGTH;j++)
        {
            s_j = j*1.0/_CHAINLENGTH;
            for(i=j;i<_CHAINLENGTH;i++)
                if(_nebinterp[i]>s_j)
                {
                    hi = _nebinterp[i]-_nebinterp[i-1];
                    w = (s_j - _nebinterp[i-1])/hi; wbar = 1-w;
                    for(m=0;m<n;m++)
                        for(k=0;k<3;k++)
                        {
                            yi = _RcStore[i][m][k]; yim1 = _RcStore[i-1][m][k];
                            get_cubicspline(hi,w,wbar,sigma(3*m+k,i),sigma(3*m+k,i-1),yi,yim1);
                            _Rc[j][m][k]=y_csp;
                        }
                    
                    break;
                }
        }

        for(j=1;j<_CHAINLENGTH;j++)
            _nebinterp[j]=j*1.0/_CHAINLENGTH;
    }

    INFO_Printf("curstep = %d\n",curstep);
    for(j=0;j<=_CHAINLENGTH;j++)
        INFO_Printf("%20.12e %25.15e %25.15e\n",_nebinterp[j],Ec[j]-Ec[0]);
    fprintf(fp,"%d ",curstep);
    for(j=0;j<=_CHAINLENGTH;j++)
        fprintf(fp,"%20.12e %25.15e ",_nebinterp[j], Ec[j]-Ec[0]);
    fprintf(fp,"\n");
    fflush(fp); 

    free(Ec);free(Lc);free(_RcStore);
    fclose(fp);    
}


