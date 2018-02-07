/*
  bubble2d.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Apr 10 2014
*/

#include "bubble2d.h"

void BUBBLE2DFrame::initparser()
{
    MDPARALLELFrame::initparser();

    /* input */
    bindvar("radius0",&_radius0,DOUBLE);
    bindvar("sigma0",&_sigma0,DOUBLE);

    bindvar("kwall",&_kwall,DOUBLE);
    bindvar("slope",&_slope,DOUBLE);
    bindvar("channel",&_channel,DOUBLE);

}

int BUBBLE2DFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    return -1;
}

void BUBBLE2DFrame::initvars()
{
    MDPARALLELFrame::initvars();
}

void BUBBLE2DFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void BUBBLE2DFrame::bubble_2d()
{
    int i,j,ipt,jpt;
    double U,r,r2;
    double r0i, r0j, f0, border, x, y;
    Vector3 sij, rij, fij, fi;
    
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();
    
    _EPOT=0;

    for(i=0;i<_NP;i++)
    {_F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    f0 = _sigma0*_radius0;
    for(ipt=0;ipt<_NP;ipt++)
    {
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(ipt>jpt) continue;
            sij=_SR[jpt]-_SR[ipt];
            sij.z = 0;
            sij.subint();
            rij=_H*sij;
            r2=rij.norm2();
            r=sqrt(r2);

            r0i=_radius0; r0j=_radius0;
            if(r<=2*_radius0)
            {
                U = 0.5*f0/(r0i+r0j)*SQR( r0i+r0j - r );
                fij = rij*( f0/(r0i+r0j)*( r0i+r0j - r )/r );

                _F[ipt]-=fij;
                _F[jpt]+=fij;
                _EPOT_IND[ipt]+=U*0.5;
                _EPOT_IND[jpt]+=U*0.5;
                _EPOT+=U;
                _VIRIAL.addnvv(1.,fij,rij);
            }
        }
    }

    /* External wall potential */
    SHtoR();
    for(ipt=0;ipt<_NP;ipt++)
    {
	x=_R[ipt].x; y=_R[ipt].y;
        border = (SQR(y) - SQR(_slope*x)) - 10;
        if (border > 0)
        {
            U = 0.5*_kwall*SQR(border);
            fi.x = -_kwall*border*(-2*SQR(_slope)*x );
            fi.y = -_kwall*border*( 2*y );
            fi.z  = 0;

            _F[ipt]+=fi;
            _EPOT_IND[ipt]=U;
            _EPOT+=U;
        }
        if (y > _channel)
        {
            U = 0.5*_kwall*SQR(y-_channel);
            fi.x = 0;
            fi.y = -_kwall*( y-_channel );
            fi.z  = 0;

            _F[ipt]+=fi;
            _EPOT_IND[ipt]=U;
            _EPOT+=U;
        }
        if (y < -_channel)
        {
            U = 0.5*_kwall*SQR(y+_channel);
            fi.x = 0;
            fi.y = -_kwall*( y+_channel );
            fi.z  = 0;

            _F[ipt]+=fi;
            _EPOT_IND[ipt]=U;
            _EPOT+=U;
        }
    }


}


void BUBBLE2DFrame::potential()
{
    bubble_2d();
}

#ifdef _TEST

/* Main Program Begins */
class BUBBLE2DFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

