/*
  meam-baskes.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Dec  3 16:32:39 2007

  FUNCTION  : MD++ with MEAM potential (using Baskes Dyn88)
*/

#include "meam-baskes.h"


void MEAMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
}    

void MEAMFrame::potential()
{
    MEAM();
}

void MEAMFrame::readMEAM()
{
    LFile::SubHomeDir(mea_.meamf,mea_.meamf);
    LFile::SubHomeDir(mea_.meafile,mea_.meafile);
    myoutput_.fid = 24; /* redirect output to fort.24 */
    inter_(); /* read MEAM parameters */

    if (twonn_.nn) INFO("MEAM-Baskes: 2NN activated");
    
    if(cneigh_.dradn<0)
        cneigh_.dradn = 0.1*interact_.rcut;

    interact_.rcutsq = interact_.rcut*interact_.rcut;
    
    _RLIST= interact_.rcut + cneigh_.dradn;
    _SKIN = cneigh_.dradn;

    cneigh_.rctsqn = _RLIST*_RLIST;

/* change cmin in meafile */    
//    for(int i=0;i<types_.ntypes;i++)
//    {
//        if (cmin0[i]>0)
//        {
//            scr_.cmin[i][i][i]=cmin0[i];
//        }
//    }
//    
//    INFO_Printf("cmin[0][0][0]=%e\n",scr_.cmin[0][0][0]);
//    INFO_Printf("cmax[0][0][0]=%e\n",scr_.cmax[0][0][0]);

    WARNING("***************************************************************************");
    WARNING("The MEAM potential is only implemented for rectangular supercells!");
    WARNING("It won't work correctly for tilted supercells (e.g. with dislocation dipole)");
    WARNING("***************************************************************************");
}

void MEAMFrame::MEAM()
{
    Vector3 h;
    int i, j;
    double eps;

    /*refreshneighborlist();*/
    
    lattice_.natoms=_NP;
    
    if(_NP > natmax)
    {
      INFO_Printf("meam-baskes: NP = %d exceeds limit %d (meam-lammps does not have this limit)\n",_NP,natmax);
      return;
    }

    SHtoR();
    for(i=0;i<_NP;i++)
    {
        particle_.rv[i][0]=_R[i].x;
        particle_.rv[i][1]=_R[i].y;
        particle_.rv[i][2]=_R[i].z;
        particle_.itype[i] = species[i]+1;
    }
    lattice_.perlen[0]=_H[0][0];
    lattice_.perlen[1]=_H[1][1];
    lattice_.perlen[2]=_H[2][2];

    eps=1e-8;
    if ( (fabs(_H[0][1])>eps)||(fabs(_H[0][2])>eps)||
         (fabs(_H[1][0])>eps)||(fabs(_H[1][2])>eps)||
         (fabs(_H[2][0])>eps)||(fabs(_H[2][1])>eps) )
    {
        _INFO(_H);
        WARNING(" Only rectangular supercell is allowed!");
    }

//    INFO_Printf("lattice_.natoms=%d\n",lattice_.natoms);
    force_();

    _EPOT=0;
    for(i=0;i<_NP;i++)
    {
        _F[i].set(forces_.f[i][0],forces_.f[i][1],forces_.f[i][2]);
        _EPOT_IND[i] = forces_.e[i];
        _EPOT += forces_.e[i];
        _VIRIAL_IND[i].set(forces_.slocal[i][0][0],forces_.slocal[i][1][0],forces_.slocal[i][2][0],
                           forces_.slocal[i][0][1],forces_.slocal[i][1][1],forces_.slocal[i][2][1],
                           forces_.slocal[i][0][2],forces_.slocal[i][1][2],forces_.slocal[i][2][2]);
    }
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            _VIRIAL[i][j]=forces_.stresst[j][i];
//            _VIRIAL[i][j]=forces_.stresst[i][j];
}

void MEAMFrame::initvars()
{
    MDPARALLELFrame::initvars();
}

void MEAMFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("meamfile",mea_.meamf,STRING);
    bindvar("meafile",mea_.meafile,STRING);
    bindvar("ntypes",&types_.ntypes,INT);
    bindvar("ename0",cmeam_.enames[0],STRING);
    bindvar("ename1",cmeam_.enames[1],STRING);
    bindvar("ename2",cmeam_.enames[2],STRING);
    bindvar("rcut",&interact_.rcut,DOUBLE);
    //bindvar("cmin0",cmin0,DOUBLE);    
    bindvar("kode0",kodes_.kodes[0],STRING);    
    bindvar("kode1",kodes_.kodes[0],STRING);    
    bindvar("kode2",kodes_.kodes[0],STRING);    
}

int MEAMFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readMEAM",readMEAM());
    bindcommand(name,"printpairpot",printpairpot());
    
    return -1;
}
   
void MEAMFrame::printpairpot()
{    
    int elti, eltj;
    double rmin, dr, rmax, r, phi, phip;
    char pairfile[200];
    FILE *fp;
    
    elti = (int) input[0] + 1;    // atom species i
    eltj = (int) input[1] + 1;    // atom species j
    rmin = input[2];
    dr   = input[3];
    rmax = input[4];

    strcpy(pairfile,"pairpot_baskes.dat");
    fp=fopen(pairfile,"w");
    if(fp==NULL)
    {
        FATAL("printpairpot: file "<<pairfile<<" open failure.");
    }
    if((rmax<rmin)||(dr<=0))
        FATAL("rmax cannot be smaller than rmin, dr must be positive");
    
    phip = 0;
    for(r=rmin;r<=rmax;r+=dr)
    {
        //phi = phif_(&r,&elti,&eltj);
        //phi = phiid_(&r,&elti);
        calphiid_(&r,&elti,&phi);
        fprintf(fp,"%21.14e  %21.14e %21.14e\n",r,phi,phip);
    }
    fclose(fp);
    INFO("results written to "<<pairfile);
}


#ifdef _TEST

/* Main Program Begins */
class MEAMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST






























































/* old code */
#if 0

inline double interp(double func[],double deriv[],double dr,int ind,double qq)
{
    double f, a, b, c, d, f1, p1, A1, A2, dr2, dr3, qq2, qq3;
//    f = func[ind] + qq*deriv[ind];
    dr2=dr*dr; dr3=dr2*dr;
    qq2=qq*qq; qq3=qq2*qq;
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    f=a+b*qq+c*qq2+d*qq3;
    return f;
}
    
    
inline double interp1(double func[],double deriv[],double dr,int ind,double qq)
{
    double fp, a, b, c, d, f1, p1, A1, A2, dr2, dr3, qq2, qq3;
//    fp = deriv[ind] + qq1/dr*(deriv[ind+1]-deriv[ind]);
    dr2=dr*dr; dr3=dr2*dr;
    qq2=qq*qq; qq3=qq2*qq;
    a = func[ind];
    b = deriv[ind];
    f1 = func[ind+1];
    p1 = deriv[ind+1];
    A1 = f1-a-b*dr;
    A2 = (p1-b)*dr;
    d = (A2-2*A1)/dr3;
    c = (3*A1-A2)/dr2;
    fp=b+2*c*qq+3*d*qq2;
    return fp;
}
#endif
