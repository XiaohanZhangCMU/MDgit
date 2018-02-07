/*
  rods.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed May  2 19:03:18 2007

  FUNCTION  :  MD simulation package of Lennard-Jones (LJ) potential
               LJ parameters depend on species
               Atoms can also form covalent bond with neighbors

  Test cases: scripts/work/ar/lj2-md.script, ljbond.tcl
*/

#include "rods.h"

void RODSFrame::initparser()
{
    MDPARALLELFrame::initparser();
    bindvar("bond_k",&BOND_K,DOUBLE);
    bindvar("repul",&REPUL,DOUBLE);
}

int RODSFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"init_structure",init_structure());
    bindcommand(name,"final_structure",final_structure());
    bindcommand(name,"connect_nodes",connect_nodes());
    bindcommand(name,"zero_com_rotation",zero_com_rotation());
    return -1;
}

void RODSFrame::Alloc()
{
    int size;
    MDPARALLELFrame::Alloc();

    size = _NP*allocmultiple;    
    Realloc(_Fext,Vector3,size);
    bindvar("Fext",_Fext,DOUBLE);
}

void RODSFrame::Alloc_Bonds()
{
    int size;
    size = NUM_BONDS;
    Realloc(BOND_R0, double,size);
    Realloc(BOND_INDEX, int,size*2);

    memset(BOND_R0, 0,sizeof(double)*size);
    memset(BOND_INDEX,  0,sizeof(int)*size*2);
}


void RODSFrame::initvars()
{
    _RLIST=1.1;
    _SKIN =1.0;
    MDPARALLELFrame::initvars();
}

void RODSFrame::init_structure()
{
    _NP = 24;
    
    Alloc();

    /* Assign positions of the nodes */
    _R[ 0].set(-1,  -1, 0);
    _R[ 1].set( 0,  -1, 0);
    _R[ 2].set( 1,  -1, 0);
    _R[ 3].set(-1,-0.6, 0);
    _R[ 4].set( 0,-0.4, 0);
    _R[ 5].set( 1,-0.6, 0);
    _R[ 6].set(-1,   0, 0);
    _R[ 7].set( 0,   0, 0);
    _R[ 8].set( 1,   0, 0);
    _R[ 9].set(-1, 0.6, 0);
    _R[10].set( 0, 0.4, 0);
    _R[11].set( 1, 0.6, 0);
    _R[12].set(-1,   1, 0);
    _R[13].set( 0,   1, 0);
    _R[14].set( 1,   1, 0);

    _R[15].set(-1,-0.6, 0);
    _R[16].set( 0,-0.4, 0);
    _R[17].set( 1,-0.6, 0);
    _R[18].set(-1,   0, 0);
    _R[19].set( 0,   0, 0);
    _R[20].set( 1,   0, 0);
    _R[21].set(-1, 0.6, 0);
    _R[22].set( 0, 0.4, 0);
    _R[23].set( 1, 0.6, 0);
    
    _H.clear(); _H[0][0]=4; _H[1][1]=4; _H[2][2]=4;

    RHtoS();

    connect_nodes();
}

void RODSFrame::final_structure()
{
    _NP = 24;
    
     /* Assign positions of the nodes */
    _R[ 0].set(-1,  -0.6, 0);
    _R[ 1].set( 0,  -0.4, 0);
    _R[ 2].set( 1,  -0.6, 0);
    _R[ 3].set(-1,-0.6, 0.4);
    _R[ 4].set( 0,-0.4, 0.6);
    _R[ 5].set( 1,-0.6, 0.4);
    _R[ 6].set(-1,   0, 0.4);
    _R[ 7].set( 0,   0, 0.6);
    _R[ 8].set( 1,   0, 0.4);
    _R[ 9].set(-1, 0.6, 0.4);
    _R[10].set( 0, 0.4, 0.6);
    _R[11].set( 1, 0.6, 0.4);
    _R[12].set(-1,   0.6, 0);
    _R[13].set( 0,   0.4, 0);
    _R[14].set( 1,   0.6, 0);

    _R[15].set(-1,-0.6, -0.4);
    _R[16].set( 0,-0.4, -0.6);
    _R[17].set( 1,-0.6, -0.4);
    _R[18].set(-1,   0, -0.4);
    _R[19].set( 0,   0, -0.6);
    _R[20].set( 1,   0, -0.4);
    _R[21].set(-1, 0.6, -0.4);
    _R[22].set( 0, 0.4, -0.6);
    _R[23].set( 1, 0.6, -0.4);
    
    _H.clear(); _H[0][0]=4; _H[1][1]=4; _H[2][2]=4;

    RHtoS();
}

void RODSFrame::connect_nodes()
{
    int i;
    Vector3 dr;
    int bond_index_data[112] = { 0,1,
                                 1,2,
                                 0,3,
                                 0,4,
                                 1,4,
                                 2,4,
                                 2,5,
                                 3,4,
                                 4,5,
                                 3,6,
                                 3,7,
                                 4,7,
                                 5,7,
                                 5,8,
                                 6,7,
                                 7,8,
                                 6,9,
                                 7,9,
                                 7,10,
                                 7,11,
                                 8,11,
                                 9,10,
                                 10,11,
                                 9,12,
                                 10,12,
                                 10,13,
                                 10,14,
                                 11,14,
                                 12,13,
                                 13,14,
                                 0,15,
                                 0,16,
                                 1,16,
                                 2,16,
                                 2,17,
                                 15,16,
                                 16,17,
                                 15,18,
                                 15,19,
                                 16,19,
                                 17,19,
                                 17,20,
                                 18,19,
                                 19,20,
                                 18,21,
                                 19,21,
                                 19,22,
                                 19,23,
                                 20,23,
                                 21,22,
                                 22,23,
                                 21,12,
                                 22,12,
                                 22,13,
                                 22,14,
                                 23,14 };
    
    NUM_BONDS = 56;
    Alloc_Bonds();
    
    /* Assign bond indices */
    for(i=0;i<NUM_BONDS;i++)
    {
        BOND_INDEX[i*2]   = bond_index_data[i*2];
        BOND_INDEX[i*2+1] = bond_index_data[i*2+1];
    }
    
    /* Compute distances of the bonds */
    SHtoR();
    for(i=0;i<NUM_BONDS;i++)
    {
        dr = _R[BOND_INDEX[i*2]]-_R[BOND_INDEX[i*2+1]];
        BOND_R0[i] = dr.norm();
        INFO_Printf("BOND_R0[%d]=%20.10e\n",i,BOND_R0[i]);
    }

}

void RODSFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void RODSFrame::rods_potential()
{
    int i,ipt,jpt;
    double U,r,r2;
    Vector3 rij, fij;
    
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();
    
    _EPOT=0;

    for(i=0;i<_NP;i++)
    {_F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    SHtoR();
    for(i=0;i<NUM_BONDS;i++)
    {
        ipt = BOND_INDEX[i*2];
        jpt = BOND_INDEX[i*2+1];
        rij = _R[jpt] - _R[ipt];

        r2=rij.norm2();
        r=sqrt(r2);
                
        U = 0.5* BOND_K * SQR(r-BOND_R0[i]);
        fij = rij * (-1.0*BOND_K*(1-BOND_R0[i]/r));

        _F[ipt]-=fij;
        _F[jpt]+=fij;
        _EPOT_IND[ipt]+=U*0.5;
        _EPOT_IND[jpt]+=U*0.5;
        _EPOT+=U;
        _VIRIAL.addnvv(1.,fij,rij);
    }

    /* Adding repulsive force between all atom pairs */
    for(ipt=0;ipt<_NP;ipt++)
        for(jpt=ipt+1;jpt<_NP;jpt++)
        {
            rij = _R[jpt] - _R[ipt];

            r2=rij.norm2();
            
            U = -0.5*REPUL* r2;
            fij = rij*REPUL;

            _F[ipt]-=fij;
            _F[jpt]+=fij;
            _EPOT_IND[ipt]+=U*0.5;
            _EPOT_IND[jpt]+=U*0.5;
            _EPOT+=U;
            _VIRIAL.addnvv(1.,fij,rij);
        }
    
    /* external force */
    for(i=0;i<_NP;i++)
    {
        if(!fixed[i])
        {
            _F[i] += _Fext[i];
            _EPOT -= dot(_Fext[i],_R[i]);
        }
    }

    /* fixed atoms */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) _F[i].clear();
    }
    
}


void RODSFrame::potential()
{
    rods_potential();
}

void RODSFrame::zero_com_rotation()
{
    int i, n;
    double Jx, Jy, Jz;
    Matrix33 hinv;
    
    for(i=0;i<_NP;i++)
    {
        _R[i]=_H*_SR[i];
        _VR[i]=_H*_VSR[i];
    }

    /* computing COM angular momentum */
    Jz = 0;  n = 0;
    for(i=0;i<_NP;i++)
        if(!fixed[i])
        {
            Jz += _R[i].x*_VR[i].y - _R[i].y*_VR[i].x;
            n++;
        }
    Jz /= n;

    Jx = 0;  n = 0;
    for(i=0;i<_NP;i++)
        if(!fixed[i])
        {
            Jx += _R[i].y*_VR[i].z - _R[i].z*_VR[i].y;
            n++;
        }
    Jx /= n;

    Jy = 0;  n = 0;
    for(i=0;i<_NP;i++)
        if(!fixed[i])
        {
            Jy += _R[i].z*_VR[i].x - _R[i].x*_VR[i].z;
            n++;
        }
    Jy /= n;
    
    for(i=0;i<_NP;i++)
        if(!fixed[i])
        {
            _VR[i].y -= Jz*_R[i].x/(_R[i].x*_R[i].x+_R[i].y*_R[i].y);
            _VR[i].x += Jz*_R[i].y/(_R[i].x*_R[i].x+_R[i].y*_R[i].y);
        }
    
#if 0 /* not working yet, need to debug */
    for(i=0;i<_NP;i++)
        if(!fixed[i])
        {
            _VR[i].z -= Jx*_R[i].y/(_R[i].y*_R[i].y+_R[i].z*_R[i].z);
            _VR[i].y += Jx*_R[i].z/(_R[i].y*_R[i].y+_R[i].z*_R[i].z);
        }

    for(i=0;i<_NP;i++)
        if(!fixed[i])
        {
            _VR[i].x -= Jy*_R[i].z/(_R[i].x*_R[i].x+_R[i].z*_R[i].z);
            _VR[i].z += Jy*_R[i].x/(_R[i].x*_R[i].x+_R[i].z*_R[i].z);
        }
#endif
    
    
    hinv = _H.inv();
    for(i=0;i<_NP;i++)
    {
        if(!fixed[i])
        {
            _SR[i]=hinv*_R[i];
            _VSR[i]=hinv*_VR[i];
        }
    }
}

void RODSFrame::plot()
{
    int i,ipt,jpt;
    double L;
    char s[100];
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    unsigned ce;
    Vector3 s1, s2, vr1, vr2, sri, srj, ri, rj;

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=max(_H[0][0],_H[1][1]);
    L=max(L,_H[2][2])*.5;

    switch(plot_color_axis){
     case(0): color_ind=_EPOT_IND; break;
     case(1): color_ind=_EPOT_IND; break;
     case(2): color_ind=_TOPOL; break;
     case(3): color_ind=_TOPOL; break;
     default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    SHtoR();
    win->Lock();
    win->Clear();
    
    /* draw atom s*/
    for(i=0;i<_NP;i++)
    {
        sri=_SR[i];
        if(plot_map_pbc==1) sri.subint();
        ri = _H*sri;
                
        if (!plot_atom_info) s[0]=0;
        else
        {   /* when an atom is clicked */
            if(plot_atom_info==1)    /* print scaled coordinates of the atom */
                sprintf(s,"%d,%6.4f,%6.4f,%6.4f",
                        i,sri.x,sri.y,sri.z);
            else if(plot_atom_info==2) /* print real coordinates of the atom */
                sprintf(s,"%d,%6.4f,%6.4f,%6.4f",
                        i,ri.x,ri.y,ri.z);
            else if(plot_atom_info==3) /* print color_ind, fixed, species info of the atom */
                sprintf(s,"%d,%8.6e %d %d",
                        i,color_ind[i], fixed[i], species[i]);
            else if(plot_atom_info==4) /* print force on the atom */
                sprintf(s,"%d,%6.4f,%6.4f,%6.4f",
                        i,_F[i].x,_F[i].y,_F[i].z);
            else
                s[0]=0;
        }
        /* otherwise color atoms according to their species */        
        ce=colors[MAXCOLORS+5+species[i]];
        if(plot_atom_info==5) /* print atom species and color code if clicked */
            sprintf(s,"%d,%d,%x",i,species[i],ce);
        
        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                       atomradius[species[i]]/L,ce,s,2);
    }
        
    /* Draw Bonds */
    for(i=0;i<NUM_BONDS;i++)
    {
        ipt = BOND_INDEX[i*2];
        jpt = BOND_INDEX[i*2+1];
        ri = _R[ipt];
        rj = _R[jpt];

        x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
        x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
        dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
        dx/=dr;dy/=dr;dz/=dr;
        win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                      y1+dy*atomradius[species[ipt]]/L,
                      z1+dz*atomradius[species[ipt]]/L,
                      x2-dx*atomradius[species[jpt]]/L,
                      y2-dy*atomradius[species[jpt]]/L,
                      z2-dz*atomradius[species[jpt]]/L,
                      colors[MAXCOLORS+1],bondradius/L,1);
    }

    /* draw frame */
#define drawsline(a,b,c,d,e,f,g,h,i) s1.set(a,b,c); s2.set(d,e,f);\
        vr1=_H*s1; vr2=_H*s2; vr1/=2*L; vr2/=2*L;\
        win->DrawLine(vr1.x,vr1.y,vr1.z,vr2.x,vr2.y,vr2.z,g,h,i);
#if 0    
    drawsline(-1,-1,-1,-1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1,-1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1,-1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1,-1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1,-1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1, 1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
#endif
    
    win->Unlock();
    win->Refresh();
}

    

    
#ifdef _TEST

/* Main Program Begins */
class RODSFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

