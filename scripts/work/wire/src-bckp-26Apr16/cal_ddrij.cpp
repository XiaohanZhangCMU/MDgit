/*
 * cal_ddrij.cpp
 * by Keonwook Kang kwkang@lanl.gov 
 * Last Modified : Wed May 18 22:36:20 MDT 2011
 *
 * Purpose:
 *    Calculate displacement difference vectors
 */

#include <stdio.h>
#include "general.h"
#include "linalg3.h"
#include "md.h"
#include "cal_ddrij.h"

//using namespace MDPP_NS;

ComputeDisregistry::ComputeDisregistry()
{
    nb = 0;
    _RLIST = 5.43;
}

ComputeDisregistry::~ComputeDisregistry()
{
  delete [] bonds;
  delete [] ddrij_sav;
}

void ComputeDisregistry::alloc()
{
    if(bonds==NULL)
        bonds = new BONDS [MAXBONDS];
    if(ddrij_sav==NULL)
        ddrij_sav = new DDRIJ [MAXBONDS];
}

void ComputeDisregistry::initparser()
{
    MDFrame::initparser();
    bindvar("slip_normal",slip_normal,DOUBLE);
}

int ComputeDisregistry::exec(const char *name)
{
    if(MDFrame::exec(name)==0) return 0;
    bindcommand(name,"identify_bonds",{alloc(); identify_bonds();});
    bindcommand(name,"write_bonds",write_bonds(finalcnfile));
    bindcommand(name,"read_bonds",read_bonds(incnfile));
    bindcommand(name,"plot_bonds",plot_bonds());
    bindcommand(name,"calddrij",calddrij());
    bindcommand(name,"retag_ddrij",retag_ddrij());
    bindcommand(name,"write_ddrij",write_ddrij(finalcnfile));
    bindcommand(name,"read_ddrij",read_ddrij(incnfile));
    bindcommand(name,"plot_ddrij",plot_ddrij());
    return -1;
}

void ComputeDisregistry::initvars()
{
    MDFrame::initvars();
}

void ComputeDisregistry::calddrij()
{
    int i, ib, ipt, jpt;
    double rc;
    Matrix33 h0, h1; 
    Vector3 sri, srj, sri0, srj0, dsrij, dsrij0;
    Vector3 ri, rj, ri0, rj0, rij, rij0, ddrij;

    rc = input[0];

    /* store the current configuration box */
    h1 = _H;
    /* store the current config */
    for(i=0;i<_NP;i++)
        _VR[i] = _SR[i];

    /* read reference structure */
    initcn.open(incnfile,LFile::O_Read);
    initcn.read(this);
    h0 = _H;

    /* calculate ddrij = drj - dri
                       = (rj - rj0) - (ri - ri0)
                       = (rj - ri) - (rj0 - ri0) 
                       = rij - rij0 */
    Realloc(ddrij_sav,DDRIJ,nb);
    for(ib=0;ib<nb;ib++)
    {
        ipt=bonds[ib].i; jpt=bonds[ib].j;

        sri0 = _SR[ipt]; srj0 = _SR[jpt];
        dsrij0 = srj0 - sri0; dsrij0.subint();
        rij0 = h0*dsrij0;

        sri = _VR[ipt]; srj = _VR[jpt];
        dsrij = srj - sri; dsrij.subint();
        rij = h1*dsrij;

        ddrij = rij - rij0;

        // debug
        //if(ib==0)
        //{
	//    INFO("rij   = [ "<<rij<<" ]");
	//    INFO("rij0  = [ "<<rij0<<" ]");
	//    INFO("ddrij = [ "<<ddrij<<" ]");
        //}

        ddrij_sav[ib].disreg = ddrij;
        ddrij_sav[ib].norm = ddrij.norm();
        if(ddrij_sav[ib].norm>rc)
	    ddrij_sav[ib].tag = 1;
	else
	    ddrij_sav[ib].tag = 0;
    }
}

void ComputeDisregistry::retag_ddrij()
{
    int ib;
    double rc;

    rc = input[0];
    for(ib=0;ib<nb;ib++)
    {
        if(ddrij_sav[ib].norm>rc)
	    ddrij_sav[ib].tag = 1;
	else
	    ddrij_sav[ib].tag = 0;
    }
}

void ComputeDisregistry::write_ddrij(char *fname)
{
    FILE *fp;
    int ib;
    
    INFO("MDFrame::write_ddrij "<<fname);

    fp=fopen(fname,"w");
    if(fp==NULL)
        FATAL("write_ddrij: open file failure");

    fprintf(fp,"%d\n",nb);
    for(ib=0;ib<nb;ib++)
        fprintf(fp,"%20.16e %20.16e %20.16e %20.16e %d\n",
                ddrij_sav[ib].disreg.x,ddrij_sav[ib].disreg.y,ddrij_sav[ib].disreg.z, 
                ddrij_sav[ib].norm, ddrij_sav[ib].tag);
    fclose(fp);
    DUMP("write_ddrij finished");
}

void ComputeDisregistry::read_ddrij(char *fname)
{
    FILE *fp;
    int ib;

    INFO("MDFrame::read_ddrij "<<fname);

    fp=fopen(fname,"r");
    if(fp==NULL)
        FATAL("read_ddrij: open file failure");

    fscanf(fp,"%d", &nb);
    INFO("Number of bonds: nb = "<<nb);
    Realloc(ddrij_sav,DDRIJ,nb);

    for(ib=0;ib<nb;ib++)
    {
      fscanf(fp,"%le %le %le %le %d",
             &ddrij_sav[ib].disreg.x,&ddrij_sav[ib].disreg.y,&ddrij_sav[ib].disreg.z, 
             &ddrij_sav[ib].norm,&ddrij_sav[ib].tag);
           
    }
    fclose(fp);
    DUMP("read_ddrij finished");
}

void ComputeDisregistry::plot_ddrij()
{
    int i;
    double L, x1, y1, z1, x2, y2, z2, dx, dy, dz, dr;
    char s[100];
    unsigned ce;
    Vector3 s1, s2, vr1, vr2, sri, srj, ri, rj, ddrij, n;

    if(win==NULL) return;
    if(!(win->alive)) return;
    win->Lock();
    win->Clear();

    n.set(slip_normal[0], slip_normal[1], slip_normal[2]);
    n/=n.norm();

    L=max(_H[0][0], _H[1][1]);
    L=max(L,_H[2][2])*.5;

    /* draw atoms */
    for(i=0;i<_NP;i++)
    {
        sri=_SR[i];
        if(plot_map_pbc==1) sri.subint();
        ri = _H*sri;

        if(plot_limits[0]==1)
            if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6]))
            {
                continue;
            }

	/* Color atoms according to their species */
        ce=colors[MAXCOLORS+5+species[i]];
        sprintf(s,"%d, %d, %6.4f %6.4f %6.4f %x",
                i,species[i],ri.x,ri.y,ri.z,ce);
        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                      atomradius[species[i]]/L,ce,s,2);
    }

    /* draw ddrij */
    for(i=0;i<nb;i++)
    {
        if(ddrij_sav[i].tag == 0) continue;
        
        ddrij=ddrij_sav[i].disreg;
        //sri = _SR[bonds[i].i]; 
        srj = _SR[bonds[i].j]; 
   
        if(plot_map_pbc==1) 
	{
	  //sri.subint(); 
          srj.subint();
        }

        if(plot_limits[0]==1)
            if((srj.x<plot_limits[1])||(srj.x>plot_limits[2])
               ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
               ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
            {
                continue;
            }

        //ri = _H*sri;
        rj = _H*srj;

        x1=rj.x/L; y1=rj.y/L; z1=rj.z/L;
        x2=(rj.x+ddrij.x)/L; y2=(rj.y+ddrij.y)/L; z2=(rj.z+ddrij.z)/L;
        dx=x2-x1; dy=y2-y1; dz=z2-z1; dr=sqrt(dx*dx+dy*dy+dz*dz);
        dx/=dr; dy/=dr; dz/=dr;
        win->DrawArrow(x1+dx*atomradius[species[i]]/L,
                       y1+dy*atomradius[species[i]]/L,
                       z1+dz*atomradius[species[i]]/L,
                       x2, y2, z2,
                       colors[MAXCOLORS+bonds[i].type+1],0.001,.2,1);
    }

    /* draw frame */
#define drawarrow(a,b,c,d,e,f,g,h,i,j) s1.set(a,b,c); s2.set(d,e,f);	\
        vr1=_H*s1; vr2=_H*s2; vr1/=2*L; vr2/=2*L;\
        win->DrawArrow(vr1.x,vr1.y,vr1.z,vr2.x,vr2.y,vr2.z,g,h,i,j);
#define drawsline(a,b,c,d,e,f,g,h,i) s1.set(a,b,c); s2.set(d,e,f);\
        vr1=_H*s1; vr2=_H*s2; vr1/=2*L; vr2/=2*L;\
        win->DrawLine(vr1.x,vr1.y,vr1.z,vr2.x,vr2.y,vr2.z,g,h,i);

    ce=win->AllocShortRGBColor(255,0,0);  // red
    drawarrow(-1,-1,-1, 1,-1,-1,ce,0,10.0/_H[0][0],0);
    ce=win->AllocShortRGBColor(0,255,0);  // green
    drawarrow(-1,-1,-1,-1, 1,-1,ce,0,10.0/_H[1][1],0);
    ce=win->AllocShortRGBColor(0,0,255);  // blue
    drawarrow(-1,-1,-1,-1,-1, 1,ce,0,10.0/_H[2][2],0);
    drawsline(-1,-1, 1,-1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1,-1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1,-1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1, 1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    if((fabs(n.x)>1e-12) || (fabs(n.y)>1e-12) || (fabs(n.z)>1e-12))
    {
        ce=win->AllocShortRGBColor(255,255,0);  // yellow
        drawarrow(-1,-1,-1,-1+n.x,-1+n.y,-1+n.z,ce,0,10.0/L,0);
    }

    win->Unlock();
    win->Refresh();
}

int ComputeDisregistry::identify_bonds()
{
    Vector3 sri, srj, srij, rij, sbond, n;
    double rc_AA, rc_BB, rc_AB, r2, rc2;
    int ipt, jpt, j, bond_type;

    rc_AA = input[0];
    rc_BB = input[1];
    rc_AB = input[2];
    n.set(slip_normal[0], slip_normal[1], slip_normal[2]);

    for(ipt=0;ipt<_NP;ipt++)
    {
        if(plot_limits[0]==1)
	{
	    if((_SR[ipt].x<plot_limits[1])||(_SR[ipt].x>plot_limits[2])
               ||(_SR[ipt].y<plot_limits[3])||(_SR[ipt].y>plot_limits[4])
               ||(_SR[ipt].z<plot_limits[5])||(_SR[ipt].z>plot_limits[6]))
	        continue;
        }
        sri = _SR[ipt];     
        for(j=0;j<nn[ipt];j++)
	{
  	    jpt = nindex[ipt][j];
            if (ipt>=jpt) continue;

            srj = _SR[jpt];
            srij = sri - srj;
            srij.subint();
            rij = _H*srij;
            r2=rij.norm2();
         
            if(species[ipt]==0 && species[jpt]==0)
	    {
                rc2 = rc_AA*rc_AA; bond_type = 0;
	    }
            else if(species[ipt]==1 && species[jpt]==1)
	    {
                rc2 = rc_BB*rc_BB; bond_type = 1;
	    }
            else
	    {
                rc2 = rc_AB*rc_AB; bond_type = 2;
	    }

            if(r2<rc2)
	    {
	        if (dot(n,rij)>0)
		{// flip the direction of the bond vector
                    bonds[nb].i=  jpt;
                    bonds[nb].j = ipt;
                }
                else
		{
	            bonds[nb].i=  ipt;
                    bonds[nb].j = jpt;
                }
                sbond = srj; sbond += (srij*.5); 
                sbond.subint();
                bonds[nb].center = _H*sbond;
                bonds[nb].type = bond_type;
 
                // debug
                //if ((nb%1000)==0)
                //INFO_Printf("Bond ID [%d]: %d %d %12.8f %12.8f %12.8f %d\n",
                //            nb, bonds[nb].i , bonds[nb].j, 
                //            bonds[nb].center.x, bonds[nb].center.y, bonds[nb].center.z,
                //            bonds[nb].type);
                nb++;
	    }

            if (nb >= MAXBONDS)
	    {
	        INFO("Too many lines (more than "<<MAXBONDS<<"). Increase MAXBONDS");
                return -1;
            }
	}
    }
    INFO(nb<<" bonds identified");
    return 0;
}

void ComputeDisregistry::write_bonds(char *fname)
{
    FILE *fp;
    int ib;
    
    INFO("MDFrame::write_bonds "<<fname);

    fp=fopen(fname,"w");
    if(fp==NULL)
        FATAL("write_bonds: open file failure");

    fprintf(fp,"%d\n",nb);
    for(ib=0;ib<nb;ib++)
        fprintf(fp,"%d %d %20.16e %20.16e %20.16e %d\n",
                bonds[ib].i,bonds[ib].j,
                bonds[ib].center.x,bonds[ib].center.y,bonds[ib].center.z,
                bonds[ib].type);
    fclose(fp);
    DUMP("write_bonds finished");
}

void ComputeDisregistry::read_bonds(char *fname)
{
    FILE *fp;
    int ib;

    INFO("MDFrame::read_bonds "<<fname);

    fp=fopen(fname,"r");
    if(fp==NULL)
        FATAL("read_bonds: open file failure");

    fscanf(fp,"%d", &nb);
    INFO("Number of bonds: nb = "<<nb);
    Realloc(bonds,BONDS,nb);

    for(ib=0;ib<nb;ib++)
    {
      fscanf(fp,"%d %d %le %le %le %d",&bonds[ib].i,&bonds[ib].j,
             &bonds[ib].center.x,&bonds[ib].center.y,&bonds[ib].center.z,
             &bonds[ib].type);             
    }
    fclose(fp);

    // debug
    //for(ib=0;ib<nb;ib++)
    //{
    //    if ((ib%1000)==0)
    //        INFO_Printf("Bond ID [%d]: %d %d %12.8f %12.8f %12.8f %d\n",
    //                    nb, bonds[ib].i , bonds[ib].j, 
    //                    bonds[ib].center.x, bonds[ib].center.y, bonds[ib].center.z,
    //                    bonds[ib].type);
    //}
    DUMP("read_bonds finished");
}

void ComputeDisregistry::plot_bonds()
{
    int i;
    double L, x1, y1, z1, x2, y2, z2, dx, dy, dz, dr, r2;
    char s[100];
    unsigned ce;
    Vector3 s1, s2, vr1, vr2, sri, srj, ri, rj;

    if(win==NULL) return;
    if(!(win->alive)) return;
    win->Lock();
    win->Clear();

    L=max(_H[0][0], _H[1][1]);
    L=max(L,_H[2][2])*.5;

    /* draw atoms */
    for(i=0;i<_NP;i++)
    {
        sri=_SR[i];
        if(plot_map_pbc==1) sri.subint();
        ri = _H*sri;

        if(plot_limits[0]==1)
            if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6]))
            {
                continue;
            }
	/* Color atoms according to their species */
        ce=colors[MAXCOLORS+5+species[i]];
        sprintf(s,"%d, %d, %6.4f %6.4f %6.4f %x",
                i,species[i],ri.x,ri.y,ri.z,ce);
        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                      atomradius[species[i]]/L,ce,s,2);
    }

    /* draw bonds */
    for(i=0;i<nb;i++)
    {
        sri = _SR[bonds[i].i]; srj = _SR[bonds[i].j];      
        if(plot_map_pbc==1) 
	{
	  sri.subint(); srj.subint();
        }
        ri = _H*sri; rj = _H*srj;

        r2=(ri-rj).norm2();
        if(r2>L) continue;
        x1=ri.x/L; y1=ri.y/L; z1=ri.z/L;
        x2=rj.x/L; y2=rj.y/L; z2=rj.z/L;
        dx=x2-x1; dy=y2-y1; dz=z2-z1; dr=sqrt(dx*dx+dy*dy+dz*dz);
        dx/=dr; dy/=dr; dz/=dr;
        win->DrawLine(x1+dx*atomradius[species[i]]/L,
                       y1+dy*atomradius[species[i]]/L,
                       z1+dz*atomradius[species[i]]/L,
                       x2-dx*atomradius[species[i]]/L,
                       y2-dy*atomradius[species[i]]/L,
                       z2-dz*atomradius[species[i]]/L,
                       colors[MAXCOLORS+bonds[i].type+1],bondradius/L,1);
    }

    /* draw frame */
#define drawsline(a,b,c,d,e,f,g,h,i) s1.set(a,b,c); s2.set(d,e,f);\
        vr1=_H*s1; vr2=_H*s2; vr1/=2*L; vr2/=2*L;\
        win->DrawLine(vr1.x,vr1.y,vr1.z,vr2.x,vr2.y,vr2.z,g,h,i);

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

    win->Unlock();
    win->Refresh();
}

void ComputeDisregistry::potential()
{
    INFO("Empty potential");
}

/* Main Program Begins */
#ifdef _TEST
class ComputeDisregistry sim;

/* The main program is defined here */
#include "main.cpp"

#endif //_TEST
