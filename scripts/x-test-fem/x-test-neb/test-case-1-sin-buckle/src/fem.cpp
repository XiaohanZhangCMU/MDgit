/*
  fem.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Mar 11 21:07:00 2016

  FUNCTION  :  FEM model of hyperelastic material with brick elements
  NOTE :       All elements must be the same shape and size for now.
*/

#include "fem.h"

void FEMFrame::initparser()
{
    MDPARALLELFrame::initparser();

    /* input */
    int _NELE;              /* number of elements */
    int _NNODE_PER_ELEMENT; /* number of nodes per element */
    int _NINT_PER_ELEMENT;  /* number of Gauss integration points per element */
    bindvar("nele",&_NELE,INT);
    bindvar("nnode_per_element",&_NNODE_PER_ELEMENT,INT);
    bindvar("nint_per_element",&_NINT_PER_ELEMENT,INT);
    bindvar("elements_file",elements_file,STRING);
    bindvar("fem_coeff_file",fem_coeff_file,STRING);
}

void FEMFrame::initvars()
{
       _RLIST=1.1; _SKIN =0.1; /* for visualization purposes only */

       sprintf(ELEMENT_TYPE,"CPE4");
       sprintf(ELEMENT_INT_TYPE,"Q4");

       DUMP(HIG"FEMFrame initvars"NOR);
       MDPARALLELFrame::initvars();
}

int FEMFrame::exec(const char *name)
{
  if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"RtoRref",RtoRref());
    bindcommand(name,"read_elements",read_elements(elements_file));
    bindcommand(name,"read_fem_coeff",read_fem_coeff(fem_coeff_file));    
    
    return -1;
}

void FEMFrame::RtoRref()
{
     int i;
     SHtoR();
     for(i=0;i<_NP;i++)
     {    
         _SRref[i] = _SR[i];
         _Rref[i]  = _R[i];
     }
}

void FEMFrame::Alloc()
{
    MDPARALLELFrame::Alloc();

    int size;
    size = _NP*allocmultiple;
    Realloc(_SRref,Vector3,size);
    Realloc(_Rref, Vector3,size);
    bindvar("SRref",_SRref,DOUBLE);
    bindvar("Rref",_Rref,DOUBLE);

    size = _NELE*_NNODE_PER_ELEMENT;
    Realloc(elements,int,size);
    Realloc(map23,int,2);
    bindvar("elements",elements,INT);
    bindvar("map23",map23,INT);

    for(int i=0;i<_NP;i++) _VSR[i].clear();
}   

void FEMFrame::Alloc_Elements()
{
    int size;
    size = _NELE*_NNODE_PER_ELEMENT;
    Realloc(elements,int,size);
    bindvar("elements",elements,DOUBLE);

}

void FEMFrame::Alloc_Element_Coeff()
{
    int size;
    size = _NINT_PER_ELEMENT;
    Realloc(gauss_weight,double,size);
    bindvar("gauss_weight",gauss_weight,DOUBLE);

    size = _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT;
    Realloc(dFdu,double,size);
    bindvar("dFdu",dFdu,DOUBLE);
}

int FEMFrame::read_elements(const char* fname)
{
    int iE, j, iN;
    char *buffer; char *pp, *q;
    int np, nframe;
    char extfname[MAXFILENAMELEN];
    double coeff;
    int readCharCount, ind;

    strcpy(extfname,fname);
    LFile::SubHomeDir(extfname,extfname);
    LFile::LoadToString(extfname,buffer,0);

    pp=buffer;

    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    sscanf(q,"%d",&_NELE);

    Alloc_Elements();
    
    for(iE=0;iE<_NELE;iE++)
    {
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
       //INFO_Printf("%s\n",q);
       for(j=0;j<_NNODE_PER_ELEMENT;j++)
       {
            sscanf(q,"%d%n",&iN,&readCharCount);
            q += readCharCount;
            ind=iE*_NNODE_PER_ELEMENT+j;
            elements[ind] = iN;
            //INFO_Printf("%d %d  ",ind,elements[ind]);    
       }
    }
    //INFO_Printf("\n");

    Free(buffer);
    return 0;
}

int FEMFrame::read_fem_coeff(const char* fname)
{
    int iA, j, k, m, iN;
    char *buffer; char *pp, *q;
    int np, nframe;
    char extfname[MAXFILENAMELEN];
    double coeff;
    int readCharCount, ind;

    strcpy(extfname,fname);
    LFile::SubHomeDir(extfname,extfname);
    LFile::LoadToString(extfname,buffer,0);

    pp=buffer;
    /* skip first line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    /* fetch a line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    sscanf(q,"%d %d %d",&_NDIM,&_NNODE_PER_ELEMENT,&_NINT_PER_ELEMENT);

    Alloc_Element_Coeff();
    
    /* skip third line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<_NDIM;j++)
    {
       /* skip lines for reference configuration */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    }

    /* skip line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
    {
       sscanf(q,"%lf%n",gauss_weight+iA,&readCharCount);
       q += readCharCount;
       //INFO_Printf("gauss_weight[%d] = %f\n",iA,gauss_weight[iA]);
    }

    /* skip line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<_NDIM;j++)
    {
       /* skip lines for Gauss point positions */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    }

    for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
    {
       /* skip a line for each Gauss point */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
       INFO_Printf("%s\n",q);

       for(j=0;j<_NDIM;j++)
       {
           for(k=0;k<_NDIM;k++)
           {
                /* skip lines for Gauss point positions */
                q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
                //INFO_Printf("%s\n",q);
                for(m=0;m<_NDIM;m++)
                {
                    for(iN=0;iN<_NNODE_PER_ELEMENT;iN++)
                    {
                        sscanf(q,"%lf%n",&coeff,&readCharCount);
                        q += readCharCount;
                        ind=(((iA*_NDIM+j)*_NDIM+k)*_NDIM+m)*_NNODE_PER_ELEMENT+iN;
                        dFdu[ind] = coeff;
                        //INFO_Printf("%d %f  ",ind,dFdu[ind]);    
                    }
                }
                //INFO_Printf("\n");
           }
       }
    }

    Free(buffer);
    return 0;
}

void FEMFrame::plot()
{
    MDPARALLELFrame::plot();

    /* draw covalent bonds between atoms */
    int iele, ipt, j, jpt;
    double L;
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    int r,g,b; double alpha; unsigned ce;
    Vector3 sri, srj, sij, ri, rj;

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=max(_H[0][0],_H[1][1]);
    L=max(L,_H[2][2])*.5;

    SHtoR();
    win->Lock();
    //win->Clear();
    
    for(iele=0;iele<_NELE;iele++)
    {
      if (_NDIM == 2)
        for(j=0;j<_NNODE_PER_ELEMENT;j++)
        {
             jpt=elements[iele*_NNODE_PER_ELEMENT+j];
             ipt=elements[iele*_NNODE_PER_ELEMENT+((j+1)%_NNODE_PER_ELEMENT)];
             sri=_SR[ipt];
             if(plot_map_pbc==1) sri.subint();
             ri = _H*sri;

             sij=_SR[jpt]-_SR[ipt];
             sij.subint();
             srj=sri+sij;
             rj = _H*srj;
                
             if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                     continue;

                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 
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
    }
    
    win->Unlock();
    win->Refresh();
}

void FEMFrame::fem_energy_force()
{
    /* no need of neighbor list */
    int i,iele,j,jpt,iA, in, ip, iq, ind;
    Vector3 dsj, drj;
    Matrix33 Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    double Eele, Eint, Jet, trace, detB, J2;

    /*
      Xiaohan: print out stress in the membrane
    */    
    Matrix33 E,I,pk2, temp1, temp2;    


    //    int map23[2] = {0, 2};

    DUMP("FEM");

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0; //put this back
    for(i=0;i<_NP;i++)
    {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      // if (species[i] == 0) 
      // 	continue;

      _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear();
    }
    _VIRIAL.clear();

    hinv = _H.inv();
    for(iele=0;iele<_NELE;iele++)
    {
        /* energy of this element */
        Eele = 0;
        for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
        {
            /* energy of this Gauss integration point */
            Eint = 0;
            /* deformation gradient */
            Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
            for(j=0;j<_NNODE_PER_ELEMENT;j++)
            {
                jpt=elements[iele*_NNODE_PER_ELEMENT+j];

                _SRref[jpt ] = hinv*_Rref[jpt];  /* put this into a separate function RrefHtoSref */
                dsj = _SR[jpt] - _SRref[jpt];
                dsj.subint();
                _H.multiply(dsj, drj);

                /* Add contribution from drj to F using dFdu */
                for(ip=0;ip<_NDIM;ip++)
                {
                   for(iq=0;iq<_NDIM;iq++)
                   {
                       for(in=0;in<_NDIM;in++)
                       {
                            ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
                            //Fdef[ip][iq] += dFdu[ind]*drj[in]; 
			    Fdef[ip][iq] += dFdu[ind]*drj[map23[in]]; 
                       }
                   }
                }
            }
            B = Fdef*Fdef.tran();
            C = Fdef.tran()*Fdef;
            Cinv = C.inv();
            FCinvtran = Fdef*Cinv.tran();
            Jet = Fdef.det();  J2 = Jet*Jet;
            detB = B.det(); 
            Eint = 0;


	    /*
	      Xiaohan: print out stress in the membrane
	     */
	    if (0){//iele == 3 || iele == 4 || iele == 11 || iele == 12 || iele == 19 || iele == 20 || iele == 27 || iele == 28 || iele == 35 || iele == 36 ) {
	      double lambda = 1.5; double mu = 1;
	      pk2.clear();
	      I[0][0]=I[1][1]=I[2][2]=1;
	      E = (C-I);
	      E *= 0.5;
	      temp1 = I;
	      temp1*=lambda*E.trace();
	      temp2 = E;
	      temp2 *= 2*mu;
	      pk2 += temp1 + temp2;
	      pk2= Fdef*pk2*Fdef.tran();

	      INFO_Printf("PK2 stress on Gaussian point %d of element %d\n", iA, iele);
	      INFO_Printf("%f %f %f %f %f %f\n", pk2[0][0], pk2[1][1], pk2[2][2], pk2[0][1], pk2[0][2], pk2[1][2]);
	    }
            /* Add contribution from F to Eint */ 
            if (_NDIM == 2)
            {
                trace = B[0][0] + B[1][1];
                Eint = 0.5*(trace + 1.0/detB - 3); /* multiply MU and V0 later */

                dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */

                //INFO_Printf("F = [ %f %f \n      %f %f ]\n",Fdef[0][0],Fdef[0][1],Fdef[1][0],Fdef[1][1]);
                //INFO_Printf("B = [ %f %f \n      %f %f ]\n",B[0][0],B[0][1],B[1][0],B[1][1]);
                //INFO_Printf("Jet =  %f   detB = %f  Eint = %f\n",Jet,detB,Eint);
                //INFO_Printf("dEdF = [ %f %f \n      %f %f ]\n",dEdF[0][0],dEdF[0][1],dEdF[1][0],dEdF[1][1]);
                /* Add contribution from drj to F using dFdu */
                for(j=0;j<_NNODE_PER_ELEMENT;j++)
                {
                    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
                    for(ip=0;ip<_NDIM;ip++)
                    {
                       for(iq=0;iq<_NDIM;iq++)
                       {
                           for(in=0;in<_NDIM;in++)
                           {
                                ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
                                if (fixed[jpt] == 0)
				  //_F[jpt][in] -= dEdF[ip][iq]*dFdu[ind]; 
				  _F[jpt][map23[in] ] -= dEdF[ip][iq]*dFdu[ind]; 
                           }
                       }
                    }
                }
            }

            /* Add contribution from Eint to Eele */
            Eele += gauss_weight[iA]*Eint;
        }
        _EPOT += Eele;
    }
    INFO_Printf("_EPOT = %lf\n", _EPOT);
}    

void FEMFrame::potential()
{
    fem_energy_force();
}

//void FEMFrame::potential_energyonly()
//{
//    fem_energy_only();
//}


#ifdef _TEST

/* Main Program Begins */
class FEMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

