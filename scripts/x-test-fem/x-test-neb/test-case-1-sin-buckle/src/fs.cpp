// Last Modified : Thu Apr 22 16:40:32 2010

// Last Modified : Thu Sep 29 2005
// Keonwook Kang, kwkang@stanford.edu

#include "fs.h"

void FSParam::initparser()
{
    bindvar("A",&const_A,DOUBLE);
    bindvar("d",&const_d,DOUBLE);
    bindvar("beta",&const_beta,DOUBLE);
    bindvar("c",&const_c,DOUBLE);
    bindvar("c0",&const_c0,DOUBLE);
    bindvar("c1",&const_c1,DOUBLE);
    bindvar("c2",&const_c2,DOUBLE);
        
    bindvar("alpha",&const_alpha,DOUBLE);
    bindvar("b0",&const_b0,DOUBLE);
    bindvar("B",&const_B,DOUBLE);
}
    
int FSParam::readfile(char *ppmfile)
{  /* could to be more elegant */

    FILE *f; char extname[200];
    LFile::SubHomeDir(ppmfile,extname);
    f=fopen(extname,"r");
    if(f==NULL)
    {
        FATAL("Potential file:"<<ppmfile<<" not found!");
        return -1;
    }
    parse(f); fclose(f);
    const_d_sq=const_d*const_d;
    const_c_sq=const_c*const_c;

    rcut=max(const_c,const_d);
    INFO("Rcut="<<rcut);
    rlist=1.1*rcut;
    skin=rlist - rcut;
    return 0;
}

void FSParam::print()
{
    char tmp[2000];
    sprintf(tmp,HIB"FS Potential Parameters"NOR
            "const_d=%f const_A=%f const_beta=%f\n"
            "const_c=%f const_c0=%f const_c1=%f const_c2=%f\n"
            "const_B=%f const_alpha=%f const_b0=%f\n",
            const_d,const_A,const_beta,const_c,const_c0,const_c1,
            const_c2,const_B,const_alpha,const_b0);
    INFO(tmp);
}
    
int FSParam::exec(const char *name)
{
    WARNING(HIC"EXEC   "NOR<<name<<"  (ignored)");
    return -1;
}


int FSParam::assignvar(int offset) 
{ /* copied from organizer.h */
    int s;
    switch(vartype[curn])
    {
    case(INT): s=read_buffer("%d",int); break;
    case(LONG): s=read_buffer("%ld",long); break;
    case(DOUBLE): s=read_buffer("%lf",double); break;
    case(STRING): s=read_buffer("%s",char); break;
    default: FATAL("unknown vartype ("<<vartype[curn]<<")");
    }
    if(s==0)WARNING("Expression syntax error: ("
                    <<varname[curn]<<" = "<<buffer<<"), variable "
                    <<varname[curn]<<" unchanged");
    switch(vartype[curn])
    {
    case(INT): INFO_buffer(ah,*,int); break;
    case(LONG): INFO_buffer(ah,*,long); break;
    case(DOUBLE): INFO_buffer(ah,*,double); break;
    case(STRING): INFO_buffer(ah,,char); break;
//        case(STRING): DUMP(varname[curn]<<"="<<((char *)varptr[curn]));break;
    default: WARNING("unknown vartype ("<<vartype[curn]<<")");
    }
    return (!s);
}

/* FSFrame Begins */
int FSFrame::readpot()
{
    int ret;
    ret=_FSPM.readfile(potfile);
    _RLIST=_FSPM.rlist;_SKIN=_FSPM.skin;
    return ret;
}
void FSFrame::initvars()
{
    MDPARALLELFrame::initvars();
   
    strcpy(incnfile,"../bcc.cn");
}
void FSFrame::initparser()
{
    MDPARALLELFrame::initparser();
}

int FSFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readpot",readpot());
#ifdef _PARALLEL
    bindcommand(name,"Broadcast_FS_Param",Broadcast_FS_Param());
#endif
    return -1;
}
 
void FSFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
    int size;
    size=_NP*allocmultiple;
    
    Realloc(rhoh,double,size);
    Realloc(af,double,size);

    bindvar("rhoh",rhoh,DOUBLE);
    bindvar("af",af,DOUBLE);
}    

void FSFrame::finnis_sinclair() 
{  /* Potential function (_SR, _H) --> (_F, _EPOT, _VIRIAL) */
    int i, j, jpt;
    //int symtype, symindex; /* symmetry operation */
    
    Vector3 sij,rij,uij,fij,fs;
    double r2ij,rrij,rhij,vij,temp0,temp1,temp2,ddtemp;
    double dr, dr2, dr3;
    
    DUMP("Finnis-Sinclair Mo");
    DUMP("NP="<<_NP);
    refreshneighborlist();

    _EPOT=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0;}
    //for(i=_NP0;i<_NP1;i++)
    for(i=0;i<_NP;i++)
        rhoh[i]=0; 
    _VIRIAL.clear();

        for(i=0;i<_NP;i++)
        {
            for(j=0;j<nn[i];j++)
            {
                jpt=nindex[i][j];
                if(i>=jpt) continue;
                sij=_SR[jpt]-_SR[i];
                sij.subint();
                rij=_H*sij;
                r2ij=rij.norm2();
                if(r2ij<_FSPM.const_d_sq)
                {
                    rrij=sqrt(r2ij);
                    dr=(rrij-_FSPM.const_d);
                    dr2=dr*dr; dr3=dr*dr2;
                    rhij=dr2+_FSPM.const_beta*dr3/_FSPM.const_d;
                    rhoh[i]+=rhij;
                    rhoh[jpt]+=rhij;
                }
            }
        }
    //for(i=_NP0;i<_NP1;i++)
    for(i=0;i<_NP;i++)
    {
        af[i]=sqrt(rhoh[i]);
    }
    //for(i=_NP0;i<_NP1;i++)
    for(i=0;i<_NP;i++)
    {
        if((fixedatomenergypartition==0)||(fixed[i]!=1))
        {
            _EPOT-=_FSPM.const_A*af[i];
            _EPOT_IND[i]-=_FSPM.const_A*af[i];
            _EPOT_RMV[i]-=_FSPM.const_A*af[i];
        }
    }

    //for(i=_NP0;i<_NP1;i++)
    for(i=0;i<_NP;i++)
    {
        if(fixedatomenergypartition)
            if(fixed[i]==1) continue;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(fixedatomenergypartition)
            {
                if(fixed[jpt]!=1)
                    if(i>=jpt) continue;
            }
            else
            {
                    if(i>=jpt) continue;
            }
                
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            rrij=sqrt(r2ij);
            dr=(rrij-_FSPM.const_d);
            dr2=dr*dr; dr3=dr*dr2;
                
            if(rrij<_FSPM.rcut)
            {
                uij=rij/rrij;
                vij=0;
                /* Two-body energy and forces */
                if(rrij<_FSPM.const_c)
                {
                    temp1=rrij-_FSPM.const_c;
                    temp2 = _FSPM.const_c0+
					_FSPM.const_c1*rrij+_FSPM.const_c2*r2ij;
                    ddtemp=temp1*temp1*temp2;
                    _EPOT_IND[i]+=ddtemp*.5;
                    _EPOT_IND[jpt]+=ddtemp*.5;

                    _EPOT_RMV[i]+=ddtemp;
                    _EPOT_RMV[jpt]+=ddtemp;

                    if(fixedatomenergypartition==0)
                    {
                        _EPOT+=ddtemp;
                        vij+=2*temp1*temp2+temp1*temp1*
                            (_FSPM.const_c1 + 2*_FSPM.const_c2*rrij);
                    }
                    else
                    {
                        if(fixed[jpt]!=1)
                        {
                            _EPOT+=ddtemp;
                            vij+=2*temp1*temp2+temp1*temp1*
                                (_FSPM.const_c1 + 2*_FSPM.const_c2*rrij);
                        }
                        else
                        {
                            _EPOT+=ddtemp*.5;
                            vij+=(2*temp1*temp2+temp1*temp1*
                                (_FSPM.const_c1 + 2*_FSPM.const_c2*rrij))*.5;
                            
                        }
                    }
                }
                /* Additional two-body term as core function */
                if(rrij<_FSPM.const_b0)
                {
                    temp0=_FSPM.const_b0-rrij;
                    temp1=temp0*temp0;
                    temp2=exp(-_FSPM.const_alpha*rrij);

                    ddtemp=_FSPM.const_B*temp0*temp1*temp2;
                    _EPOT_IND[i]+=ddtemp*.5;
                    _EPOT_IND[jpt]+=ddtemp*.5;

                    _EPOT_RMV[i]+=ddtemp;
                    _EPOT_RMV[jpt]+=ddtemp;
                    if(fixedatomenergypartition==0)
                    {
                        _EPOT+=ddtemp;
                        vij-=_FSPM.const_B*temp1*temp2*
                            (3+_FSPM.const_alpha*temp0);
                    }
                    else
                    {
                        if(fixed[jpt]!=1)
                        {
                            _EPOT+=ddtemp;
                            vij-=_FSPM.const_B*temp1*temp2*
                                (3+_FSPM.const_alpha*temp0);
                        }
                        else
                        {
                            _EPOT+=ddtemp*.5;
                            vij-=_FSPM.const_B*temp1*temp2*
                                (3+_FSPM.const_alpha*temp0)*.5;
                        }
                    }
                }
                /* Many-body energy and forces */
                if(rrij<_FSPM.const_d)
                {
                    if(fixedatomenergypartition==0)
                    {
                        vij-=_FSPM.const_A*
                            ( 1./af[i]+1./af[jpt] )/2*(
                                2* dr +
                                _FSPM.const_beta * 3 * dr2/_FSPM.const_d);
                    }
                    else
                    {
                        if(fixed[jpt]!=1)
                            vij-=_FSPM.const_A*
                                ( 1./af[i]+1./af[jpt] )/2*(
                                    2* dr +
                              _FSPM.const_beta * 3 * dr2/_FSPM.const_d);
                        else
                            vij-=_FSPM.const_A*
                                ( 1./af[i] )/2*(
                                    2* dr +
                              _FSPM.const_beta * 3 * dr2/_FSPM.const_d);
                    }
                }
                fij=uij*vij;
                _F[i]+=fij;
                _F[jpt]-=fij;

                _VIRIAL.addnvv(-vij*rrij,uij,uij);
            }
        }
    }
    /* zero out the forces on fixed atoms */
    /* no longer doing this */

//    INFO_Printf("EPOT=%20.14e (end of finnis_sinclair)\n",_EPOT);

#if 0 /* disable symmetry operation on force */
//    INFO_Printf("EPOT=%20.14e (before accumulation)\n",_EPOT);
//    DUMP("FS Collect");
    /* _EPOT, _F, _EPOT_IND, _VIRIAL will be summed */
    /* apply symmetry operation on force */
    /* also need to symmetrize position and velocity in calcprop */
/*    symtype = (int)symmetryspec[0];
    symindex = (int)symmetryspec[1]; symindex--;
    if((symtype>0)&&(symindex>=0))
    {
        if(symtype==1)
        {
//            INFO("symmetrize force type="<<symtype<<", index="<<symindex);
            for(i=0;i<_NP;i++)
            {
                if(i>image[i]) continue;
//                INFO_Printf("_F[%d]=(%e %e %e) _F[%d]=(%e %e %e)\n",
//                            i,_shmF[i].x,_shmF[i].y,_shmF[i].z,
//                            image[i],_shmF[image[i]].x,_shmF[image[i]].y,_shmF[image[i]].z);
                if(i==image[i])
                {
                    _shmF[i][symindex]=0;
                    continue;
                }
                fs=_shmF[image[i]]+_shmF[i];
                fs[symindex]=_shmF[image[i]][symindex]-_shmF[i][symindex];
                fs*=0.5;
                _shmF[i]=fs;
                _shmF[image[i]]=fs;
                _shmF[image[i]][symindex]-=fs[symindex]*2;
//                INFO_Printf("_F[%d]=(%e %e %e) _F[%d]=(%e %e %e)\n",
//                            i,_shmF[i].x,_shmF[i].y,_shmF[i].z,
//                            image[i],_shmF[image[i]].x,_shmF[image[i]].y,_shmF[image[i]].z);
            }
        }
        else
        {
            FATAL("symmetry type("<<symtype<<", "<<symindex<<") not implemented");
        }
    }

    StrMemCpy(&(_EPOT),_shmEPOTptr,sizeof(double));
    StrMemCpy((double *)_F,(double *)_shmF,3*_NP*sizeof(double));
    StrMemCpy(_EPOT_IND,_shmEPOT_IND,_NP*sizeof(double));
    StrMemCpy((double *)(_VIRIAL[0]),
              (double *)((*_shmVIRIALptr)[0]),9*sizeof(double));
*/
    /* debug print out EPOT_IND */
//    temp0=0;
//    for(i=0;i<10;i++)
//    {
//        INFO_Printf("EPOT_IND[%d]=%20.14f\n",i,_EPOT_IND[i]);
//        temp0+=_EPOT_IND[i];
//    }
//    INFO_Printf("EPOT=%20.14e  sum{EPOT_IND}=%20.14e\n",_EPOT,temp0);
#endif
}

void FSFrame::finnis_sinclair_energyonly() 
{  /* Potential function (_SR, _H) --> (_F, _EPOT, _VIRIAL) */
    int i, j, jpt;
    
    Vector3 sij,rij,uij,fij,fs;
    double r2ij,rrij,rhij,vij,temp0,temp1,temp2,ddtemp;
    double dr, dr2, dr3;
    
    refreshneighborlist();

    _EPOT=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0;}
    for(i=0;i<_NP;i++)
        rhoh[i]=0; 
    _VIRIAL.clear();

        for(i=0;i<_NP;i++)
        {
            for(j=0;j<nn[i];j++)
            {
                jpt=nindex[i][j];
                if(i>=jpt) continue;
                sij=_SR[jpt]-_SR[i];
                sij.subint();
                rij=_H*sij;
                r2ij=rij.norm2();
                if(r2ij<_FSPM.const_d_sq)
                {
                    rrij=sqrt(r2ij);
                    dr=(rrij-_FSPM.const_d);
                    dr2=dr*dr; dr3=dr*dr2;
                    rhij=dr2+_FSPM.const_beta*dr3/_FSPM.const_d;
                    rhoh[i]+=rhij;
                    rhoh[jpt]+=rhij;
                }
            }
        }
    for(i=0;i<_NP;i++)
    {
        af[i]=sqrt(rhoh[i]);
    }
    for(i=0;i<_NP;i++)
    {
        if((fixedatomenergypartition==0)||(fixed[i]!=1))
        {
            _EPOT-=_FSPM.const_A*af[i];
            _EPOT_IND[i]-=_FSPM.const_A*af[i];
            _EPOT_RMV[i]-=_FSPM.const_A*af[i];
        }
    }

    for(i=0;i<_NP;i++)
    {
        if(fixedatomenergypartition)
            if(fixed[i]==1) continue;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(fixedatomenergypartition)
            {
                if(fixed[jpt]!=1)
                    if(i>=jpt) continue;
            }
            else
            {
                    if(i>=jpt) continue;
            }
                
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r2ij=rij.norm2();
            rrij=sqrt(r2ij);
            dr=(rrij-_FSPM.const_d);
            dr2=dr*dr; dr3=dr*dr2;
                
            if(rrij<_FSPM.rcut)
            {
                uij=rij/rrij;
                vij=0;
                /* Two-body energy and forces */
                if(rrij<_FSPM.const_c)
                {
                    temp1=rrij-_FSPM.const_c;
                    temp2 = _FSPM.const_c0+
					_FSPM.const_c1*rrij+_FSPM.const_c2*r2ij;
                    ddtemp=temp1*temp1*temp2;
                    _EPOT_IND[i]+=ddtemp*.5;
                    _EPOT_IND[jpt]+=ddtemp*.5;

                    _EPOT_RMV[i]+=ddtemp;
                    _EPOT_RMV[jpt]+=ddtemp;

                    if(fixedatomenergypartition==0)
                    {
                        _EPOT+=ddtemp;
                    }
                    else
                    {
                        if(fixed[jpt]!=1)
                        {
                            _EPOT+=ddtemp;
                        }
                        else
                        {
                            _EPOT+=ddtemp*.5;
                        }
                    }
                }
                /* Additional two-body term as core function */
                if(rrij<_FSPM.const_b0)
                {
                    temp0=_FSPM.const_b0-rrij;
                    temp1=temp0*temp0;
                    temp2=exp(-_FSPM.const_alpha*rrij);

                    ddtemp=_FSPM.const_B*temp0*temp1*temp2;
                    _EPOT_IND[i]+=ddtemp*.5;
                    _EPOT_IND[jpt]+=ddtemp*.5;

                    _EPOT_RMV[i]+=ddtemp;
                    _EPOT_RMV[jpt]+=ddtemp;

                    if(fixedatomenergypartition==0)
                    {
                        _EPOT+=ddtemp;
                    }
                    else
                    {
                        if(fixed[jpt]!=1)
                        {
                            _EPOT+=ddtemp;
                        }
                        else
                        {
                            _EPOT+=ddtemp*.5;
                        }
                    }
                }
            }
        }
    }
}



void FSFrame::potential()
{
    finnis_sinclair();
}

void FSFrame::potential_energyonly()
{
    finnis_sinclair_energyonly();
}

#ifdef _PARALLEL
void FSFrame::Broadcast_FS_Param()
{
    double *buf_double;
    int nparam, i, j;

    /* master asking slaves to call the same function */
    if(myDomain==0) Master_to_Slave("Broadcast_FS_Param");    

    /* packing parameters */
    nparam = 17;

    buf_double = (double *) malloc(nparam*sizeof(double));
    if(myDomain==0)
    {
        buf_double[0] = _FSPM.const_d;
        buf_double[1] = _FSPM.const_A;
        buf_double[2] = _FSPM.const_beta;
        buf_double[3] = _FSPM.const_c;
        buf_double[4] = _FSPM.const_c0;
        buf_double[5] = _FSPM.const_c1;
        buf_double[6] = _FSPM.const_c2;
        buf_double[7] = _FSPM.const_d_sq;
        buf_double[8] = _FSPM.const_c_sq;
        buf_double[9] = _FSPM.rcut;
        buf_double[10] = _FSPM.const_B;
        buf_double[11] = _FSPM.const_alpha;
        buf_double[12] = _FSPM.const_b0;
        buf_double[13]  = _RLIST;
        buf_double[14]  = _SKIN;
        buf_double[15] = (double) _NNM;
        buf_double[16] = (double) _NIC;
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
        _FSPM.const_d = buf_double[0];
        _FSPM.const_A = buf_double[1];
        _FSPM.const_beta = buf_double[2];
        _FSPM.const_c = buf_double[3];
        _FSPM.const_c0 = buf_double[4];
        _FSPM.const_c1 = buf_double[5];
        _FSPM.const_c2 = buf_double[6];
        _FSPM.const_d_sq = buf_double[7];
        _FSPM.const_c_sq = buf_double[8];
        _FSPM.rcut = buf_double[9];
        _FSPM.const_B = buf_double[10];
        _FSPM.const_alpha = buf_double[11];
        _FSPM.const_b0 = buf_double[12];
        _RLIST = buf_double[13];
        _SKIN  = buf_double[14]; 
        _NNM = (int) buf_double[15];
        _NIC = (int) buf_double[16];
    }

    free(buf_double);

}
#endif//_PARALLEL

#ifdef _TEST

/* Main Program Begins */
class FSFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
