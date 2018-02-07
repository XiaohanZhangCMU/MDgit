/*
  ysz.cpp
  by Eunseok Lee euniv@stanford.edu
  Last Modified : Wed Jan 13 14:26:26 2010
  by Eunseok Lee euniv@stanford.edu
  Last Modified : Sep 2  
  Last Modified : 4/12/2009
  
  FUNCTION  :  KMC simulation of YSZ
*/

#include "ysz9.h"


void YSZFrame::initparser()
{
    KMCFrame::initparser();

    /* input */
    bindvar("linklistfile",linklistfile,STRING);
    bindvar("ebarrierfile",ebarrierfile,STRING);
    bindvar("potfile",potfile,STRING);
    bindvar("catpotfile",catpotfile,STRING);
    bindvar("catspeciesfile",catspeciesfile,STRING);
    bindvar("nnYlistfile",nnYlistfile,STRING);
    bindvar("snYlistfile",snYlistfile,STRING);
    
    bindvar("savecorr",&savecorr,INT);
    bindvar("savecorrfreq",&savecorrfreq,INT);
    bindvar("corrfiles",corrfiles,STRING);
    bindvar("corrfilel",corrfilel,STRING);
    bindvar("saveac",&saveac,INT);
    bindvar("saveacfreq",&saveacfreq,INT);
    bindvar("acfile",acfile,STRING);
    bindvar("savediff",&savediff,INT);
    bindvar("savedifffreq",&savedifffreq,INT);
    bindvar("difffile",difffile,STRING);
    bindvar("savehist",&savehist,INT);
    bindvar("savehistfreq",&savehistfreq,INT);   
    bindvar("histfile",histfile,STRING);
    bindvar("saveatomeye",&saveatomeye,INT);
    bindvar("saveatomeyefreq",&saveatomeyefreq,INT);
    bindvar("atomeyefile",atomeyefile,STRING);
    bindvar("savecontdata",&savecontdata,INT);
    bindvar("savecontfreq",&savecontfreq,INT);
    
    bindvar("extU0",&extU0,DOUBLE);
    bindvar("extw0",&extw0,DOUBLE);
    bindvar("timelimit",&timelimit,DOUBLE);

    bindvar("extdt",&extdt,DOUBLE);
    bindvar("ncycle",&ncycle,INT);
    
    bindvar("num_jump",&num_jump,INT);
    bindvar("num_loop",&num_loop,INT);
    bindvar("ncorrsteps",&ncorrsteps,INT);
    bindvar("ncorrbins",&ncorrbins,INT);
    bindvar("tdurs",&tdurs,DOUBLE);
    bindvar("corrdts",&corrdts,DOUBLE);
    bindvar("ncorrstepl",&ncorrstepl,INT);
    bindvar("ncorrbinl",&ncorrbinl,INT);
    bindvar("tdurl",&tdurl,DOUBLE);
    bindvar("corrdtl",&corrdtl,DOUBLE);

    bindvar("CTL_Ebarrier",&CTL_Ebarrier,INT);
    bindvar("CTL_CoulombE",&CTL_CoulombE,INT);
    bindvar("CTL_CatPot",&CTL_CatPot,INT);

    bindvar("EnnY",&EnnY,DOUBLE);
    bindvar("EsnY",&EsnY,DOUBLE);
    bindvar("Ebconst",&Ebconst,DOUBLE);
    
    bindvar("_T",&_T,DOUBLE);
    bindvar("epsilon_r",&epsilon_r,DOUBLE);

    bindvar("nAA",&nAA,INT);
    bindvar("nCC",&nCC,INT);
    bindvar("nCA",&nCA,INT);
    bindvar("nCCC",&nCCC,INT);
    bindvar("nCCA",&nCCA,INT);
    bindvar("nCAA",&nCAA,INT);
    bindvar("nAAA",&nAAA,INT);
    bindvar("neci",&neci,INT);
    bindvar("cluster_rcut",&cluster_rcut,DOUBLE);
    bindvar("cluster_dcut",&cluster_dcut,DOUBLE);
    bindvar("cemfilesdir",cemfilesdir,STRING);
    bindvar("clusterCC_file",clusterCC_file,STRING);
    bindvar("clusterAA_file",clusterAA_file,STRING);
    bindvar("clusterCA_file",clusterCA_file,STRING);
    bindvar("clusterCCC_file",clusterCCC_file,STRING);
    bindvar("clusterCCA_file",clusterCCA_file,STRING);
    bindvar("clusterCAA_file",clusterCAA_file,STRING);
    bindvar("clusterAAA_file",clusterAAA_file,STRING);
    bindvar("eci_opt_file",eci_opt_file,STRING);
    bindvar("maxsize_cluster_buffer",&maxsize_cluster_buffer,INT);

    bindvar("ntransient",&ntransient,INT);
    bindvar("ebmaxratio",&ebmaxratio,DOUBLE);
    bindvar("ebaddratio",&ebaddratio,DOUBLE);
    bindvar("ebcorrtable_file",ebcorrtable_file,STRING);
    bindvar("ebxcorrtable_file",ebxcorrtable_file,STRING);
    bindvar("nebcorrtable",&nebcorrtable,INT);

    bindvar("EbZZ_Car",&EbZZ_Car,DOUBLE);
    bindvar("EbZY_Car",&EbZY_Car,DOUBLE);
    bindvar("EbYY_Car",&EbYY_Car,DOUBLE);
}

int YSZFrame::exec(const char *name)
{
    if(KMCFrame::exec(name)==0) return 0;
    bindcommand(name,"AllocLinkList",AllocLinkList());
    bindcommand(name,"BuildLinkList",BuildLinkList());
    bindcommand(name,"ReadLinkList",ReadLinkList());    
    bindcommand(name,"BuildEnergyBarrierList",BuildEnergyBarrierList());
    bindcommand(name,"ReadEnergyBarrierList",ReadEnergyBarrierList());
    bindcommand(name,"BuildCoulombEList",BuildCoulombEList());
    bindcommand(name,"BuildCatPotList",BuildCatPotList());
    bindcommand(name,"BuildCatSpeciesList",BuildCatSpeciesList());
    bindcommand(name,"BuildNNYList",BuildNNYList());
    bindcommand(name,"BuildSNYList",BuildSNYList());
    bindcommand(name,"ReadCoulombEList",ReadCoulombEList());
    bindcommand(name,"ReadCatPotList",ReadCatPotList());
    bindcommand(name,"ReadCatSpeciesList",ReadCatSpeciesList());
    bindcommand(name,"ReadNNYList",ReadNNYList());
    bindcommand(name,"ReadSNYList",ReadSNYList());
    bindcommand(name,"CreateConfig_Vacancy_Only",CreateConfig_Vacancy_Only());
    bindcommand(name,"kmcrun",kmcrun());
    bindcommand(name,"kmcruneq",kmcruneq());
    bindcommand(name,"kmcrunAC",kmcrunAC());
    bindcommand(name,"kmcrun_continued",kmcrun_continued());

    bindcommand(name,"PrintLinkList",PrintLinkList());
    bindcommand(name,"PrintCationList",PrintCationList());
    bindcommand(name,"PrintAnionList",PrintAnionList());

    bindcommand(name,"Add_Yttrium",Add_Yttrium());
    bindcommand(name,"Add_Yttrium_Manually",Add_Yttrium_Manually());
    bindcommand(name,"Add_Yttrium_Random",Add_Yttrium_Random());
    bindcommand(name,"Add_Yttrium_On_Surface",Add_Yttrium_On_Surface());
    bindcommand(name,"Add_Yttrium_On_Surface_X",Add_Yttrium_On_Surface_X());
    bindcommand(name,"Add_Vacancy", Add_Vacancy());
    bindcommand(name,"Add_Vacancy_On_Surface", Add_Vacancy_On_Surface());

    bindcommand(name,"opencorrfiles",cfs.open(corrfiles));
    bindcommand(name,"opencorrfilel",cfl.open(corrfilel));
    bindcommand(name,"openacfile",af.open(acfile));
    bindcommand(name,"opendifffile",df.open(difffile));
    bindcommand(name,"openhistfile",hf.open(histfile));
    bindcommand(name,"openatomeyefile",mf.open(atomeyefile));
    bindcommand(name,"CorrFcnX",CorrFcnX(DOUBLE));
    bindcommand(name,"CorrFcnXn",CorrFcnXn(INT,INT));
    bindcommand(name,"set_corr_filecounter",set_corr_filecounter(input[0]));
    bindcommand(name,"set_diff_filecounter",set_diff_filecounter(input[0]));
    bindcommand(name,"set_hist_filecounter",set_hist_filecounter(input[0]));
    bindcommand(name,"set_atomeye_filecounter",set_atomeye_filecounter(input[0]));

    bindcommand(name,"ReadClusters",ReadClusters());
    
    return -1;
}

void YSZFrame::initvars()
{
    strcpy(corrfiles,"corrs.dat");
    strcpy(corrfilel,"corrl.dat");
    strcpy(acfile,"ac.dat");
    strcpy(difffile,"diff.dat");
    strcpy(histfile,"hist.dat");
    strcpy(atomeyefile,"atomeye.dat");
    
    KMCFrame::initvars();
}

void YSZFrame::Alloc()
{
    KMCFrame::Alloc();
}

void YSZFrame::AllocLinkList()
{
    NCation = _NP / 3;
    NAnion  = _NP - NCation;

    NCellX = (int) (latticesize[0][3]*4);
    NCellY = (int) (latticesize[1][3]*4);
    NCellZ = (int) (latticesize[2][3]*4);

    INFO_Printf("AllocLinkList: NP=%d NCation=%d NAnion=%d\n", _NP,NCation,NAnion);
    INFO_Printf("               NCellX = %d  NCellY = %d  NCellZ = %d\n",NCellX,NCellY,NCellZ);
    
    if(NAnion!=NCellX*NCellY*NCellZ/8)
    {
        FATAL("NAnion = "<<NAnion<<" NCellX = "<<NCellX
              <<" NCellY = "<<NCellY<<" NCellZ = "<<NCellZ
              <<"NAnion != NCellX * NCellY * NCellZ / 8");
    }

    Realloc(_RI,int,_NP*3);
    Realloc(linklist,int,_NP*LINKLEN);
    Realloc(anionlist,int,NAnion);
    Realloc(cationlist,int,NCation);
    Realloc(vaclist,int,NAnion);    // resize to num_vac later
    Realloc(CoulombE,double,NCellX*NCellY*NCellZ);
    Realloc(CatPot,double,_NP);
    Realloc(Ebarrier,double,_NP*NUMJUMPDIR);
    Realloc(EbarrierN,double,_NP*NUMJUMPDIR);
    
    memset(_RI,0,sizeof(int)*_NP*3);
    memset(linklist,0,sizeof(int)*_NP*LINKLEN);
    memset(anionlist,0,sizeof(int)*NAnion);
    memset(cationlist,0,sizeof(int)*NCation);
    memset(vaclist,0,sizeof(int)*NCation);
    memset(CoulombE,0,sizeof(double)*NCellX*NCellY*NCellZ);
    memset(CatPot,0,sizeof(double)*_NP);
    memset(Ebarrier,0,sizeof(double)*_NP*NUMJUMPDIR); /* NAnion=NCellX*NCellY*NCellZ/8 */
    memset(EbarrierN,0,sizeof(double)*_NP*NUMJUMPDIR);   /* _NP by 6 */
    
    bindvar("_RI",_RI,INT);
    bindvar("linklist",linklist,INT);
    bindvar("anionlist",anionlist,INT);
    bindvar("cationlist",cationlist,INT);
    bindvar("CoulombE",CoulombE,DOUBLE);    
    bindvar("CatPot",CatPot,DOUBLE);    
    bindvar("Ebarrier",Ebarrier,DOUBLE);    
    bindvar("EbarrierN",EbarrierN,DOUBLE);    
/*    bindvar("corr",corr,DOUBLE);    
    bindvar("corrtime",corrtime,DOUBLE);    
    bindvar("corrdx",corrdx,DOUBLE);    */
}

    

void YSZFrame::BuildLinkList()
{
    int ipt, jpt, j, n, ncation, nanion, nx, ny, nz;
    double r, rnn, rsn, rcut1, rcut2;
    Vector3 dr, ds;
    Matrix33 hinv;
    FILE *fp;
    
    INFO_Printf("buildlinklist\n");

    /* nearest neighbor */
    rnn = sqrt(3)*0.25*latticeconst[0];  rcut1 = rnn*1.05;
    rsn = 0.5*latticeconst[0];           rcut2 = rsn*1.05;
    INFO_Printf("latticeconst = %f  rnn=%f  rsn = %f\n",latticeconst[0],rnn,rsn);

    INFO_Printf("NCellX = %d  NCellY = %d  NCellZ = %d\n",NCellX,NCellY,NCellZ);

    if(NAnion!=NCellX*NCellY*NCellZ/8)
    {
        FATAL("NAnion = "<<NAnion<<" NCellX = "<<NCellX
              <<" NCellY = "<<NCellY<<" NCellZ = "<<NCellZ
              <<"NAnion != NCellX * NCellY * NCellZ / 8");
    }
    
    refreshneighborlist();
    memset(linklist,0,sizeof(int)*_NP*LINKLEN);

    hinv = _H.inv();

    /* build anionlist and cationlist */
    ncation = nanion = 0;
    for(ipt=0;ipt<_NP;ipt++)
    {
        /* going through all atoms
         * whenever see a cation, append id to cationlist
         * the sequence of cation sites is not ordered
         */
        if((species[ipt]==0)||(species[ipt]==2)) /* Zr or Y */
        {
            cationlist[ncation]=ipt;
            ncation++;
        }

        /* whenever see an anion, insert id in the 
         * the sequence of anion sites is ordered by (x,y,z)
         * only works with perfect crystal created by makecrystal
         */
        if((species[ipt]==1)||(species[ipt]==3)) /* O or Vo */
        {
            nx = (int) round((_SR[ipt].x+0.5)*(NCellX/2)-0.5); /* nx goes from 0 to NCellX/2-1 */
            ny = (int) round((_SR[ipt].y+0.5)*(NCellY/2)-0.5); /* ny goes from 0 to NCellY/2-1 */
            nz = (int) round((_SR[ipt].z+0.5)*(NCellZ/2)-0.5); /* nz goes from 0 to NCellZ/2-1 */
            n = nx*(NCellY/2)*(NCellZ/2)+ny*(NCellZ/2)+nz;
            anionlist[n]=ipt;
            nanion++;
        }        
    }
    if(ncation!=NCation)
        ERROR("BuildLinkList: ncation="<<ncation<<" NCation="<<NCation);
    if(nanion!=NAnion)
        ERROR("BuildLinkList: nanion="<<nanion<<" NAnion="<<NAnion);    
    
    /* find cell type for all Anion sites and identify Cation neighbors */
    for(ipt=0;ipt<_NP;ipt++)
    {
        /* species[ipt]: 0 - Zr, 1 - O */
        if(species[ipt]!=1) continue; /* only consider Oxygen sites */

        /* find all atoms within the cut-off radius rcut1 */
        n=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            ds = _SR[jpt]-_SR[ipt]; /* applying PBC */
            ds.subint();
            dr = _H*ds;
            r = dr.norm(); /* distance between atom ipt and jpt */
            if(r<=rcut1)   /* rcut is 1.05 of nearest neighbor distance (between Zr and O) */
            {
                if(species[jpt]!=0) /* jpt must be Zr site */
                    ERROR("nearest neighbor not a cation site, consider change latticeconst");
                n++;
                if(n>4)
                    ERROR("nearest neighbor more than 4, consider change latticeconst");

                /* format of each Anion line in the link list:
                 *  0    oxygen(0)/vacancy(1)
                 *  1-4  nearest cation sites,
                 *  5-10 nearest anion sites,
                 *  11   subcell type
                 */
                
                if(ds.x*ds.y*ds.z<0)
                {
                    if(linklist[ipt*LINKLEN+11]==0)
                        linklist[ipt*LINKLEN+11] = 1;     /* anion cell type 1 */
                    else if(linklist[ipt*LINKLEN+11]!=1)  /* all four neighbors should claim the same cell type */
                        ERROR("anion cell type inconsistent!");

                    if((ds.x<0)&&(ds.y<0)&&(ds.z<0))      /* neighbor no. 1 (---) */
                        linklist[ipt*LINKLEN+1] = jpt;
                    else if((ds.x<0)&&(ds.y>0)&&(ds.z>0)) /* neighbor no. 2 (-++) */
                        linklist[ipt*LINKLEN+2] = jpt;
                    else if((ds.x>0)&&(ds.y<0)&&(ds.z>0)) /* neighbor no. 3 (+-+) */
                        linklist[ipt*LINKLEN+3] = jpt;
                    else if((ds.x>0)&&(ds.y>0)&&(ds.z<0)) /* neighbor no. 4 (++-) */
                        linklist[ipt*LINKLEN+4] = jpt;
                }
                else
                    if(linklist[ipt*LINKLEN+11]==0)
                        linklist[ipt*LINKLEN+11] = 2;      /* anion cell type 2 */
                    else if(linklist[ipt*LINKLEN+11]!=2)   /* all four neighbors should claim the same cell type */
                        ERROR("anion cell type inconsistent!");

                    if((ds.x<0)&&(ds.y<0)&&(ds.z>0))      /* neighbor no. 1 (--+) */
                        linklist[ipt*LINKLEN+1] = jpt;
                    else if((ds.x<0)&&(ds.y>0)&&(ds.z<0)) /* neighbor no. 2 (-+-) */
                        linklist[ipt*LINKLEN+2] = jpt;
                    else if((ds.x>0)&&(ds.y<0)&&(ds.z<0)) /* neighbor no. 3 (+--) */
                        linklist[ipt*LINKLEN+3] = jpt;
                    else if((ds.x>0)&&(ds.y>0)&&(ds.z>0)) /* neighbor no. 4 (+++) */
                        linklist[ipt*LINKLEN+4] = jpt;                    
            }
        }
        if(n!=4)
            ERROR("number of nearest neighbor not equal to 4");        
    }

    /* find Anion neighbors of each Anion */
    for(ipt=0;ipt<_NP;ipt++)
    {
        /* species[ipt]: 0 - Zr, 1 - O */
        if(species[ipt]!=1) continue; /* only consider Oxygen sites */

        /* find all atoms within the cut-off radius rcut */
        n=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            ds = _SR[jpt]-_SR[ipt];
            ds.subint();
            dr = _H*ds;
            r = dr.norm();
            if((r<=rcut2)&&(r>rcut1)) /* rcut2 is 1.05 of second nearest neighbor distance */
            {
                if(species[jpt]!=1) /* jpt must be Oxygen site */
                    ERROR("2nd nearest neighbor not a cation site, consider change latticeconst");
                n++;
                if(n>6)
                    ERROR("2nd nearest neighbor more than 6, consider change latticeconst");

                /* format of each Anion line in the link list:
                 *  0    oxygen(0)/vacancy(1)
                 *  1-4  nearest cation sites,
                 *  5-10 nearest anion sites,
                 *  11   subcell type
                 */
                
                if((fabs(ds.x)>fabs(ds.y))&&(fabs(ds.x)>fabs(ds.y)))
                {
                    if(ds.x>0) /* Anion neighbor no. 1 (+x) */
                        linklist[ipt*LINKLEN+5] = jpt;
                    else       /* Anion neighbor no. 2 (-x) */
                        linklist[ipt*LINKLEN+6] = jpt;
                }
                if((fabs(ds.y)>fabs(ds.x))&&(fabs(ds.y)>fabs(ds.z)))
                {
                    if(ds.y>0) /* Anion neighbor no. 3 (+y) */
                        linklist[ipt*LINKLEN+7] = jpt;
                    else       /* Anion neighbor no. 4 (-y) */
                        linklist[ipt*LINKLEN+8] = jpt;
                }
                if((fabs(ds.z)>fabs(ds.x))&&(fabs(ds.z)>fabs(ds.y)))
                {
                    if(ds.z>0) /* Anion neighbor no. 5 (+z) */
                        linklist[ipt*LINKLEN+9] = jpt;
                    else       /* Anion neighbor no. 6 (-z) */
                        linklist[ipt*LINKLEN+10]= jpt;
                }
            }
        }
        if(n!=6)
            ERROR("number of nearest neighbor not equal to 6");
    }

    /* loop through all Cation sites */
    for(ipt=0;ipt<_NP;ipt++)
    {
        /* species[ipt]: 0 - Zr, 1 - O */
        if(species[ipt]!=0) continue; /* only consider Zr sites */

        /* find all atoms within the cut-off radius rcut1 */
        n=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            ds = _SR[jpt]-_SR[ipt];
            ds.subint();
            dr = _H*ds;
            r = dr.norm();
            if(r<=rcut1) /* rcut is 1.05 of nearest neighbor distance */
            {
                if(species[jpt]!=1) /* jpt must be O site */
                    ERROR("nearest neighbor not a Anion site, consider change latticeconst");
                n++;
                if(n>8)
                    ERROR("nearest neighbor more than 8, consider change latticeconst");
        
                /* format of each Cation line in the link list:
                 *  0   Zr(0)/Y(1),
                 *  1-8 nearest Anion neighbors
                 */
                if     ((ds.x>0)&&(ds.y>0)&&(ds.z>0)) /* neighbor no. 1 (+x,+y,+z) */
                    linklist[ipt*LINKLEN+1] = jpt;
                else if((ds.x<0)&&(ds.y>0)&&(ds.z>0)) /* neighbor no. 2 (-x,+y,+z) */
                    linklist[ipt*LINKLEN+2] = jpt;
                else if((ds.x<0)&&(ds.y<0)&&(ds.z>0)) /* neighbor no. 3 (-x,-y,+z) */
                    linklist[ipt*LINKLEN+3] = jpt;
                else if((ds.x>0)&&(ds.y<0)&&(ds.z>0)) /* neighbor no. 4 (+x,-y,-z) */
                    linklist[ipt*LINKLEN+4] = jpt;
                else if((ds.x>0)&&(ds.y>0)&&(ds.z<0)) /* neighbor no. 5 (+x,+y,-z) */
                    linklist[ipt*LINKLEN+5] = jpt;
                else if((ds.x<0)&&(ds.y>0)&&(ds.z<0)) /* neighbor no. 6 (-x,+y,-z) */
                    linklist[ipt*LINKLEN+6] = jpt;
                else if((ds.x<0)&&(ds.y<0)&&(ds.z<0)) /* neighbor no. 7 (-x,-y,-z) */
                    linklist[ipt*LINKLEN+7] = jpt;
                else if((ds.x>0)&&(ds.y<0)&&(ds.z<0)) /* neighbor no. 8 (+x,-y,-z) */
                    linklist[ipt*LINKLEN+8] = jpt;
                
            }
        }
        if(n!=8)
            ERROR("number of nearest neighbor not equal to 8");        
    }
    
    for (ipt=0;ipt<_NP;ipt++) {
        *(_RI + 3*ipt + 0) = (int) round((_R[ipt].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*ipt + 1) = (int) round((_R[ipt].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*ipt + 2) = (int) round((_R[ipt].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    /* save Link list into a file */
    fp=fopen("LinkList.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildLinkList: file open error");
        return;
    }
    
    for(ipt=0;ipt<_NP;ipt++)
    {
        fprintf(fp,"%6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n",
                    linklist[ipt*LINKLEN+0],
                    linklist[ipt*LINKLEN+1],
                    linklist[ipt*LINKLEN+2],
                    linklist[ipt*LINKLEN+3],
                    linklist[ipt*LINKLEN+4],
                    linklist[ipt*LINKLEN+5],
                    linklist[ipt*LINKLEN+6],
                    linklist[ipt*LINKLEN+7],
                    linklist[ipt*LINKLEN+8],
                    linklist[ipt*LINKLEN+9],
                    linklist[ipt*LINKLEN+10],                        
                    linklist[ipt*LINKLEN+11]);
    }
    fclose(fp);
    INFO_Printf("written into LinkList.dat\n");
    
    /* Future KMC simulations can read this file */
}

int YSZFrame::ReadLinkList()
{
    FILE *fp;
    char string[500];
    int ipt;
    
    LFile::SubHomeDir(linklistfile,linklistfile);
    
    fp=fopen(linklistfile,"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadLinkList file ("<<linklistfile<<") not found!");
        return -1;
    }

    //AllocLinkList();
    
    for(ipt=0;ipt<_NP;ipt++)
    {
        fgets(string,500,fp);
        sscanf(string,"%d %d %d %d %d %d %d %d %d %d %d %d",
               & linklist[ipt*LINKLEN+0],
               & linklist[ipt*LINKLEN+1],
               & linklist[ipt*LINKLEN+2],
               & linklist[ipt*LINKLEN+3],
               & linklist[ipt*LINKLEN+4],
               & linklist[ipt*LINKLEN+5],
               & linklist[ipt*LINKLEN+6],
               & linklist[ipt*LINKLEN+7],
               & linklist[ipt*LINKLEN+8],
               & linklist[ipt*LINKLEN+9],
               & linklist[ipt*LINKLEN+10],                        
               & linklist[ipt*LINKLEN+11]);
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadLinkList done\n");
    return 0;
}

void YSZFrame::PrintLinkList()
{
    int n, i, ipt;
    n = (int)input[0];
    if(n<=0)
        ERROR("PrintLinkList n="<<n);
    
    for(i=1;i<=n;i++)
    {
        ipt=(int)input[i];
        INFO_Printf("ipt = %6d species = %1d linklist = [%6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d]\n",
                    ipt,species[ipt],
                    linklist[ipt*LINKLEN+0],
                    linklist[ipt*LINKLEN+1],
                    linklist[ipt*LINKLEN+2],
                    linklist[ipt*LINKLEN+3],
                    linklist[ipt*LINKLEN+4],
                    linklist[ipt*LINKLEN+5],
                    linklist[ipt*LINKLEN+6],
                    linklist[ipt*LINKLEN+7],
                    linklist[ipt*LINKLEN+8],
                    linklist[ipt*LINKLEN+9],
                    linklist[ipt*LINKLEN+10],                        
                    linklist[ipt*LINKLEN+11]);
    }
}

void YSZFrame::PrintCationList()
{
    int i, ipt;
    for(i=0;i<NCation;i++)
    {
        ipt=cationlist[i];
        INFO_Printf("i = %6d ipt = %6d species = %2d SR = [%9.2e %9.2e %9.2e]\n",
                    i,ipt,species[ipt],_SR[ipt].x,_SR[ipt].y,_SR[ipt].z);
    }
}

void YSZFrame::PrintAnionList()
{
    int i, ipt;
    for(i=0;i<NAnion;i++)
    {
        ipt=anionlist[i];
        INFO_Printf("i = %6d ipt = %6d species = %2d SR = [%9.2e %9.2e %9.2e]\n",
                    i,ipt,species[ipt],_SR[ipt].x,_SR[ipt].y,_SR[ipt].z);
    }
}

void YSZFrame::BuildEnergyBarrierList()
{
    int i, j, k, anionid, barrier_group, cations[6];
    int idir[NUMJUMPDIR] = {+1,-1,+2,-2,+3,-3}; /* +x,-x,+y,-y,+z,-z */
    FILE *fp;
    
    /* fills up Ebarrier array */
//    for(i=0;i<57;i++) Ebarrier_tab[i]=0;
    for(i=0;i<64;i++) Ebarrier_tab[i]=0;

    /* data base from ab initio calculations
     * the paper, J. Appl. Phys. 98, 103513 (2005)
     * has 42 entries.
     * Format: (e.g.)
     *   Zr Zr Zr Y Y Y -> 0*1 + 0*2 + 0*4 + 1*8 + 1*16 + 1*32
     */
    double Eb1, Eb2, Eb3;

    Eb1 = EbZZ_Car;    // original: 0.58
    Eb2 = EbZY_Car;
    Eb3 = EbYY_Car;

    
    Ebarrier_tab[0]  = Eb1;
    Ebarrier_tab[1]  = Eb1;  
    Ebarrier_tab[2]  = Eb1;
    Ebarrier_tab[4]  = Eb2;
    Ebarrier_tab[8]  = Eb2;
    Ebarrier_tab[16] = Eb1;
    Ebarrier_tab[32] = Eb1;
    Ebarrier_tab[3]  = Eb1;
    Ebarrier_tab[5]  = Eb2;
    Ebarrier_tab[9]  = Eb2;
    Ebarrier_tab[17] = Eb1;
    Ebarrier_tab[33] = Eb1;
    Ebarrier_tab[6]  = Eb2;
    Ebarrier_tab[10] = Eb2;
    Ebarrier_tab[18] = Eb1;
    Ebarrier_tab[34] = Eb1;
    Ebarrier_tab[12] = Eb3;
    Ebarrier_tab[20] = Eb2;
    Ebarrier_tab[36] = Eb2;
    Ebarrier_tab[24] = Eb2;
    Ebarrier_tab[40] = Eb2;
    Ebarrier_tab[48] = Eb1;
    Ebarrier_tab[7]  = Eb2;
    Ebarrier_tab[11] = Eb2;
    Ebarrier_tab[19] = Eb1;
    Ebarrier_tab[35] = Eb1;
    Ebarrier_tab[13] = Eb3;
    Ebarrier_tab[21] = Eb2;
    Ebarrier_tab[37] = Eb2;
    Ebarrier_tab[25] = Eb2;
    Ebarrier_tab[41] = Eb2;
    Ebarrier_tab[49] = Eb1;
    Ebarrier_tab[14] = Eb3;
    Ebarrier_tab[22] = Eb2;
    Ebarrier_tab[38] = Eb2;
    Ebarrier_tab[26] = Eb2;
    Ebarrier_tab[42] = Eb2;
    Ebarrier_tab[50] = Eb1;
    Ebarrier_tab[28] = Eb3;
    Ebarrier_tab[44] = Eb3;
    Ebarrier_tab[52] = Eb2;
    Ebarrier_tab[56] = Eb2;
    Ebarrier_tab[15] = Eb3;
    Ebarrier_tab[23] = Eb2;
    Ebarrier_tab[39] = Eb2;
    Ebarrier_tab[27] = Eb2;
    Ebarrier_tab[43] = Eb2;
    Ebarrier_tab[51] = Eb1;
    Ebarrier_tab[29] = Eb3;
    Ebarrier_tab[45] = Eb3;
    Ebarrier_tab[53] = Eb2;
    Ebarrier_tab[57] = Eb2;
    Ebarrier_tab[30] = Eb3;
    Ebarrier_tab[46] = Eb3;
    Ebarrier_tab[54] = Eb2;
    Ebarrier_tab[58] = Eb2;
    Ebarrier_tab[60] = Eb3;
    Ebarrier_tab[31] = Eb3;
    Ebarrier_tab[47] = Eb3;
    Ebarrier_tab[55] = Eb2;
    Ebarrier_tab[59] = Eb2;
    Ebarrier_tab[61] = Eb3;
    Ebarrier_tab[62] = Eb3;
    Ebarrier_tab[63] = Eb3;

    /* construct the entire Ebarrier array for all possible jumps */
    /* Ebarrier has size of _NP by NUMJUMPDIR, because anionlist */
    /* is not stored after generating Ebarrier.dat*/
    for(i=0;i<NAnion;i++) {
        anionid = anionlist[i];
        for(j=0;j<NUMJUMPDIR;j++)
        {
            /* anionlist[i] is the position of atom i in the array [0:NAnion-1] */
            MigCellCationSites_sorted(anionid,idir[j],cations);
            barrier_group = 0;
            for(k=0;k<6;k++)
            {
                if(species[cations[k]]==0) /* Zr */
                {
                    barrier_group += 0;
                }
                else /* Y */
                {
//                    barrier_group += 1<<k;
                    barrier_group += (int) pow(2,k);
                }
            }
            Ebarrier[anionid*NUMJUMPDIR+j] = Ebarrier_tab[barrier_group];
            if(Ebarrier[anionid*NUMJUMPDIR+j]<=0)
            {
                ERROR("Ebarrier = "<<Ebarrier[i*NUMJUMPDIR+j]
                      <<" site("<<i<<") jumpdir("<<idir<<")");
            }
        }
    }

    /* save Ebarrier data into a file */
    fp=fopen("Ebarrier.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildEnergyBarrierList: file open error");
        return;
    }
    fprintf(fp,"%d %d %d %d %d %d\n",NCellX,NCellY,NCellZ,NAnion,NUMJUMPDIR,NCation);
//    INFO_Printf("%d %d %d %d %d %d\n",NCellX,NCellY,NCellZ,NAnion,NUMJUMPDIR,NCation);
    for(i=0;i<_NP;i++)
        fprintf(fp,"%13.4e %13.4e %13.4e %13.4e %13.4e %13.4e\n",
                Ebarrier[i*NUMJUMPDIR+0],Ebarrier[i*NUMJUMPDIR+1],
                Ebarrier[i*NUMJUMPDIR+2],Ebarrier[i*NUMJUMPDIR+3],
                Ebarrier[i*NUMJUMPDIR+4],Ebarrier[i*NUMJUMPDIR+5]);

    fclose(fp);
    INFO_Printf("written into Ebarrier.dat\n");
    
    /* Future KMC simulations can read this file and directly work on
     *  a simple cubic lattice, without the need to deal with the YSZ
     *  crystal structure
     */
}

int YSZFrame::ReadEnergyBarrierList()
{
    FILE *fp;
    char string[500];
    int i, numjumpdir;
    
    LFile::SubHomeDir(ebarrierfile,ebarrierfile);
    
    fp=fopen(ebarrierfile,"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadEnergyBarrierList file ("<<ebarrierfile<<") not found!");
        return -1;
    }

    fgets(string,500,fp);
    sscanf(string,"%d %d %d %d %d %d\n",&NCellX,&NCellY,&NCellZ,&NAnion,&numjumpdir,&NCation);
    INFO_Printf("%d %d %d %d %d %d\n",NCellX,NCellY,NCellZ,NAnion,numjumpdir,NCation);
    
    if(numjumpdir!=NUMJUMPDIR)
    {
        ERROR("YSZFrame::ReadEnergyBarrierList numjumpdir = "<<numjumpdir
              <<" ! = NUMJUMPDIR ("<<NUMJUMPDIR);
    }

    //AllocLinkList();
    
    for(i=0;i<_NP;i++)
    {
        fgets(string,500,fp);
        sscanf(string,"%lf %lf %lf %lf %lf %lf",
              Ebarrier+(i*NUMJUMPDIR+0),Ebarrier+(i*NUMJUMPDIR+1),
              Ebarrier+(i*NUMJUMPDIR+2),Ebarrier+(i*NUMJUMPDIR+3),
              Ebarrier+(i*NUMJUMPDIR+4),Ebarrier+(i*NUMJUMPDIR+5) );
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadEnergyBarrierList done\n");
    return 0;
}

void YSZFrame::CreateConfig_Vacancy_Only()
{    /* This function is to generate an array, vaclist
        , which stores the only vacancy IDs among all Anions*/   
    int i, j;
    j = 0;
    for (i=0;i<_NP;i++) {
        if (species[i] == 3) {
            vaclist[j] = i;
//            INFO_Printf("vaclist[%d] = %d\n",j,i);
            j++;
        }
    }

    Realloc(vaclist,int,j);
    num_vac = j;            /*num_jump should be same as nV in Add_Vacancy*/
}

void YSZFrame::MigCellCationSites_unsorted(int ianion, int idir, int cations[6])
{
    int i, j, inbr, n, candidate;

    /* To Do:
     *       need to rotate the cations to standard form 
     */

    for(i=0;i<4;i++)
        cations[i] = linklist[ianion*LINKLEN+1+i];
    cations[5] = cations[6] = -1;

    switch(idir) {
    case(+1) : inbr = linklist[ianion*LINKLEN+5]; break;
    case(-1) : inbr = linklist[ianion*LINKLEN+6]; break;
    case(+2) : inbr = linklist[ianion*LINKLEN+7]; break;
    case(-2) : inbr = linklist[ianion*LINKLEN+8]; break;
    case(+3) : inbr = linklist[ianion*LINKLEN+9]; break;
    case(-3) : inbr = linklist[ianion*LINKLEN+10]; break;
    default: inbr = -1;
        ERROR("MigCellCationSites_unsorted: unrecognized idir="<<idir);
        break;
    }

    n=4;
    for(i=0;i<4;i++)
    {
        candidate = linklist[inbr*LINKLEN+1+i];
        /* find whether candidate already appears in cations list */
        for(j=0;j<n;j++)
            if(cations[j]==candidate) break;
        if(j==n)
        {
            cations[n]=candidate;
            n++;
        }
    }
    if(n!=6) ERROR("MigCellCationSites_unsorted: failed to find 6 cations, instead n="<<n);

/*
    INFO_Printf("MigCellCationSites_sorted:   cations key = %d\n",
                cations[0]+cations[1]+cations[2]+cations[3]+cations[4]+cations[5]);
    INFO_Printf("MigCellCationSites_unsorted: cations = [%d %d %d %d %d %d]\n", cations[0],cations[1],
                 cations[2],cations[3],cations[4],cations[5]);
    INFO_Printf(" species = [%d %d %d %d %d]\n",species[cations[0]],species[cations[1]],
                 species[cations[2]],species[cations[3]],species[cations[4]],species[cations[5]]);
*/
}

void YSZFrame::MigCellCationSites_sorted(int ianion, int idir, int cations[6])
{ /* Cations rotated to standard form */

  int subcelltype, inbr, *index1, *index2;

  switch(idir) {
    case(+1) : inbr = linklist[ianion*LINKLEN+5]; break;
    case(-1) : inbr = linklist[ianion*LINKLEN+6]; break;
    case(+2) : inbr = linklist[ianion*LINKLEN+7]; break;
    case(-2) : inbr = linklist[ianion*LINKLEN+8]; break;
    case(+3) : inbr = linklist[ianion*LINKLEN+9]; break;
    case(-3) : inbr = linklist[ianion*LINKLEN+10]; break;
    default: inbr = -1;
        ERROR("MigCellCationSites_sorted: unrecognized idir="<<idir);
        break;
  }

  index1 = linklist + ianion*LINKLEN;
  index2 = linklist + inbr*LINKLEN;
  subcelltype = linklist[ianion*LINKLEN+11];
  
  if (subcelltype == 1) /* subcell type 1 */
    switch (idir) {
        case +1:
          cations[0] = *(index1+1); cations[1] = *(index1+2); cations[2] = *(index1+3); cations[3] = *(index1+4);
          cations[4] = *(index2+3); cations[5] = *(index2+4); break;
        case -1:
          cations[0] = *(index1+3); cations[1] = *(index1+4); cations[2] = *(index1+1); cations[3] = *(index1+2);
          cations[4] = *(index2+1); cations[5] = *(index2+2); break;
        case +2:
          cations[0] = *(index1+1); cations[1] = *(index1+3); cations[2] = *(index1+2); cations[3] = *(index1+4);
          cations[4] = *(index2+2); cations[5] = *(index2+4); break;
        case -2:
          cations[0] = *(index1+2); cations[1] = *(index1+4); cations[2] = *(index1+1); cations[3] = *(index1+3);
          cations[4] = *(index2+1); cations[5] = *(index2+3); break;
        case +3:
          cations[0] = *(index1+1); cations[1] = *(index1+4); cations[2] = *(index1+2); cations[3] = *(index1+3);
          cations[4] = *(index2+1); cations[5] = *(index2+4); break;
        case -3:
          cations[0] = *(index1+2); cations[1] = *(index1+3); cations[2] = *(index1+1); cations[3] = *(index1+4);
          cations[4] = *(index2+2); cations[5] = *(index2+3); break;
        default: FATAL("MigCellCationSites_sorted: unrecognized idir="<<idir); break;
    }
  else if (subcelltype == 2) /* subcell type 2 */
    switch (idir) {
        case +1:
          cations[0] = *(index1+1); cations[1] = *(index1+2); cations[2] = *(index1+3); cations[3] = *(index1+4);
          cations[4] = *(index2+3); cations[5] = *(index2+4); break;
        case -1:
          cations[0] = *(index1+3); cations[1] = *(index1+4); cations[2] = *(index1+1); cations[3] = *(index1+2);
          cations[4] = *(index2+1); cations[5] = *(index2+2); break;
        case +2:
          cations[0] = *(index1+1); cations[1] = *(index1+3); cations[2] = *(index1+2); cations[3] = *(index1+4);
          cations[4] = *(index2+2); cations[5] = *(index2+4); break;
        case -2:
          cations[0] = *(index1+2); cations[1] = *(index1+4); cations[2] = *(index1+1); cations[3] = *(index1+3);
          cations[4] = *(index2+1); cations[5] = *(index2+3); break;
        case +3:
          cations[0] = *(index1+2); cations[1] = *(index1+3); cations[2] = *(index1+1); cations[3] = *(index1+4);
          cations[4] = *(index2+2); cations[5] = *(index2+3); break;
        case -3:
          cations[0] = *(index1+1); cations[1] = *(index1+4); cations[2] = *(index1+2); cations[3] = *(index1+3);
          cations[4] = *(index2+1); cations[5] = *(index2+4); break;
    default: FATAL("MigCellCationSites_sorted: unrecognized idir="<<idir); break;
    }
  
/*
    INFO_Printf("MigCellCationSites_sorted:   cations key = %d  subcell=%d\n",
                cations[0]+cations[1]+cations[2]+cations[3]+cations[4]+cations[5],
                subcelltype);
    INFO_Printf("MigCellCationSites_sorted:   cations = [%d %d %d %d %d %d]  subcell=%d\n", cations[0],cations[1],
                 cations[2],cations[3],cations[4],cations[5],linklist[ianion*LINKLEN+11]);
    INFO_Printf(" species = [%d %d %d %d %d]\n",species[cations[0]],species[cations[1]],
                 species[cations[2]],species[cations[3]],species[cations[4]],species[cations[5]]);
*/
}


void YSZFrame::Add_Yttrium() {
    /* There cannot be more than three Yttrium atoms in one migration unit cell
     * MC_NY array stores the number of Y atoms in each migration unit cell
     */
    int i, j, k, n, nY, nallow, nsum, *Yallow;
    int ianion, cations[6];
    int idir[NUMJUMPDIR] = {+1,-1,+2,-2,+3,-3}; /* +x,-x,+y,-y,+z,-z */

    for (i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

/*  This part is to add Yttrium at specific position for comparison */
/*
    int Yx,Yy,Yz,catid,catx,caty,catz;
    nY = 10;
    int Ypos[30] = {0,0,0,
                    0,2,2,
                    4,4,4,
                    4,8,0,
                    10,2,8,
                    12,6,6,
                    12,12,12,
                    14,18,20,
                    20,4,8,
                    24,20,24};
    for (i=0;i<10;i++) {
        Yx = Ypos[3*i+0];
        Yy = Ypos[3*i+1];
        Yz = Ypos[3*i+2];
        for (j=0;j<NCation;j++) {
            catid = *(cationlist+j);
            catx = _RI[catid*3+0];
            caty = _RI[catid*3+1];
            catz = _RI[catid*3+2];
            if (catx == Yx && caty == Yy && catz == Yz) {
                species[catid] = 2;
                INFO_Printf("Added Y at (%d,%d,%d)\n",catx,caty,catz);
                break;
            }
        }
    }
*/   
    
#if 1 
    nY = (int) input[0];
//    INFO_Printf("Add_Yttrium: nY=%d  NCation = %d\n", nY, NCation);
    if(nY<=0)
    {
        ERROR("Add_Yttrium: nY="<<nY<<" NCation="<<NCation);
        return;
    }

    /* species:
     *         Zr: 0
     *         O : 1
     *         Y : 2
     *         Vo: 3
     */
    Yallow = 0;
    Realloc(Yallow,int,_NP);

    memset(Yallow,0,sizeof(int)*_NP);

    nallow = 0;
    for(i=0;i<_NP;i++) /* establish allowable sites for Y */
    {
        if(species[i]==0) /* only Zr atoms can change to Y */
        {
            Yallow[i]=1;
            nallow++;
        }
    }

    for(i=0;i<nY;i++)
    {
        if(nallow<=0) /* nallow is the number of allowable sites to add Y */
        {
            INFO("  "<<i<<" instead of "<<nY<<" Yttrium atoms added");
            INFO("NP = "<<_NP<<" NCation = "<<NCation<<" NAnion = "<<NAnion);
            ERROR("Add_Yttrium: nallow="<<nallow);
            break;
        }
        n = (int) floor(drand48()*nallow); /* n is a random number from 0 to nallow-1 */
        k = 0;
        for(j=0;j<_NP;j++)
        {
            if (!Yallow[j]) continue; /* skip sites that are not allowed */
            if (n==k) break;
            k++;
        }
        /* insert Y at site j */
        Yallow[j]=0;
        species[j]=2; /* 0-Zr, 1-O, 2-Y, 3-Vo */
        linklist[j*LINKLEN+0] = 0;
        INFO_Printf("Add_Yttrium: select atom %6d(%3d,%3d,%3d) (n=%6d nallow=%6d)\n",j,_RI[3*j],_RI[3*j+1],_RI[3*j+2],n,nallow);

        /* decide how many sites remain allowed to be changed to Y
         * a site is not allowed if it is already occupied by Y
         * or if changing it to Y would make some migration unit cell have 3 Y atoms
         */
        for(n=0;n<8;n++)
        {/* loop through all 8 Anion sites connected to the chosen site */
            ianion = linklist[j*LINKLEN+1+n];
            for(k=0;k<6;k++)
            {/* loop through all 6 possible jump directions */

                /* compute total number of Y atoms in the nearest neighbor of
                 * ianion and its anion neighbor
                 */
                /* MigCellCationSites_unsorted(ianion,idir[k],cations); */
                MigCellCationSites_sorted  (ianion,idir[k],cations);

                nsum = (species[cations[0]]==2) + (species[cations[1]]==2)
                     + (species[cations[2]]==2) + (species[cations[3]]==2)
                     + (species[cations[4]]==2) + (species[cations[5]]==2) ;
                
                /* if there are already 3 Y atoms, then all the Cation sites
                 * neighboring ianion and ianionnbr are not allowed to change to Y
                 */
                if(nsum==3)
                {
                    Yallow[cations[0]] = Yallow[cations[1]] = Yallow[cations[2]]
                 =  Yallow[cations[3]] = Yallow[cations[4]] = Yallow[cations[5]] = 0;
                }
                else if(nsum>3)
                {
                    ERROR("Add_Yttrium:  more than 3 Y atoms in migration cell ");
                }
            }
        }
    
        /* determine how many sites are still allowabled */
        nallow = 0;
        for(j=0;j<_NP;j++)
            nallow += Yallow[j];
    }
    Free(Yallow);
#endif
}

void YSZFrame::Add_Yttrium_Manually() {
    /* There cannot be more than three Yttrium atoms in one migration unit cell
     * MC_NY array stores the number of Y atoms in each migration unit cell
     */
    int i, j, k, n, nY, nsum;
    int ianion, cations[6];
    int idir[NUMJUMPDIR] = {+1,-1,+2,-2,+3,-3}; /* +x,-x,+y,-y,+z,-z */
    
    for (i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }
    INFO_Printf("NCation=%d\n",NCation);
/*  This part is to add Yttrium at specific position for comparison */
//
//    int Yx,Yy,Yz,catid,catx,caty,catz;    
//    nY = 1000;
//    int Ypos[3000];
//    n = 0;
    // uniform distribution
//    for (i=0;i<10;i++)
//        for (j=0;j<10;j++)
//            for (k=0;k<10;k++) {
//                Ypos[3*n+0] = i*4;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }
    // from this line, Y block in z direction
//    for (i=0;i<10;i++)
//        for (j=0;j<10;j++)
//            for (k=0;k<2;k++) {
//                Ypos[3*n+0] = i*4;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//                Ypos[3*n+0] = i*4;
//                Ypos[3*n+1] = j*4+2;
//                Ypos[3*n+2] = k*4+2;
//                n++;
//                Ypos[3*n+0] = i*4+2;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4+2;
//                n++;
//                Ypos[3*n+0] = i*4+2;
//                Ypos[3*n+1] = j*4+2;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }          
//    for (i=0;i<10;i++)
//        for (j=0;j<10;j++)
//            for (k=2;k<3;k++) {
//                Ypos[3*n+0] = i*4;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//                Ypos[3*n+0] = i*4+2;
//                Ypos[3*n+1] = j*4+2;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }
//    //until this line
    // from this line, pillars in x direction
//    int nlayerid[6]={0,1,2,5,6,7};
//    for (i=0;i<6;i++)
//        for (j=0;j<6;j++)
//            for (k=0;k<6;k++) {
//                Ypos[3*n+0] = nlayerid[i]*4;
//                Ypos[3*n+1] = nlayerid[j]*4;
//                Ypos[3*n+2] = nlayerid[k]*4;
//                n++;
//                Ypos[3*n+0] = nlayerid[i]*4;
//                Ypos[3*n+1] = nlayerid[j]*4+2;
//                Ypos[3*n+2] = nlayerid[k]*4+2;
//                n++;
//                Ypos[3*n+0] = nlayerid[i]*4+2;
//                Ypos[3*n+1] = nlayerid[j]*4;
//                Ypos[3*n+2] = nlayerid[k]*4+2;
//                n++;
//                Ypos[3*n+0] = nlayerid[i]*4+2;
//                Ypos[3*n+1] = nlayerid[j]*4+2;
//                Ypos[3*n+2] = nlayerid[k]*4;
//                n++;
//            }
//    int nlayerid1[4]={0,1,5,6};
//    int nlayerid2[2]={3,8};
//    for (i=0;i<2;i++)
//        for (j=0;j<4;j++)
//            for (k=0;k<4;k++) {
//                Ypos[3*n+0] = nlayerid2[i]*4;
//                Ypos[3*n+1] = nlayerid1[j]*4;
//                Ypos[3*n+2] = nlayerid1[k]*4;
//                n++;
//                Ypos[3*n+0] = nlayerid2[i]*4;
//                Ypos[3*n+1] = nlayerid1[j]*4+2;
//                Ypos[3*n+2] = nlayerid1[k]*4+2;
//                n++;
//                Ypos[3*n+0] = nlayerid2[i]*4+2;
//                Ypos[3*n+1] = nlayerid1[j]*4;
//                Ypos[3*n+2] = nlayerid1[k]*4+2;
//                n++;
//                Ypos[3*n+0] = nlayerid2[i]*4+2;
//                Ypos[3*n+1] = nlayerid1[j]*4+2;
//                Ypos[3*n+2] = nlayerid1[k]*4;
//                n++;
//            }
//    for (i=0;i<2;i++)
//        for (j=2;j<3;j++)
//            for (k=2;k<3;k++) {
//                Ypos[3*n+0] = nlayerid2[i]*4;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//                Ypos[3*n+0] = nlayerid2[i]*4;
//                Ypos[3*n+1] = j*4+2;
//                Ypos[3*n+2] = k*4+2;
//                n++;
//                Ypos[3*n+0] = nlayerid2[i]*4+2;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4+2;
//                n++;
//                Ypos[3*n+0] = nlayerid2[i]*4+2;
//                Ypos[3*n+1] = j*4+2;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }   
//    //until this line
    // uniform cluster distribution
//    for (i=0;i<5;i++)
//        for (j=0;j<3;j++)
//            for (k=0;k<1;k++) {
//                Ypos[3*n+0] = i*4; //(0,0,0) case 
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }
//    for (i=0;i<5;i++)
//        for (j=0;j<2;j++)
//            for (k=1;k<2;k++) {
//                Ypos[3*n+0] = i*4; //(0,0,0) case 
//                Ypos[3*n+1] = j*8;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }
//    for (i=0;i<5;i++)
//        for (j=0;j<3;j++)
//            for (k=2;k<3;k++) {
//                Ypos[3*n+0] = i*4; //(0,0,0) case 
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }
//    for (i=0;i<5;i++)
//        for (j=0;j<2;j++)
//            for (k=0;k<2;k++) {
//                Ypos[3*n+0] = i*4+2; //(2,2,0) case
//                Ypos[3*n+1] = j*4+2;
//                Ypos[3*n+2] = k*8;
//                n++;
//            }          

//    for (i=0;i<10;i++)
//        for (j=0;j<5;j++)
//            for (k=0;k<10;k++) {
//                Ypos[6*n+0] = i*4; //(0,0,0) case 
//                Ypos[6*n+1] = j*8;
//                Ypos[6*n+2] = k*4;
//                Ypos[6*n+3] = i*4+2; //(2,2,0) case
//                Ypos[6*n+4] = j*8+2;
//                Ypos[6*n+5] = k*4;
//                n++;
//            }          
//    for (i=0;i<10;i++)
//        for (j=0;j<10;j++)
//            for (k=0;k<10;k++) {
//                Ypos[3*n+0] = i*4;
//                Ypos[3*n+1] = j*4;
//                Ypos[3*n+2] = k*4;
//                n++;
//            }
//    
//    for (i=0;i<nY;i++) {
//        Yx = Ypos[3*i+0];
//        Yy = Ypos[3*i+1];
//        Yz = Ypos[3*i+2];
//        for (j=0;j<NCation;j++) {
//            catid = *(cationlist+j);
//            catx = _RI[catid*3+0];
//            caty = _RI[catid*3+1];
//            catz = _RI[catid*3+2];
//            if (catx == Yx && caty == Yy && catz == Yz) {
//                species[catid] = 2; /* 0-Zr, 1-O, 2-Y, 3-Vo */
//                linklist[catid*LINKLEN+0] = 0;
//                INFO_Printf("Added %d-th Y at (%d,%d,%d)\n",i,catx,caty,catz);
//                // from this line, we don't need any more. It's about 3Y rule
//                for(n=0;n<8;n++)
//                {/* loop through all 8 Anion sites connected to the chosen site */
//                    ianion = linklist[catid*LINKLEN+1+n];
//                    for(k=0;k<6;k++)
//                    {/* loop through all 6 possible jump directions */
//                        
//                        /* compute total number of Y atoms in the nearest neighbor of
//                         * ianion and its anion neighbor
//                         */
//                        /* MigCellCationSites_unsorted(ianion,idir[k],cations); */
//                        MigCellCationSites_sorted  (ianion,idir[k],cations);
//                        
//                        nsum = (species[cations[0]]==2) + (species[cations[1]]==2)
//                            + (species[cations[2]]==2) + (species[cations[3]]==2)
//                            + (species[cations[4]]==2) + (species[cations[5]]==2) ;
//                        
//                        /* if there are already 3 Y atoms, then all the Cation sites
//                         * neighboring ianion and ianionnbr are not allowed to change to Y
//                         */
//                        if (nsum>3)
//                            ERROR("Add_Yttrium:  more than 3 Y atoms in migration cell ");
//                    }
//                }
//                // until this line
//                break;
//            }
//            
//        }
//    }
//
    // for cluster distribution without 3Y rule   
    int Yx,Yy,Yz,catid,catx,caty,catz;    
    nY = (int) input[0];
    int Ypos[3*NCation];
    n = 0;
    double rNN0 = 0.0;
    double rNN1 = 0.05;
    double distNN;
    while (n<NCation) {
        for (i=0;i<NCation;i++) {
            catid = *(cationlist+i);
            catx = _RI[catid*3+0];
            caty = _RI[catid*3+1];
            catz = _RI[catid*3+2];
            distNN = sqrt((catx-NCellX*0.5)*(catx-NCellX*0.5)*1.0+(caty-NCellY*0.5)*(caty-NCellY*0.5)*1.0+(catz-NCellZ*0.5)*(catz-NCellZ*0.5)*1.0);
            if (distNN >= rNN0 && distNN < rNN1) {
                Ypos[3*n+0] = catx;
                Ypos[3*n+1] = caty;
                Ypos[3*n+2] = catz;
                n++;
                INFO_Printf("%d-th NN.id=%d, info, distNN=%e,rNN0=%e,rNN1=%e\n",n,catid,distNN,rNN0,rNN1);
                if (n == NCation)
                    break;
            }            
        }
        rNN0 = rNN0 + 0.05;
        rNN1 = rNN1 +0.05;
    }
    INFO_Printf("n=%d, NCation=%d\n",n,NCation);
    int mY = 0;
//    int errctr;
    for (i=0;i<NCation;i++) {
        Yx = Ypos[3*i+0];
        Yy = Ypos[3*i+1];
        Yz = Ypos[3*i+2];
        for (j=0;j<NCation;j++) {
            catid = *(cationlist+j);
            catx = _RI[catid*3+0];
            caty = _RI[catid*3+1];
            catz = _RI[catid*3+2];
            if (catx == Yx && caty == Yy && catz == Yz) {
                species[catid] = 2; /* 0-Zr, 1-O, 2-Y, 3-Vo */
                linklist[catid*LINKLEN+0] = 0;
//                errctr = 0;
//                for(n=0;n<8;n++)
//                {/* loop through all 8 Anion sites connected to the chosen site */
//                    ianion = linklist[catid*LINKLEN+1+n];
//                    for(k=0;k<6;k++)
//                    {/* loop through all 6 possible jump directions */
//                        
//                        /* compute total number of Y atoms in the nearest neighbor of
//                         * ianion and its anion neighbor
//                         */
//                        /* MigCellCationSites_unsorted(ianion,idir[k],cations); */
//                        MigCellCationSites_sorted  (ianion,idir[k],cations);
//                        
//                        nsum = (species[cations[0]]==2) + (species[cations[1]]==2)
//                            + (species[cations[2]]==2) + (species[cations[3]]==2)
//                            + (species[cations[4]]==2) + (species[cations[5]]==2) ;
//                        
//                        /* if there are already 3 Y atoms, then all the Cation sites
//                         * neighboring ianion and ianionnbr are not allowed to change to Y
//                         */
//                        if (nsum>3) {
//                            ERROR("Add_Yttrium:  more than 3 Y atoms in migration cell ");
//                            errctr++;
//                        }
//                    }
//                }
//                if (errctr > 0)
//                    species[catid] = 0;
//                else {
//                    mY++;
//                    INFO_Printf("Added %d-th Y at (%d,%d,%d)\n",mY,catx,caty,catz);
//                }
                mY++;
                INFO_Printf("Added %d-th Y at (%d,%d,%d)\n",mY,catx,caty,catz);
                INFO_Printf("%d,%d,%d,%d;\n",mY,catx,caty,catz);
                break;
            }
        }
        if (mY == nY)
            break;
    }
}

void YSZFrame::Add_Yttrium_Random() {
    /* There cannot be more than three Yttrium atoms in one migration unit cell
     * MC_NY array stores the number of Y atoms in each migration unit cell
     */
    int i, j, k, n, nY, mY, catid, catx, caty, catz;
//    int ianion, cations[6];
//    int idir[NUMJUMPDIR] = {+1,-1,+2,-2,+3,-3}; /* +x,-x,+y,-y,+z,-z */

    for (i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    nY = (int) input[0];
//    INFO_Printf("Add_Yttrium: nY=%d  NCation = %d\n", nY, NCation);
    if(nY<=0)
    {
        ERROR("Add_Yttrium: nY="<<nY<<" NCation="<<NCation);
        return;
    }

    /* species:
     *         Zr: 0
     *         O : 1
     *         Y : 2
     *         Vo: 3
     */
    mY = 0;
    while (mY<nY) {
        j = (int) floor(drand48()*NCation);
        catid = *(cationlist+j);
        if (species[catid] == 0) {
            species[catid]=2; /* 0-Zr, 1-O, 2-Y, 3-Vo */
            linklist[catid*LINKLEN+0] = 0;
            catx = _RI[catid*3+0];
            caty = _RI[catid*3+1];
            catz = _RI[catid*3+2];
	    INFO_Printf("Added Yttrium on Catlist-%d ,_NPlist-%d at (%d,%d,%d)\n",j,catid,catx,caty,catz);
            mY++;
//            INFO_Printf("%d,%d,%d,%d;\n",mY+1000-nY,catx,caty,catz);
        }
    }
}

void YSZFrame::Add_Yttrium_On_Surface() {
    int i, j, k, n, nY, nallow;

    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    nY = (int) input[0];
    INFO_Printf("Add_Yttrium: nY=%d  NCation = %d\n", nY, NCation);
    if(nY<=0)
    {
        ERROR("Add_Yttrium: nY="<<nY<<" NCation="<<NCation);
        return;
    }
    
    int ZLimMin  = (int) input[1];
    int ZLimMax = (int) input[2];

    /* species:
     *         Zr: 0
     *         O : 1
     *         Y : 2
     *         Vo: 3
     */
    
    nallow = 0;
    for(j=0;j<_NP;j++)   /* only Zr sites bet. ZLims allowed to be vacant */
        if (_RI[3*j+2] > ZLimMin && _RI[3*j+2] < ZLimMax)
            nallow += (species[j]==0);  /* 0-Zr */    
    
    for(i=0;i<nY;i++)
    {
        if(nallow<=0) /* nallow is the number of allowable sites to add Y */
        {
            INFO("  "<<i<<" instead of "<<nY<<" yttriums added");
            INFO("NP = "<<_NP<<" NCation = "<<NCation<<" NAnion = "<<NAnion);
            ERROR("Add_Yttrium: nallow="<<nallow);
            break;
        }
        
        n = (int) floor(drand48()*nallow);        
        k = 0;
        for(j=0;j<_NP;j++)
        {
            if (species[j]!=0) continue;
            if (_RI[3*j+2] <= ZLimMin || _RI[3*j+2] >= ZLimMax) continue;
            if (n==k)
                break;
            k++;
        }
        /* insert Vo at j */
        species[j]=2; /* 0-Zr, 1-O, 2-Y, 3-Vo */
        linklist[j*LINKLEN+0] = 0;
        INFO_Printf("Add_Yttrium: select atom %d for %d-th Y (%d,%d,%d)\n",j,i+1,_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);

        nallow--;
    }    
}

void YSZFrame::Add_Yttrium_On_Surface_X() {
    int i, j, k, n, nY, nallow;

    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    nY = (int) input[0];
    INFO_Printf("Add_Yttrium: nY=%d  NCation = %d\n", nY, NCation);
    if(nY<=0)
    {
        ERROR("Add_Yttrium: nY="<<nY<<" NCation="<<NCation);
        return;
    }

    int XLimMin  = (int) input[1];
    int XLimMax = (int) input[2];

    /* species:
     *         Zr: 0
     *         O : 1
     *         Y : 2
     *         Vo: 3
     */

    nallow = 0;
    for(j=0;j<_NP;j++)   /* only Zr sites bet. ZLims allowed to be vacant */
        if (_RI[3*j+0] > XLimMin && _RI[3*j+0] < XLimMax)
            nallow += (species[j]==0);  /* 0-Zr */

    for(i=0;i<nY;i++)
    {
        if(nallow<=0) /* nallow is the number of allowable sites to add Y */
        {
            INFO("  "<<i<<" instead of "<<nY<<" yttriums added");
            INFO("NP = "<<_NP<<" NCation = "<<NCation<<" NAnion = "<<NAnion);
            ERROR("Add_Yttrium: nallow="<<nallow);
            break;
        }

        n = (int) floor(drand48()*nallow);
        k = 0;
        for(j=0;j<_NP;j++)
        {
            if (species[j]!=0) continue;
            if (_RI[3*j+0] <= XLimMin || _RI[3*j+0] >= XLimMax) continue;
            if (n==k)
                break;
            k++;
        }
        /* insert Vo at j */
        species[j]=2; /* 0-Zr, 1-O, 2-Y, 3-Vo */
        linklist[j*LINKLEN+0] = 0;
        INFO_Printf("Add_Yttrium: select atom %d for %d-th Y (%d,%d,%d)\n",j,i+1,_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);

        nallow--;
    }
}

void YSZFrame::Add_Vacancy() {
    int i, j, k, n, nV, nallow;

    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    nV = (int) input[0];
    INFO_Printf("Add_Vacancy: nV=%d  NAnion = %d\n", nV, NAnion);
    if(nV<=0)
    {
        ERROR("Add_Vacancy: nV="<<nV<<" NAnion="<<NAnion);
        return;
    }

    nallow = 0;
    for(j=0;j<_NP;j++)   /* only O sites can be turned to Vacancy */
        nallow += (species[j]==1);  /* 1-O */

    for(i=0;i<nV;i++)
    {
        if(nallow<=0) /* nallow is the number of allowable sites to add Y */
        {
            INFO("  "<<i<<" instead of "<<nV<<" vacancies added");
            INFO("NP = "<<_NP<<" NCation = "<<NCation<<" NAnion = "<<NAnion);
            ERROR("Add_Vacancy: nallow="<<nallow);
            break;
        }
        
        n = (int) floor(drand48()*nallow);
        k = 0;
        for(j=0;j<_NP;j++)
        {
            if (species[j]!=1) continue;
            if (n==k) break;
            k++;
        }
        /* insert Vo at j */
        species[j]=3; /* 0-Zr, 1-O, 2-Y, 3-Vo */
        linklist[j*LINKLEN+0] = 1;
        INFO_Printf("Add_Vacancy: select atom %6d(%d,%d,%d) (n=%6d nallow=%6d)\n",j,_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],n,nallow);

        nallow--;
    }
}

void YSZFrame::Add_Vacancy_On_Surface() {
    int i, j, k, n, nV, nallow;
    
    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    nV = (int) input[0];
    INFO_Printf("Add_Vacancy: nV=%d  NAnion = %d\n", nV, NAnion);
    if(nV<=0)
    {
        ERROR("Add_Vacancy: nV="<<nV<<" NAnion="<<NAnion);
        return;
    }

    int ZLimMin  = (int) input[1];
    int ZLimMax = (int) input[2];
    
    nallow = 0;
    for(j=0;j<_NP;j++)   /* only O sites bet. ZLims allowed to be vacant */
        if (_RI[3*j+2] > ZLimMin && _RI[3*j+2] < ZLimMax)
            nallow += (species[j]==1);  /* 1-O */    
    
    for(i=0;i<nV;i++)
    {
        if(nallow<=0) /* nallow is the number of allowable sites to add V */
        {
            INFO("  "<<i<<" instead of "<<nV<<" vacancies added");
            INFO("NP = "<<_NP<<" NCation = "<<NCation<<" NAnion = "<<NAnion);
            ERROR("Add_Vacancy: nallow="<<nallow);
            break;
        }
        
        n = (int) floor(drand48()*nallow);        
        k = 0;
        for(j=0;j<_NP;j++)
        {
            if (species[j]!=1) continue;
            if (_RI[3*j+2] <= ZLimMin || _RI[3*j+2] >= ZLimMax) continue;
            if (n==k)
                break;
            k++;
        }
        /* insert Vo at j */
        species[j]=3; /* 0-Zr, 1-O, 2-Y, 3-Vo */
        linklist[j*LINKLEN+0] = 1;
        INFO_Printf("Add_Vacancy: select atom %6d for %d-th V (%d,%d,%d)\n",j,i+1,_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);

        nallow--;
    }
}

void YSZFrame::BuildCoulombEList()
{
#if 0
    int i; double V0;
    FILE *fp;

    INFO_Printf("calling CoulombPotential\n");
    CoulombPotential(NCellX,NCellY,NCellZ,CoulombE);
//    Ewald2DPotential(NCellX,NCellY,NCellZ,CoulombE);
    INFO_Printf("finished CoulombPotential\n");

    V0 = (2*2) * 1.602e-19 / ( 4.0 * M_PI * epsilon_0 * epsilon_r * latticeconst[0]*0.25e-10 );
//    printf("epsilon_0 = %e epsilon_r = %e, latticeconst = %e, V0=%e\n",epsilon_0,epsilon_r,latticeconst[0],V0);
    for(i=0;i<NCellX*NCellY*NCellZ;i++) {
        CoulombE[i] *= V0;
//        printf("Vc[%d]=%e\n",i,CoulombE[i]);
    }

    /* save CoulombE data into a file */
    fp=fopen("CoulombE.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildCoulombElist: file open error");
        return;
    }
    fprintf(fp,"%d %d %d %g\n",NCellX,NCellY,NCellZ,latticeconst[0]/4);
    INFO_Printf("%d %d %d %g\n",NCellX,NCellY,NCellZ,latticeconst[0]/4);
    for(i=0;i<NCellX*NCellY*NCellZ/4;i++)
        fprintf(fp,"%21.14e %21.14e %21.14e %21.14e\n",
                CoulombE[i*4+0],CoulombE[i*4+1],CoulombE[i*4+2],CoulombE[i*4+3]);
    fclose(fp);
    INFO_Printf("written into CoulombE.dat\n");
    
    /* Future KMC simulations can read this file and directly work on
     *  a simple cubic lattice, without the need to deal with the YSZ
     *  crystal structure
     */
#endif
}

void YSZFrame::BuildEwald2DList()
{
#if 0
    int i; double V0;
    FILE *fp;

    INFO_Printf("calling Ewald2D Summation for CoulombPotential\n");
    Ewald2DPotential(NCellX,NCellY,NCellZ,CoulombE);
    INFO_Printf("finished Ewald2D Summation for CoulombPotential\n");

    V0 = (2*2) * 1.602e-19 / ( 4.0 * M_PI * epsilon_0 * epsilon_r * latticeconst[0]*0.25e-10 ); 
        
    for(i=0;i<2*NCellX*NCellY*NCellZ;i++) CoulombE[i] *= V0;

    /* save CoulombE data into a file */
    fp=fopen("CoulombEwald.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildEwald2Dlist: file open error");
        return;
    }
    fprintf(fp,"%d %d %d %g\n",NCellX,NCellY,NCellZ,latticeconst[0]/4);
//    INFO_Printf("%d %d %d %g\n",NCellX,NCellY,NCellZ,latticeconst[0]/4);
    for(i=0;i<NCellX*NCellY*NCellZ/2;i++)
        fprintf(fp,"%21.14e %21.14e %21.14e %21.14e\n",
                CoulombE[i*4+0],CoulombE[i*4+1],CoulombE[i*4+2],CoulombE[i*4+3]);
    fclose(fp);
    INFO_Printf("written into CoulombEwald.dat\n");   
#endif
}

void YSZFrame::BuildCatSpeciesList()
{
//  Find the number of Y in the nearest neighbor sites
//  require "BuildLinkList" to be precedent
    int ipt, jpt;
    FILE *fp;

    fp=fopen("CatSpeciesList.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildCatSpeciesList: file open error");
        return;
    }
    for(ipt=0;ipt<_NP;ipt++) {
        if (species[ipt] == 2)
            fprintf(fp,"%d\n",species[ipt]);
        else
            fprintf(fp,"%d\n",0);
    }
    fclose(fp);
    INFO_Printf("written into CatSpeciesList.dat\n");    
}

void YSZFrame::BuildNNYList()
{
//  Find the number of Y in the nearest neighbor sites
//  require "BuildLinkList" to be precedent
    int ipt, jpt;
    FILE *fp;
    Realloc(nnYlist,int,_NP);    
    memset(nnYlist,0,sizeof(int)*_NP);
    for (ipt=0;ipt<_NP;ipt++) {
        if(species[ipt]!=1) continue; 
        for (jpt=1;jpt<5;jpt++)
            nnYlist[ipt] += (int) species[linklist[ipt*LINKLEN+jpt]]/2;
    }
    fp=fopen("NNYList.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildNNYList: file open error");
        return;
    }
    for(ipt=0;ipt<_NP;ipt++)
        fprintf(fp,"%d\n",nnYlist[ipt]);
    fclose(fp);
    INFO_Printf("written into NNYList.dat\n");    
}

void YSZFrame::BuildSNYList()
{
//  Find the number of Y in the nearest neighbor sites
//  require "BuildLinkList" to be precedent
    int ipt, jpt, j, n;
    double r, rnn, rsn, rcut1, rcut2;
    Vector3 dr, ds;
    Matrix33 hinv;
    FILE *fp;

    /* nearest neighbor */
    rnn = sqrt(3)*0.25*latticeconst[0];  rcut1 = rnn*1.05;
    rsn = sqrt(11)*0.25*latticeconst[0]; rcut2 = rsn*1.05;
//    INFO_Printf("latticeconst = %f  rnn=%f  rsn = %f\n",latticeconst[0],rnn,rsn);

    Realloc(snYlist,int,_NP);    
    memset(snYlist,0,sizeof(int)*_NP);
#if 0
    for (ipt=0;ipt<_NP;ipt++) {
        if(species[ipt]!=1) continue;
	n = 0;
        for(jpt=0;jpt<_NP;jpt++)
        {
            ds = _SR[jpt]-_SR[ipt];
            ds.subint();
            dr = _H*ds;
            r = dr.norm();
            if(r<rcut2 && r>rcut1) {
	        if (species[jpt] == 2)
		    n++;
	    }
        }
	snYlist[ipt] = n;
    }
#endif
#if 1 
    for (ipt=0;ipt<_NP;ipt++) {
        if(species[ipt]!=1) continue;
	n = 0;
        for(j=0;j<nn[ipt];j++)
        {
	    jpt = nindex[ipt][j];
            ds = _SR[jpt]-_SR[ipt];
            ds.subint();
            dr = _H*ds;
            r = dr.norm();
            if(r<rcut2 && r>rcut1) {
	        if (species[jpt] == 2)
		  n++;
	    }
        }
	snYlist[ipt] = n;
    }
#endif
    fp=fopen("SNYList.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildSNYList: file open error");
        return;
    }
    for(ipt=0;ipt<_NP;ipt++)
        fprintf(fp,"%d\n",snYlist[ipt]);
    fclose(fp);
    INFO_Printf("written into SNYList.dat\n");   
}

void YSZFrame::BuildCatPotList()
{
#if 0
    int i,j;
//    double V0;
    int index1,index2,index3;
    FILE *fp;
//    run following comments only when the BuildCoulombEList is not called before
    
//    INFO_Printf("calling CoulombPotential\n");
//    CoulombPotential(NCellX,NCellY,NCellZ,CoulombE);
//    Ewald2DPotential(NCellX,NCellY,NCellZ,CoulombE);
//    INFO_Printf("finished CoulombPotential\n");

//    V0 = (2*2) * 1.602e-19 / ( 4.0 * M_PI * epsilon_0 * epsilon_r * latticeconst[0]*0.25e-10 );
//    printf("epsilon_0 = %e epsilon_r = %e, latticeconst = %e, V0=%e\n",epsilon_0,epsilon_r,latticeconst[0],V0);
//    for(i=0;i<NCellX*NCellY*NCellZ;i++) {
//        CoulombE[i] *= V0;
//    }
    for(i=0;i<_NP;i++) {
        if (species[i] == 2)    /* 0-Zr, 1-O, 2-Y, 3-Vo */
            for (j=0;j<_NP;j++) {               
                index1 = (int) mymod(*(_RI+3*i+0) - *(_RI+3*j+0) + NCellX/2,NCellX);
                index2 = (int) mymod(*(_RI+3*i+1) - *(_RI+3*j+1) + NCellY/2,NCellY);
                index3 = (int) mymod(*(_RI+3*i+2) - *(_RI+3*j+2) + NCellZ/2,NCellZ);
                CatPot[j] += -CoulombE[NCellY*NCellZ*index1 + NCellZ*index2 + index3]*0.5;
            }
//        if (species[i] == 0)   /* 0-Zr, 1-O, 2-Y, 3-Vo */
//            for (j=0;j<_NP;j++) {               
//                index1 = (int) mymod(*(_RI+3*i+0) - *(_RI+3*j+0) + NCellX/2,NCellX);
//                index2 = (int) mymod(*(_RI+3*i+1) - *(_RI+3*j+1) + NCellY/2,NCellY);
//                index3 = (int) mymod(*(_RI+3*i+2) - *(_RI+3*j+2) + NCellZ/2,NCellZ);
//                CatPot[j] += -CoulombE[NCellY*NCellZ*index1 + NCellZ*index2 + index3]*2.0;
//            }
    }
    /* save CoulombE data into a file */
    fp=fopen("CatPotList.dat","w");
    if(fp==NULL)
    {
        ERROR("BuildCatPotlist: file open error");
        return;
    }
    fprintf(fp,"%d %d %d %g\n",NCellX,NCellY,NCellZ,latticeconst[0]/4);
    INFO_Printf("%d %d %d %g\n",NCellX,NCellY,NCellZ,latticeconst[0]/4);
    for(i=0;i<_NP/4;i++)
        fprintf(fp,"%21.14e %21.14e %21.14e %21.14e\n",
                CatPot[i*4+0],CatPot[i*4+1],CatPot[i*4+2],CatPot[i*4+3]);
    fclose(fp);
    
    /* Future KMC simulations can read this file and directly work on
     *  a simple cubic lattice, without the need to deal with the YSZ
     *  crystal structure
     */
#endif
}

int YSZFrame::ReadCoulombEList()
{
    FILE *fp;
    char string[500];
    int i; double gridsize;
    
    LFile::SubHomeDir(potfile,potfile);
    
    fp=fopen(potfile,"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadCoulombEList file ("<<potfile<<") not found!");
        return -1;
    }

    fgets(string,500,fp);
    sscanf(string,"%d %d %d %lf\n",&NCellX,&NCellY,&NCellZ,&gridsize);
    INFO_Printf("%d %d %d %g\n",NCellX,NCellY,NCellZ,gridsize);

    latticeconst[0] = gridsize*4;
    INFO_Printf("latticeconst = %g\n",latticeconst[0]);
    
    //AllocLinkList();
    
    for(i=0;i<NCellX*NCellY*NCellZ/4;i++)
    {
        fgets(string,500,fp);
        sscanf(string,"%lf %lf %lf %lf",
               CoulombE+(i*4+0),CoulombE+(i*4+1),CoulombE+(i*4+2),CoulombE+(i*4+3));
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadCoulombEList done\n");
    return 0;


}

int YSZFrame::ReadCatSpeciesList()
{
    FILE *fp;
    char string[500];
    int ipt, tmpcatspecies;

    fp=fopen(catspeciesfile,"r");
    
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadCatSpeciesList file ("<<catspeciesfile<<") not found!");
        return -1;
    }
  
    for(ipt=0;ipt<_NP;ipt++) {
        fgets(string,500,fp);
        sscanf(string,"%d",&tmpcatspecies);
        if (tmpcatspecies == 2)
            species[ipt] = tmpcatspecies;
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadCatSpeciesList done\n");
    return 0;
}

int YSZFrame::ReadNNYList()
{
    FILE *fp;
    char string[500];
    int ipt, tmpnnYnum;

    Realloc(nnYlist,int,_NP);    
    memset(nnYlist,0,sizeof(int)*_NP);
    LFile::SubHomeDir(nnYlistfile,nnYlistfile);

    fp=fopen(nnYlistfile,"r");
    
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadNNYList file ("<<nnYlistfile<<") not found!");
        return -1;
    }
  
    for(ipt=0;ipt<_NP;ipt++) {
        fgets(string,500,fp);
        sscanf(string,"%d",&tmpnnYnum);
	nnYlist[ipt] = tmpnnYnum;
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadNNYList done\n");
    return 0;
}

int YSZFrame::ReadSNYList()
{
    FILE *fp;
    char string[500];
    int ipt, tmpsnYnum;

    Realloc(snYlist,int,_NP);    
    memset(snYlist,0,sizeof(int)*_NP);
    LFile::SubHomeDir(snYlistfile,snYlistfile);
    
    fp=fopen(snYlistfile,"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadSNYList file ("<<snYlistfile<<") not found!");
        return -1;
    }
  
    for(ipt=0;ipt<_NP;ipt++) {
        fgets(string,500,fp);
        sscanf(string,"%d",&tmpsnYnum);
	snYlist[ipt] = tmpsnYnum;
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadSNYList done\n");
    return 0;
}

int YSZFrame::ReadCatPotList()
{
    FILE *fp;
    char string[500];
    int i; double gridsize;
    
    LFile::SubHomeDir(catpotfile,catpotfile);
    
    fp=fopen(catpotfile,"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadCatPotList file ("<<catpotfile<<") not found!");
        return -1;
    }

    fgets(string,500,fp);
    sscanf(string,"%d %d %d %lf\n",&NCellX,&NCellY,&NCellZ,&gridsize);
    INFO_Printf("%d %d %d %g\n",NCellX,NCellY,NCellZ,gridsize);

    latticeconst[0] = gridsize*4;
    INFO_Printf("latticeconst = %g\n",latticeconst[0]);
    
    //AllocLinkList();
    
    for(i=0;i<_NP/4;i++)
    {
        fgets(string,500,fp);
        sscanf(string,"%lf %lf %lf %lf",
               CatPot+(i*4+0),CatPot+(i*4+1),CatPot+(i*4+2),CatPot+(i*4+3));
    }
    fclose(fp);
    INFO_Printf("YSZFrame::ReadCoulombEList done\n");
    return 0;
}


/********************************************************************/
/* kmc steps */
/********************************************************************/
void YSZFrame::kmcrun() {
    int i, j, k, vac_target, nei, selectnum, si, sj, newpos, mvac;
    double ExtV1, ExtV2, InteractE1, InteractE2, prob, totrate, num_rand1, num_rand2, time;
    double dx, dy, dz;
    double E_eci1, E_eci2;
//    FILE *fp_enposxyz, *fp_posx;
    FILE *fp_singledata;
    FILE *fp_vaclist, *fp_vaclist0;
    FILE *fp_rx, *fp_ry, *fp_rz;
    FILE *fp_corrs, *fp_corrtimes, *fp_corrdxs;
    FILE *fp_corrl, *fp_corrtimel, *fp_corrdxl;
    FILE *fp_treal;
    FILE *fp_histEb0;
    
    /* allocate memory */
    Realloc(corrs,double,ncorrsteps);
    Realloc(corrtimes,double,ncorrbins);
//    Realloc(corrdxs,double,num_vac*ncorrbins);
    Realloc(corrdxs,double,ncorrbins);
    Realloc(corrl,double,ncorrstepl);
    Realloc(corrtimel,double,ncorrbinl);
//    Realloc(corrdxl,double,num_vac*ncorrbinl);
    Realloc(corrdxl,double,ncorrbinl);
    Realloc(histEb0,double,3);
    Realloc(histEb,double,1000);
    vaclist0 = (int *) calloc(num_vac,sizeof(int));
    memset(corrs,0,sizeof(double)*ncorrsteps);
    memset(corrtimes,0,sizeof(double)*ncorrbins);
//    memset(corrdxs,0,sizeof(double)*num_vac*ncorrbins);
    memset(corrdxs,0,sizeof(double)*ncorrbins);
    memset(corrl,0,sizeof(double)*ncorrstepl);
    memset(corrtimel,0,sizeof(double)*ncorrbinl);
//    memset(corrdxl,0,sizeof(double)*num_vac*ncorrbinl);
    memset(corrdxl,0,sizeof(double)*ncorrbinl);
    memset(vaclist0,0,sizeof(int)*num_vac);
    memset(histEb0,0,sizeof(double)*3);
    memset(histEb,0,sizeof(double)*1000);
    
    /* transfer from real position(double) to lattice position(int) */
    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    /* make CoulombE bet. vacancies be zero*/
    if (CTL_CoulombE == 0) {
        for (i=0;i<NCellX*NCellY*NCellZ;i++)
	    *(CoulombE + i) = 0.0;
	INFO_Printf("Coulomb interaction bet V-V is turned off\n");
    }

    /* make Energy barrier be same */
    if (CTL_Ebarrier == 0) {
        for (i=0;i<_NP*NUMJUMPDIR;i++)
	    *(Ebarrier + i) = Ebconst;
	INFO_Printf("Energy barrier is set to same numbers;Eb=%e\n",Ebconst);
    }

    /* make Coulomb Potential from Cations be zero*/
    if (CTL_CatPot == 0) {
        for (i=0;i<_NP;i++)
	    *(CatPot + i) = 0.0;
	INFO_Printf("Coulomb interaction bet V-Cat is OFF\n");
    }

//    INFO_Printf("Set EnnY=%e, EsnY=%e\n",EnnY,EsnY);
    
    double kT = KB * _T;
    int unitmove = 2;       /* each vacancy moves by unitmove(int), 2 */
    double unitensmove = (double) unitmove / num_vac; /* average unitmove in view of all vacancies */
    int problength = NUMJUMPDIR*num_vac;

    double rate[problength], cumrate[problength], cumprob[problength];
    double PDInteractE[num_vac][NUMJUMPDIR];

    rx = (int *) calloc(num_vac,sizeof(int));
    ry = (int *) calloc(num_vac,sizeof(int));
    rz = (int *) calloc(num_vac,sizeof(int));
    for (i=0;i<num_vac;i++) {
        rx[i] = _RI[3*vaclist[i]+0];
        ry[i] = _RI[3*vaclist[i]+1];
        rz[i] = _RI[3*vaclist[i]+2];
    }
    enpos[0]=enpos[1]=enpos[2]=0;
    for (i=0;i<num_vac;i++) {
        enpos[0] += rx[i];
        enpos[1] += ry[i];
        enpos[2] += rz[i];
    }
    enpos[0] = enpos[0]/num_vac;
    enpos[1] = enpos[1]/num_vac;
    enpos[2] = enpos[2]/num_vac;
    
//    fprintf(fp_enposxyz,"%6d %6.6e %6.6e %6.6e\n",0,enpos[0],enpos[1],enpos[2]);
//    for (i=0;i<num_vac;i++) 
//        fprintf(fp_posx,"%d ",rx[i]);
//    fprintf(fp_posx,"\n");
//    fprintf(fp_treal,"%6.10e\n",0.0);
    treal = 0.0;

    /* Initial probability calculation: At each kmc step, only probability difference
       is passed to the next kmc step and updated by multiplying the probability
       difference to the previous probability  */
    INFO_Printf("Calculating Initial probability\n");
    E_eci1 = cal_energy_by_CEM_initial();
    INFO_Printf("Total energy by ECI at initial state = %e\n",E_eci1);
    for (i=0;i<num_vac;i++) {
        vac_target = vaclist[i];
        vaclist0[i] = vac_target;
//        InteractE1 = CalInteract(target);
//        E_eci1 = cal_energy_by_CEM();
        for (j=0;j<NUMJUMPDIR;j++) {
            nei = linklist[vac_target*LINKLEN+j+5];
//            printf("nei for %dth vac(%d) in %d direction is %d\n",i,vac_target,j,nei);
            if (species[nei] == 1) {
//                vaclist[i] = nei;                
//                InteractE2 = CalInteract(nei);
//                species[nei] = 3;
//                species[vac_target] = 1;
//                E_eci2 = cal_energy_by_CEM();
//                PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT);
                PDInteractE[i][j] = 1.0;
//                vaclist[i] = vac_target;
//                species[vac_target] = 3;
//                species[nei] = 1;
            }
            else if (species[nei] == 3) {
                PDInteractE[i][j] = 0.0;
            }
            else {
                INFO_Printf("Wrong Species Error(1)! \n");
                return;
            }
//            if (PDInteractE[i][j] == 0.0)
//                printf("Initial zero:i=%d,j=%d,species=%d,%d(%d,%d,%d)->%d(%d,%d,%d)\n",i,j,species[nei],vac_target,_RI[3*vac_target],_RI[3*vac_target+1],_RI[3*vac_target+2],nei,_RI[3*nei],_RI[3*nei+1],_RI[3*nei+2]);
        }
    }
    INFO_Printf("Initial probability calculation is done\n");

    /* start MC steps */
    iids = 0; eids = 0;  iidl = 0; eidl = 0;/*iid and eid is initialized to zero*/
    for (k=0;k<num_loop;k++) {
        for (curstep=0;curstep<num_jump;curstep++) {
            /* update the probability rate for current kmc step */
            for (i=0;i<num_vac;i++) {
                vac_target = vaclist[i];
                ExtV1 = 0;
                E_eci1 = 0.0;
                for (j=0;j<NUMJUMPDIR;j++) {
                    nei = linklist[vac_target*LINKLEN+j+5];
                    ExtV2 = 0; 
                    if (species[nei] == 1) {
                        prob = trialfreq*exp(-(Ebarrier[vac_target*NUMJUMPDIR+j]+0.5*(ExtV2-ExtV1))/kT)*exp(-(CatPot[nei]-CatPot[vac_target])/2/kT)*PDInteractE[i][j];
                        prob *= exp(-EnnY*(nnYlist[nei]-nnYlist[vac_target])/2/kT)*exp(-EsnY*(snYlist[nei]-snYlist[vac_target])/2/kT);
                        /*** cluster_expansion_method correction ***/
                        vaclist[i] = nei;
                        E_eci2 = cal_energy_by_CEM_update(vac_target,nei);
//                        INFO_Printf("i=%d, j=%d, dE_eci=%e\n",i,j,E_eci2-E_eci1);
//                        INFO_Printf("V(%d,%d,%d)->V(%d,%d,%d), dE_eci=%e\n",_RI[3*vac_target+0],_RI[3*vac_target+1],_RI[3*vac_target+2],_RI[3*nei+0],_RI[3*nei+1],_RI[3*nei+2],E_eci2-E_eci1);
                        vaclist[i] = vac_target;
                        prob *= exp(-(E_eci2-E_eci1)/2/kT);                        
                        /*******************************************/
//                    printf("PDInteractE[%d][%d]=%e\n",i,j,PDInteractE[i][j]);
//                    prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]/kT);
//                    if ( PDInteractE[i][j] == 0 ) {
//                        printf("curstep=%d, prob[%d][%d]=0, vac_target=%d(%d,%d,%d),nei=%d(%d,%d,%d),species=%d\n",curstep,i,j,vac_target,_RI[3*vac_target],_RI[3*vac_target+1],_RI[3*vac_target+2],nei,_RI[3*nei],_RI[3*nei+1],_RI[3*nei+2],species[nei]);
//                        return;
//                    }                        
                    }
                    else if (species[nei] == 3)
                        prob = 0.0;
                    else {
                        INFO_Printf("Wrong Species Error(2)! \n");
                        return;
                    }
                    rate[j+6*i] = prob;
                }
            }

            /* generate cumulated probability list */
            cumrate[0] = rate[0];
            for (i=1;i<problength;i++) {
                cumrate[i] = cumrate[i-1] + rate[i];
            }
            totrate = cumrate[problength-1];
            for (i=0;i<problength;i++) 
                cumprob[i] = cumrate[i] / totrate;
            
            /* event selection */
            num_rand1 = drand48();
            time = -1*log(num_rand1) / totrate;
            num_rand2 = drand48();
            selectnum = -1;
            for (i=0;i<problength;i++)
                if (cumprob[i] > num_rand2) {
                    selectnum = i;
                    break;
                }
            si = (int) floor(selectnum/NUMJUMPDIR);  /* determine vacancy ID */
            sj = selectnum - NUMJUMPDIR*si;          /* determine moving direction */
            treal += time;                           /* update the real time */
            
            if((curstep%printfreq)==0)
                INFO_Printf("%dth kmcstep: vac:%d, dir:%d \n",curstep+1,si,sj);
            
            /* update vacancy list, vaclist, and update the ensemble average position
               of total vacancy */
            vac_target = vaclist[si];
            newpos = linklist[vac_target*LINKLEN+sj+5];
            if (species[vac_target] == 1) {
                INFO_Printf("Wrong Species Error(3)!\n");
                return;
            }
            species[vac_target] = 1;
            if (species[newpos] == 3) {
                INFO_Printf("Wrong Species Error(4)!\n");
                return;
            }
            species[newpos] = 3;
            vaclist[si] = newpos;
            /* update cluster spins */
            cluster_spin[vac_target] = 1;
            cluster_spin[newpos] = -1;
            
//        printf("curstep=%d, moved #%d:(%d,%d,%d) -> #%d:(%d,%d,%d)\n",curstep,vac_target,_RI[3*vac_target],_RI[3*vac_target+1],_RI[3*vac_target+2],newpos,_RI[3*newpos],_RI[3*newpos+1],_RI[3*newpos+2]);
            switch (sj) {
            case 0: {rx[si] += unitmove; dx =  unitensmove; dy = 0.0; dz = 0.0; break;}
            case 1: {rx[si] -= unitmove; dx = -unitensmove; dy = 0.0; dz = 0.0; break;}
            case 2: {ry[si] += unitmove; dx = 0.0; dy =  unitensmove; dz = 0.0; break;}
            case 3: {ry[si] -= unitmove; dx = 0.0; dy = -unitensmove; dz = 0.0; break;}
            case 4: {rz[si] += unitmove; dx = 0.0; dy = 0.0; dz =  unitensmove; break;}
            case 5: {rz[si] -= unitmove; dx = 0.0; dy = 0.0; dz = -unitensmove; break;}
            default: dx = dy = dz = 0.0; FATAL("unrecognized sj"<<sj);
            }
            enpos[0] += dx;
            enpos[1] += dy;
            enpos[2] += dz;
            CorrFcnX(dx);
//        CorrFcnXn(si,sj);
#if 0
            /* update the probability difference */
            for (i=0;i<num_vac;i++) {
                if (i == si) {
                    InteractE1 = CalInteract(newpos);
                    for (j=0;j<NUMJUMPDIR;j++) {
                        nei = linklist[newpos*LINKLEN+5+j];
                        if (species[nei] == 1) {
                            vaclist[i] = nei;
                            InteractE2 = CalInteract(nei);
                            PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                            vaclist[i] = newpos;
                        }
                        else if(species[nei] == 3)
                            PDInteractE[i][j] = 0;                
                        else {
                            INFO_Printf("Wrong Species Error(5)!\n");
                            return;
                        }
                    }
                }
                else {
                    mvac = vaclist[i];
                    for (j=0;j<NUMJUMPDIR;j++) {
                        nei = linklist[mvac*LINKLEN+5+j];
                        if (nei == vac_target) {
                            InteractE1 = CalInteract(mvac);
                            vaclist[i] = nei;
                            InteractE2 = CalInteract(nei);
                            PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                            vaclist[i] = mvac;
                        }
                        else                        
                            PDInteractE[i][j] *= PDiffInteract(mvac,j,vac_target,newpos);
                    }
                }
            }
#endif
//        fprintf(fp_enposxyz,"%6d %6.6e %6.6e %6.6e\n",curstep+1,enpos[0],enpos[1],enpos[2]);
//        for (i=0;i<num_vac;i++)
//            fprintf(fp_posx,"%d ",rx[i]);
//        fprintf(fp_posx,"\n");       
//        fprintf(fp_treal,"%6.10e\n",treal);


            if(savecorr)
            {
                if(((curstep+1)%savecorrfreq)==0)
                {
                    cfs.write(this);
                    cfl.write(this);
                    fp_treal = fopen("treal.dat","w");
                    fprintf(fp_treal,"%6.10e\n",treal);
                    fclose(fp_treal);
                }
            }   
            
            if(savecn)
            {
                if(((curstep+1)%savecnfreq)==0&&(curstep!=0))
                    intercn.write(this);            
            }
            
            if(savediff)
            {
                if(((curstep+1)%savedifffreq)==0) {
                    df.write(this);
                }
            }

            if(savehist)
            {
                if(((curstep+1)%savehistfreq)==0) {
                    hf.write(this);
                }
                if(((curstep+1)%savehistfreq)==0) {
                    fp_histEb0 = fopen("histEb0.dat","a");
                    fprintf(fp_histEb0,"%e %e %e\n",histEb0[0],histEb0[1],histEb0[2]);
                }
            }

            if(saveatomeye)
            {
                if(((curstep+1)%saveatomeyefreq)==0) {
                    mf.write(this);
                }
            }
            
            if(savecontdata)
            {
                if (((curstep+1)%savecontfreq)==0) {
                    fp_singledata = fopen("singledata.dat","w");
                    fp_corrs = fopen("corrs.dat","w");
                    fp_corrtimes = fopen("corrtimes.dat","w");
                    fp_corrdxs = fopen("corrdxs.dat","w");
                    fp_corrl = fopen("corrl.dat","w");
                    fp_corrtimel = fopen("corrtimel.dat","w");
                    fp_corrdxl = fopen("corrdxl.dat","w");
                    fp_rx = fopen("rx.dat","w");
                    fp_ry = fopen("ry.dat","w");
                    fp_rz = fopen("rz.dat","w");
                    fp_vaclist = fopen("vaclist.dat","w");
                    fp_vaclist0 = fopen("vaclist0.dat","w");
                    
                    fprintf(fp_singledata,"%d %d %d %d %d %d %e\n",k,curstep,iids,eids,iidl,eidl,treal);
                    for (i=0;i<ncorrsteps;i++) {
                        fprintf(fp_corrs,"%e\n",corrs[i]);
                    }
                    for (i=0;i<ncorrbins;i++) {
                        fprintf(fp_corrtimes,"%e\n",corrtimes[i]);
                        fprintf(fp_corrdxs,"%e\n",corrdxs[i]);       
                    }
                    for (i=0;i<ncorrstepl;i++) {
                        fprintf(fp_corrl,"%e\n",corrl[i]);
                    }
                    for (i=0;i<ncorrbinl;i++) {
                        fprintf(fp_corrtimel,"%e\n",corrtimel[i]);
                        fprintf(fp_corrdxl,"%e\n",corrdxl[i]);       
                    }
                    for (i=0;i<num_vac;i++) {
                        fprintf(fp_rx,"%d\n",rx[i]);
                        fprintf(fp_ry,"%d\n",ry[i]);
                        fprintf(fp_rz,"%d\n",rz[i]);
                        fprintf(fp_vaclist,"%d\n",vaclist[i]);
                        fprintf(fp_vaclist0,"%d\n",vaclist0[i]);
                    }
                    fclose(fp_singledata);
                    fclose(fp_corrs);
                    fclose(fp_corrtimes);
                    fclose(fp_corrdxs);
                    fclose(fp_corrl);
                    fclose(fp_corrtimel);
                    fclose(fp_corrdxl);
                    fclose(fp_rx);
                    fclose(fp_ry);
                    fclose(fp_rz);
                    fclose(fp_vaclist);
                    fclose(fp_vaclist0);
                    INFO_Printf("Saved all data\n");
                }
            }
        }       
    }

//    INFO_Printf("iid=%d, eid=%d, ncorrstep=%d\n",iid,eid,ncorrstep);
//    fclose(fp_enposxyz);
//    fclose(fp_posx);
    
}

/****************************************************************/
/* kmcrun, keep the part only after transient period            */
/****************************************************************/
void YSZFrame::kmcruneq() {
    int i, j, k, l, vac_target, nei, selectnum, si, sj, newpos, mvac;
    double ExtV1, ExtV2, InteractE1, InteractE2, prob, totrate, num_rand1, num_rand2, time;
    double dx, dy, dz;
    double E_eci1, E_eci2, dE_eci;
    FILE *fp_singledata;
    FILE *fp_vaclist, *fp_vaclist0;
    FILE *fp_rx, *fp_ry, *fp_rz;
    FILE *fp_corrs, *fp_corrtimes, *fp_corrdxs;
    FILE *fp_corrl, *fp_corrtimel, *fp_corrdxl;
    FILE *fp_treal;
    FILE *fp_histEb0;
    
    /* allocate memory */
    Realloc(corrs,double,ncorrsteps);
    Realloc(corrtimes,double,ncorrbins);
    Realloc(corrdxs,double,ncorrbins);
    Realloc(corrl,double,ncorrstepl);
    Realloc(corrtimel,double,ncorrbinl);
    Realloc(corrdxl,double,ncorrbinl);
    Realloc(histEb0,double,3);
    Realloc(histEb,double,1000);
    vaclist0 = (int *) calloc(num_vac,sizeof(int));
    memset(corrs,0,sizeof(double)*ncorrsteps);
    memset(corrtimes,0,sizeof(double)*ncorrbins);
    memset(corrdxs,0,sizeof(double)*ncorrbins);
    memset(corrl,0,sizeof(double)*ncorrstepl);
    memset(corrtimel,0,sizeof(double)*ncorrbinl);
    memset(corrdxl,0,sizeof(double)*ncorrbinl);
    memset(vaclist0,0,sizeof(int)*num_vac);
    memset(histEb0,0,sizeof(double)*3);
    memset(histEb,0,sizeof(double)*1000);

    /* transfer from real position(double) to lattice position(int) */
    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    /* make CoulombE bet. vacancies be zero*/
    if (CTL_CoulombE == 0) {
        for (i=0;i<NCellX*NCellY*NCellZ;i++)
	    *(CoulombE + i) = 0.0;
	INFO_Printf("Coulomb interaction bet V-V is turned off\n");
    }

    /* make Energy barrier be same */
    if (CTL_Ebarrier == 0) {
        for (i=0;i<_NP*NUMJUMPDIR;i++)
	    *(Ebarrier + i) = Ebconst;
	INFO_Printf("Energy barrier is set to same numbers;Eb=%e\n",Ebconst);
    }

    /* make Coulomb Potential from Cations be zero*/
    if (CTL_CatPot == 0) {
        for (i=0;i<_NP;i++)
	    *(CatPot + i) = 0.0;
	INFO_Printf("Coulomb interaction bet V-Cat is OFF\n");
    }
    
    double kT = KB * _T;
    int unitmove = 2;       /* each vacancy moves by unitmove(int), 2 */
    double unitensmove = (double) unitmove / num_vac; /* average unitmove in view of all vacancies */
    int problength = NUMJUMPDIR*num_vac;

    double rate[problength], cumrate[problength], cumprob[problength];
    double PDInteractE[num_vac][NUMJUMPDIR];

    E_eci1 = cal_energy_by_CEM_initial();
    INFO_Printf("Total energy by ECI at initial state = %e\n",E_eci1);
    INFO_Printf("Transient period starts and proceedes to %d steps\n",ntransient);
    /* transient period */
    for (curstep=0;curstep<ntransient;curstep++) {
        /* update the probability rate for current kmc step */
        for (i=0;i<num_vac;i++) {
            vac_target = vaclist[i];
            ExtV1 = 0;
            E_eci1 = 0.0;
            for (j=0;j<NUMJUMPDIR;j++) {
                nei = linklist[vac_target*LINKLEN+j+5];
                ExtV2 = 0; 
                if (species[nei] == 1) {
                    /*** cluster_expansion_method correction ***/
                    E_eci2 = cal_energy_by_CEM_update(vac_target,nei);
                    dE_eci = E_eci2 - E_eci1;
                    if (dE_eci < Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[0])
                        prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]*ebcorrtable[0]/kT);
                    else if (dE_eci >= Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[nebcorrtable-1])
                        prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]*ebcorrtable[nebcorrtable-1]/kT);
                    else
                        for (k=0;k<nebcorrtable-1;k++) 
                            if (dE_eci >= Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[k] && dE_eci < Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[k+1]) {
                                prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]*ebcorrtable[k]/kT);
                                break;
                            }
                    /*******************************************/
                }
                else if (species[nei] == 3)
                    prob = 0.0;
                else {
                    INFO_Printf("Wrong Species Error(2)! \n");
                    return;
                }
                if (prob > 1.0*trialfreq)
                    prob = 1.0*trialfreq;
                rate[j+6*i] = prob;
            }
        }
        
        /* generate cumulated probability list */
        cumrate[0] = rate[0];
        for (i=1;i<problength;i++) {
            cumrate[i] = cumrate[i-1] + rate[i];
        }
        totrate = cumrate[problength-1];
        for (i=0;i<problength;i++) 
            cumprob[i] = cumrate[i] / totrate;
        
        /* event selection */
        num_rand1 = drand48();
        time = -1*log(num_rand1) / totrate;
        num_rand2 = drand48();
        selectnum = -1;
        for (i=0;i<problength;i++)
            if (cumprob[i] > num_rand2) {
                selectnum = i;
                break;
            }
        si = (int) floor(selectnum/NUMJUMPDIR);  /* determine vacancy ID */
        sj = selectnum - NUMJUMPDIR*si;          /* determine moving direction */
        
        if((curstep%printfreq)==0)
            INFO_Printf("%dth kmcstep: vac:%d, dir:%d \n",curstep+1,si,sj);
        
        /* update vacancy list, vaclist, and update the ensemble average position
           of total vacancy */
        vac_target = vaclist[si];
        newpos = linklist[vac_target*LINKLEN+sj+5];
        if (species[vac_target] == 1) {
            INFO_Printf("Wrong Species Error(3)!\n");
            return;
        }
        species[vac_target] = 1;
        if (species[newpos] == 3) {
            INFO_Printf("Wrong Species Error(4)!\n");
            return;
        }
        species[newpos] = 3;
        vaclist[si] = newpos;
        /* update cluster spins */
        cluster_spin[vac_target] = 1;
        cluster_spin[newpos] = -1;

    }    
    INFO_Printf("Transient period finished. Eq simulation starts\n");

    rx = (int *) calloc(num_vac,sizeof(int));
    ry = (int *) calloc(num_vac,sizeof(int));
    rz = (int *) calloc(num_vac,sizeof(int));
    for (i=0;i<num_vac;i++) {
        rx[i] = _RI[3*vaclist[i]+0];
        ry[i] = _RI[3*vaclist[i]+1];
        rz[i] = _RI[3*vaclist[i]+2];
        vaclist0[i] = vaclist[i];
    }
    enpos[0]=enpos[1]=enpos[2]=0;
    for (i=0;i<num_vac;i++) {
        enpos[0] += rx[i];
        enpos[1] += ry[i];
        enpos[2] += rz[i];
    }
    enpos[0] = enpos[0]/num_vac;
    enpos[1] = enpos[1]/num_vac;
    enpos[2] = enpos[2]/num_vac;
    
    treal = 0.0;

    /* start equilibrium MC steps */
    iids = 0; eids = 0;  iidl = 0; eidl = 0;/*iid and eid is initialized to zero*/
    for (k=0;k<num_loop;k++) {
        for (curstep=0;curstep<num_jump;curstep++) {
            /* update the probability rate for current kmc step */
            for (i=0;i<num_vac;i++) {
                vac_target = vaclist[i];
                ExtV1 = 0;
                E_eci1 = 0.0;
                for (j=0;j<NUMJUMPDIR;j++) {
                    nei = linklist[vac_target*LINKLEN+j+5];
                    ExtV2 = 0; 
                    if (species[nei] == 1) {
                        /*** cluster_expansion_method correction ***/
                        E_eci2 = cal_energy_by_CEM_update(vac_target,nei);
                        dE_eci = E_eci2 - E_eci1;
                        if (dE_eci < Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[0])
                            prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]*ebcorrtable[0]/kT);
                        else if (dE_eci >= Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[nebcorrtable-1])
                            prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]*ebcorrtable[nebcorrtable-1]/kT);
                        else
                            for (l=0;l<nebcorrtable-1;l++) 
                                if (dE_eci >= Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[l] && dE_eci < Ebarrier[vac_target*NUMJUMPDIR+j]*ebxcorrtable[l+1]) {
                                    prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]*ebcorrtable[l]/kT);
                                    break;
                                }
                        /*******************************************/
                    }
                    else if (species[nei] == 3)
                        prob = 0.0;
                    else {
                        INFO_Printf("Wrong Species Error(2)! \n");
                        return;
                    }
                    if (prob > 1.0*trialfreq)
                        prob = 1.0*trialfreq;
                    rate[j+6*i] = prob;
                }
            }

            /* generate cumulated probability list */
            cumrate[0] = rate[0];
            for (i=1;i<problength;i++) {
                cumrate[i] = cumrate[i-1] + rate[i];
            }
            totrate = cumrate[problength-1];
            for (i=0;i<problength;i++) 
                cumprob[i] = cumrate[i] / totrate;
            
            /* event selection */
            num_rand1 = drand48();
            time = -1*log(num_rand1) / totrate;
            num_rand2 = drand48();
            selectnum = -1;
            for (i=0;i<problength;i++)
                if (cumprob[i] > num_rand2) {
                    selectnum = i;
                    break;
                }
            si = (int) floor(selectnum/NUMJUMPDIR);  /* determine vacancy ID */
            sj = selectnum - NUMJUMPDIR*si;          /* determine moving direction */
            treal += time;                           /* update the real time */
            
            if((curstep%printfreq)==0)
                INFO_Printf("%dth kmcstep: vac:%d, dir:%d \n",curstep+1,si,sj);
            
            /* update vacancy list, vaclist, and update the ensemble average position
               of total vacancy */
            vac_target = vaclist[si];
            newpos = linklist[vac_target*LINKLEN+sj+5];
            if (species[vac_target] == 1) {
                INFO_Printf("Wrong Species Error(3)!\n");
                return;
            }
            species[vac_target] = 1;
            if (species[newpos] == 3) {
                INFO_Printf("Wrong Species Error(4)!\n");
                return;
            }
            species[newpos] = 3;
            vaclist[si] = newpos;
            /* update cluster spins */
            cluster_spin[vac_target] = 1;
            cluster_spin[newpos] = -1;

//            switch(Ebarrier[vac_target*NUMJUMPDIR+sj]) {
//            case(EbZZ_Car) : histEb0[0] += 1.0; break;
//            case(EbZY_Car) : histEb0[1] += 1.0; break;
//            case(EbYY_Car) : histEb0[2] += 1.0; break;
//            }

            if (sj < 2) {
                if (Ebarrier[vac_target*NUMJUMPDIR+sj]==EbZZ_Car)
                    histEb0[0] += 1.0;
                else if (Ebarrier[vac_target*NUMJUMPDIR+sj]==EbZY_Car)
                    histEb0[1] += 1.0;
                else if (Ebarrier[vac_target*NUMJUMPDIR+sj]==EbYY_Car)
                    histEb0[2] += 1.0;
                else
                    ERROR("No Eb0 matched! Eb0="<<Ebarrier[vac_target*NUMJUMPDIR+sj]);        
                
                for (i=0;i<1000;i++) 
                    if (rate[sj+6*si] > trialfreq*exp(-(0.02*i-10)/kT)) {
                        histEb[i] += 1.0;
                        break;
                    }
            }
            switch (sj) {
            case 0: {rx[si] += unitmove; dx =  unitensmove; dy = 0.0; dz = 0.0; break;}
            case 1: {rx[si] -= unitmove; dx = -unitensmove; dy = 0.0; dz = 0.0; break;}
            case 2: {ry[si] += unitmove; dx = 0.0; dy =  unitensmove; dz = 0.0; break;}
            case 3: {ry[si] -= unitmove; dx = 0.0; dy = -unitensmove; dz = 0.0; break;}
            case 4: {rz[si] += unitmove; dx = 0.0; dy = 0.0; dz =  unitensmove; break;}
            case 5: {rz[si] -= unitmove; dx = 0.0; dy = 0.0; dz = -unitensmove; break;}
            default: dx = dy = dz = 0.0; FATAL("unrecognized sj"<<sj);
            }
            enpos[0] += dx;
            enpos[1] += dy;
            enpos[2] += dz;
            CorrFcnX(dx);

            if(savecorr)
            {
                if(((curstep+1)%savecorrfreq)==0)
                {
                    cfs.write(this);
                    cfl.write(this);
                    fp_treal = fopen("treal.dat","w");
                    fprintf(fp_treal,"%6.10e\n",treal);
                    fclose(fp_treal);
                }
            }   
            
            if(savecn)
            {
                if(((curstep+1)%savecnfreq)==0&&(curstep!=0))
                    intercn.write(this);            
            }
            
            if(savediff)
            {
                if(((curstep+1)%savedifffreq)==0) {
                    df.write(this);
                }
            }

            if(savehist)
            {
                if(((curstep+1)%savehistfreq)==0) {
                    hf.write(this);
                }
                if(((curstep+1)%savehistfreq)==0) {
                    fp_histEb0 = fopen("histEb0.dat","a");
                    fprintf(fp_histEb0,"%e %e %e\n",histEb0[0],histEb0[1],histEb0[2]);
                }
            }

            if(saveatomeye)
            {
                if(((curstep+1)%saveatomeyefreq)==0) {
                    mf.write(this);
                }
            }
            
            if(savecontdata)
            {
                if (((curstep+1)%savecontfreq)==0) {
                    fp_singledata = fopen("singledata.dat","w");
                    fp_corrs = fopen("corrs.dat","w");
                    fp_corrtimes = fopen("corrtimes.dat","w");
                    fp_corrdxs = fopen("corrdxs.dat","w");
                    fp_corrl = fopen("corrl.dat","w");
                    fp_corrtimel = fopen("corrtimel.dat","w");
                    fp_corrdxl = fopen("corrdxl.dat","w");
                    fp_rx = fopen("rx.dat","w");
                    fp_ry = fopen("ry.dat","w");
                    fp_rz = fopen("rz.dat","w");
                    fp_vaclist = fopen("vaclist.dat","w");
                    fp_vaclist0 = fopen("vaclist0.dat","w");
                    
                    fprintf(fp_singledata,"%d %d %d %d %d %d %e\n",k,curstep,iids,eids,iidl,eidl,treal);
                    for (i=0;i<ncorrsteps;i++) {
                        fprintf(fp_corrs,"%e\n",corrs[i]);
                    }
                    for (i=0;i<ncorrbins;i++) {
                        fprintf(fp_corrtimes,"%e\n",corrtimes[i]);
                        fprintf(fp_corrdxs,"%e\n",corrdxs[i]);       
                    }
                    for (i=0;i<ncorrstepl;i++) {
                        fprintf(fp_corrl,"%e\n",corrl[i]);
                    }
                    for (i=0;i<ncorrbinl;i++) {
                        fprintf(fp_corrtimel,"%e\n",corrtimel[i]);
                        fprintf(fp_corrdxl,"%e\n",corrdxl[i]);       
                    }
                    for (i=0;i<num_vac;i++) {
                        fprintf(fp_rx,"%d\n",rx[i]);
                        fprintf(fp_ry,"%d\n",ry[i]);
                        fprintf(fp_rz,"%d\n",rz[i]);
                        fprintf(fp_vaclist,"%d\n",vaclist[i]);
                        fprintf(fp_vaclist0,"%d\n",vaclist0[i]);
                    }
                    fclose(fp_singledata);
                    fclose(fp_corrs);
                    fclose(fp_corrtimes);
                    fclose(fp_corrdxs);
                    fclose(fp_corrl);
                    fclose(fp_corrtimel);
                    fclose(fp_corrdxl);
                    fclose(fp_rx);
                    fclose(fp_ry);
                    fclose(fp_rz);
                    fclose(fp_vaclist);
                    fclose(fp_vaclist0);
                    INFO_Printf("Saved all data\n");
                }
            }                
        }       
    }

    
}


/****************************************************************/
/*kmcsrun with external potential */
/****************************************************************/

void YSZFrame::kmcrunAC() {
    int i, j, vac_target, nei, selectnum, si, sj, newpos, mvac;
    double ExtV1, ExtV2, InteractE1, InteractE2, prob, totrate, num_rand1, num_rand2, time, oldtreal;
    double dx, dy, dz;
    double E_eci1, E_eci2;
//    FILE *fp_enposxyz;
//    FILE *fp_treal;
    int nr, cid, nid;
    
    nextdt = (int) round(2*M_PI/extdt/extw0);

    /* allocate memory */
    Realloc(enposwext,double,nextdt);
    memset(enposwext,0,sizeof(double)*nextdt);
//    fp_enposxyz = fopen("enposxyz.dat","w");
//    fp_treal = fopen("treal.dat","w");

    /* transfer from real position(double) to lattice position(int) */
    for (int i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    /* make CoulombE be zero*/
    if (CTL_CoulombE == 0) {
        for (i=0;i<NCellX*NCellY*NCellZ;i++)
	    *(CoulombE + i) = 0.0;
	INFO_Printf("Coulomb Energy is turned off\n");
    }

    /* make Energy barrier be same */
    if (CTL_Ebarrier == 0) {
        for (i=0;i<_NP*NUMJUMPDIR;i++)
	    *(Ebarrier + i) = 0.1;
	INFO_Printf("Energy barrier is set to same numbers\n");
    }

        /* make Coulomb Potential from Cations be zero*/
    if (CTL_CatPot == 0) {
        for (i=0;i<_NP;i++)
	    *(CatPot + i) = 0.0;
	INFO_Printf("Coulomb Potential due to Cations is OFF\n");
    }

    double kT = KB * _T;
    int unitmove = 2;       /* each vacancy moves by unitmove(int), 2 */
    double unitensmove = (double) unitmove / num_vac; /* average unitmove in view of all vacancies */
    int problength = NUMJUMPDIR*num_vac;

    double rate[problength], cumrate[problength], cumprob[problength];
    double PDInteractE[num_vac][NUMJUMPDIR];

    rx = (int *) calloc(num_vac,sizeof(int));
    ry = (int *) calloc(num_vac,sizeof(int));
    rz = (int *) calloc(num_vac,sizeof(int));
    for (i=0;i<num_vac;i++) {
        rx[i] = _RI[3*vaclist[i]+0];
        ry[i] = _RI[3*vaclist[i]+1];
        rz[i] = _RI[3*vaclist[i]+2];
    }
    enpos[0]=enpos[1]=enpos[2]=0;
    for (i=0;i<num_vac;i++) {
        enpos[0] += rx[i];
        enpos[1] += ry[i];
        enpos[2] += rz[i];
    }
    enpos[0] = enpos[0]/num_vac;
    enpos[1] = enpos[1]/num_vac;
    enpos[2] = enpos[2]/num_vac;
    enposwext[0] = enpos[0];

    treal = 0.0;
    
//    fprintf(fp_enposxyz,"%6d %6.6e %6.6e %6.6e\n",0,enpos[0],enpos[1],enpos[2]);
//    fprintf(fp_treal,"%6.10e\n",0.0);

    /* Initial probability calculation: At each kmc step, only probability difference
       is passed to the next kmc step and updated by multiplying the probability
       difference to the previous probability  */
    INFO_Printf("Initial probability calculation\n");
    for (i=0;i<num_vac;i++) {
        vac_target = vaclist[i];
        InteractE1 = CalInteract(vac_target);
        for (j=0;j<NUMJUMPDIR;j++) {
            nei = linklist[vac_target*LINKLEN+j+5];
            if (species[nei] == 1) {
                vaclist[i] = nei;
                InteractE2 = CalInteract(nei);
                PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                vaclist[i] = vac_target;
            }
            else if (species[nei] == 3) {
                PDInteractE[i][j] = 0.0;
            }
            else {
                INFO_Printf("Wrong Species Error(1)! \n");
                return;
            }
        }
    }
    

    /* start MC steps */
    cid = 0; nr = 0;  /* cid and nr is initialized to zero */
    while (1!=0) {
        /* update the probability rate for current kmc step */
        for (i=0;i<num_vac;i++) {
            vac_target = vaclist[i];
            ExtV1 = 0;
            E_eci1 = cal_energy_by_CEM();
            for (j=0;j<NUMJUMPDIR;j++) {
                nei = linklist[vac_target*LINKLEN+j+5];
                if (j==0 || j==1)
                    ExtV2 = extU0*sin(extw0*treal)*(4*j-2)/NCellX;
                else
                    ExtV2 = 0;
                if (species[nei] == 1) {
                    prob = trialfreq*exp(-(Ebarrier[vac_target*NUMJUMPDIR+j]+0.5*(ExtV2-ExtV1))/kT)*exp(-(CatPot[nei]-CatPot[vac_target])/2/kT)*PDInteractE[i][j];
                    prob *= exp(-EnnY*(nnYlist[nei]-nnYlist[vac_target])/2/kT)*exp(-EsnY*(snYlist[nei]-snYlist[vac_target])/2/kT);
                    /*** cluster_expansion_method correction ***/
                    vaclist[i] = nei;
                    species[nei] = 3;
                    species[vac_target] = 1;
                    E_eci2 = cal_energy_by_CEM();
                    vaclist[i] = vac_target;
                    species[nei] = 1;
                    species[vac_target] = 3;
                    /*******************************************/
		}
                else if (species[nei] == 3)
                    prob = 0.0;
                else {
                    INFO_Printf("Wrong Species Error(2)! \n");
                    return;
                }
                rate[j+6*i] = prob;
            }
        }

        /* generate cumulated probability list */
        cumrate[0] = rate[0];
        for (i=1;i<problength;i++) {
            cumrate[i] = cumrate[i-1] + rate[i];
        }
        totrate = cumrate[problength-1];
        for (i=0;i<problength;i++)
            cumprob[i] = cumrate[i] / totrate;

        /* event selection */
        num_rand1 = drand48();
        time = -1*log(num_rand1) / totrate;
        oldtreal = treal;
        if (time < timelimit) {
            num_rand2 = drand48();
            selectnum = -1;
            for (i=0;i<problength;i++)
                if (cumprob[i] > num_rand2) {
                    selectnum = i;
                    break;
                }
            si = (int) floor(selectnum/NUMJUMPDIR);  /* determine vacancy ID */
            sj = selectnum - NUMJUMPDIR*si;          /* determine moving direction */
            
            /* update vacancy list, vaclist, and update the ensemble average position
               of total vacancy */
            vac_target = vaclist[si];
            newpos = linklist[vac_target*LINKLEN+sj+5];
            if (species[vac_target] == 1) {
                INFO_Printf("Wrong Species Error(3)!\n");
                return;
            }
            species[vac_target] = 1;
            if (species[newpos] == 3) {
                INFO_Printf("Wrong Species Error(4)!\n");
                return;
            }
            species[newpos] = 3;
            vaclist[si] = newpos;
            switch (sj) {
            case 0: {rx[si] += unitmove; dx =  unitensmove; dy = 0.0; dz = 0.0; break;}
            case 1: {rx[si] -= unitmove; dx = -unitensmove; dy = 0.0; dz = 0.0; break;}
            case 2: {ry[si] += unitmove; dx = 0.0; dy =  unitensmove; dz = 0.0; break;}
            case 3: {ry[si] -= unitmove; dx = 0.0; dy = -unitensmove; dz = 0.0; break;}
            case 4: {rz[si] += unitmove; dx = 0.0; dy = 0.0; dz =  unitensmove; break;}
            case 5: {rz[si] -= unitmove; dx = 0.0; dy = 0.0; dz = -unitensmove; break;}
            default: dx = dy = dz = 0.0; FATAL("unrecoganized sj"<<sj);
            }
            enpos[0] += dx;
            enpos[1] += dy;
            enpos[2] += dz;

            for (i=0;i<num_vac;i++) {
                if (i == si) {
                    InteractE1 = CalInteract(newpos);
                    for (j=0;j<NUMJUMPDIR;j++) {
                        nei = linklist[newpos*LINKLEN+5+j];
                        if (species[nei] == 1) {
                            vaclist[i] = nei;
                            InteractE2 = CalInteract(nei);
                            PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                            vaclist[i] = newpos;
                        }
                        else if(species[nei] == 3)
                            PDInteractE[i][j] = 0;                
                        else {
                            INFO_Printf("Wrong Species Error(5)!\n");
                            return;
                        }
                    }
                }
                else {
                    mvac = vaclist[i];
                    for (j=0;j<NUMJUMPDIR;j++) {
                        nei = linklist[mvac*LINKLEN+5+j];
                        if (nei == vac_target) {
                            InteractE1 = CalInteract(mvac);
                            vaclist[i] = nei;
                            InteractE2 = CalInteract(nei);
                            PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                            vaclist[i] = mvac;
                        }
                        else                        
                            PDInteractE[i][j] *= PDiffInteract(mvac,j,vac_target,newpos);
                    }
                }
            }
            treal += time;
//            INFO_Printf("time is updated\n");
        }
        else {
            treal += timelimit;
//            INFO_Printf("timelimit is added, treal=%e,oldtreal=%e\n",treal,oldtreal);
        }
#if 0
	double cycleid = (double) nr*nextdt;
        nid = (int) floor(treal/extdt-cycleid-cid);
	if( ! ((nid>=0)&&(nid<nextdt)) )
	{
            FATAL("nr="<<nr<<",nid="<<nid<<",treal="<<treal<<",oldtreal="<<oldtreal<<",cid="<<cid);
        }
        int numcycle = (int) floor(nid/nextdt);
        if (numcycle > 0)
            ERROR("numcycle is not 0, but"<<numcycle);
        /* update the probability difference */
        for (k=0;k<numcycle;k++) {
            for (i=cid+1;i<nextdt;i++) {
                enposwext[i] += enpos[0];
            }
            nr++;
            if (nr%printfreq == 0)
                INFO_Printf("%d cycles were completed!(but oversized)\n",nr);
            if (nr%saveacfreq == 0)
                af.write(this);
            if (nr == ncycle) {
                INFO_Printf("final cycle was completed!\n");                    
                return;
            }
            for (i=0;i<cid+1;i++) {
                enposwext[i] += enpos[0];
            }
        }
        nid = nid - numcycle*nextdt;
        
        if (nid != 0) {
            if (cid+nid > nextdt-1) {         /*previous "if (cid+nid > nextdt-1)" was fatal bug??? */
                for (i=cid+1;i<nextdt;i++) {
                    enposwext[i] += enpos[0];
                }
                nr++;
                if (nr%printfreq == 0)
                    INFO_Printf("%d cycles were completed! treal=%e oldtreal=%e nid=%d cid=%d\n",nr,treal,oldtreal,nid,cid);
                if (nr%saveacfreq == 0)
                    af.write(this);
                if (nr == ncycle) {
                    INFO_Printf("final cycle was completed!\n");                    
                    return;
                }
                for (i=0;i<cid+nid-nextdt+1;i++) {
                    enposwext[i] += enpos[0];
                }
                cid = cid + nid - nextdt;
            }
            else {
                for (i=cid+1;i<cid+nid+1;i++) {
                    enposwext[i] += enpos[0];
                }
                cid = cid + nid;
            }            
        }
#endif
#if 0
        nid = (int) floor(treal/extdt-cid);
	if( ! ((nid>=0)&&(nid<nextdt)) )
	{
            FATAL("nr="<<nr<<",nid="<<nid<<",treal="<<treal<<",time="<<time<<",oldtreal="<<oldtreal<<",cid="<<cid);
        }
        int jid = 0;
        if (nid != 0) 
            for (i=1;i<nid+1;i++) {
                cid++;
                jid = cid%nextdt;
                enposwext[jid] += enpos[0];
                if (jid == nextdt-1) {
                    nr++;
                    if (nr%printfreq == 0)
                        INFO_Printf("%d cycles were completed! treal=%e nid=%d cid=%d jid=%d\n",nr,treal,nid,cid,jid);
                    if (nr%saveacfreq == 0)
                        af.write(this);
                    if (nr == ncycle) {
                        INFO_Printf("final cycle was completed!\n");                    
                        return;
                    }
                }                  
            }
#endif
#if 1
        nid = (int) floor(treal/extdt-cid);
	if( ! ((nid>=0)&&(nid<nextdt)) )
	{
            FATAL("nr="<<nr<<",nid="<<nid<<",treal="<<treal<<",time="<<time<<",oldtreal="<<oldtreal<<",cid="<<cid);
        }
        int jid = 0;
        if (nid != 0) {
            for (i=1;i<nid+1;i++) {
                cid++;
                jid = (int) mymod(cid,nextdt);
                if (jid == 0) {
                    nr++;
                    treal -= extdt*nextdt;
                    if (nr%printfreq == 0)
                        INFO_Printf("%d cycles were completed! treal=%e oldtreal=%e nid=%d cid=%d jid=%d\n",nr,treal,oldtreal,nid,cid,jid);
                    if (nr%saveacfreq == 0)
                        af.write(this);
                    if (nr == ncycle) {
                        INFO_Printf("final cycle was completed!\n");                    
                        return;
                    }
                }                  
                enposwext[jid] += enpos[0];
            }
            cid = (int) mymod(cid,nextdt);
        }
#endif
    }

//    INFO_Printf("iid=%d, eid=%d, ncorrstep=%d\n",iid,eid,ncorrstep);
//    fclose(fp_enposxyz);
//    fclose(fp_treal);
    
}

void YSZFrame::kmcrun_continued() {
    int i, j, k, vac_target, nei, selectnum, si, sj, newpos, mvac;
    double ExtV1, ExtV2, InteractE1, InteractE2, prob, totrate, num_rand1, num_rand2, time;
    double dx, dy, dz;

    FILE *fp_singledata;
    FILE *fp_vaclist, *fp_vaclist0;
    FILE *fp_rx, *fp_ry, *fp_rz;
    FILE *fp_corrs, *fp_corrtimes, *fp_corrdxs;
    FILE *fp_corrl, *fp_corrtimel, *fp_corrdxl;
    FILE *fp_treal;
    FILE *fp_histEb0;

    /* allocate memory */
    Realloc(corrs,double,ncorrsteps);
    Realloc(corrtimes,double,ncorrbins);
//    Realloc(corrdxs,double,num_vac*ncorrbins);
    Realloc(corrdxs,double,ncorrbins);
    Realloc(corrl,double,ncorrstepl);
    Realloc(corrtimel,double,ncorrbinl);
//    Realloc(corrdxl,double,num_vac*ncorrbinl);
    Realloc(corrdxl,double,ncorrbinl);
    Realloc(histEb0,double,3);
    Realloc(histEb,double,1000);
    vaclist0 = (int *) calloc(num_vac,sizeof(int));
    memset(histEb0,0,sizeof(double)*3);
    memset(histEb,0,sizeof(double)*1000);

    rx = (int *) calloc(num_vac,sizeof(int));
    ry = (int *) calloc(num_vac,sizeof(int));
    rz = (int *) calloc(num_vac,sizeof(int));
    
/*  Loading the data from the previous simulation */
    char string[500];
    double tmpdata;
    int tmpdata2;
    int k_c, curstep_c, iids_c, eids_c, iidl_c, eidl_c;
    
    fp_singledata = fopen("singledata.dat","r");
    fp_corrs = fopen("corrs.dat","r");
    fp_corrtimes = fopen("corrtimes.dat","r");
    fp_corrdxs = fopen("corrdxs.dat","r");
    fp_corrl = fopen("corrl.dat","r");
    fp_corrtimel = fopen("corrtimel.dat","r");
    fp_corrdxl = fopen("corrdxl.dat","r");
    fp_rx = fopen("rx.dat","r");
    fp_ry = fopen("ry.dat","r");
    fp_rz = fopen("rz.dat","r");
    fp_vaclist = fopen("vaclist.dat","r");
    fp_vaclist0 = fopen("vaclist0.dat","r");
    
    fgets(string,500,fp_singledata);
    sscanf(string,"%d %d %d %d %d %d %le",&k_c,&curstep_c,&iids_c,&eids_c,&iidl_c,&eidl_c,&treal);
    
    for (i=0;i<ncorrsteps;i++) {
        fgets(string,500,fp_corrs);
        sscanf(string,"%le",&tmpdata);
	corrs[i] = tmpdata;
    }
    for (i=0;i<ncorrbins;i++) {
        fgets(string,500,fp_corrtimes);
        sscanf(string,"%le",&tmpdata);
	corrtimes[i] = tmpdata;
        fgets(string,500,fp_corrdxs);
        sscanf(string,"%le",&tmpdata);
	corrdxs[i] = tmpdata;       
    }
    for (i=0;i<ncorrstepl;i++) {
        fgets(string,500,fp_corrl);
        sscanf(string,"%le",&tmpdata);
	corrl[i] = tmpdata;
    }
    for (i=0;i<ncorrbinl;i++) {
        fgets(string,500,fp_corrtimel);
        sscanf(string,"%le",&tmpdata);
	corrtimel[i] = tmpdata;
        fgets(string,500,fp_corrdxl);
        sscanf(string,"%le",&tmpdata);
	corrdxl[i] = tmpdata;       
    }
    for (i=0;i<num_vac;i++) {
        fgets(string,500,fp_rx);
        sscanf(string,"%d",&tmpdata2);
        rx[i] = tmpdata2;
        fgets(string,500,fp_ry);
        sscanf(string,"%d",&tmpdata2);
        ry[i] = tmpdata2;
        fgets(string,500,fp_rz);
        sscanf(string,"%d",&tmpdata2);
        rz[i] = tmpdata2;
        fgets(string,500,fp_vaclist);
        sscanf(string,"%d",&tmpdata2);
        vaclist[i] = tmpdata2;
        fgets(string,500,fp_vaclist0);
        sscanf(string,"%d",&tmpdata2);
        vaclist0[i] = tmpdata2;
    }
    fclose(fp_singledata);
    fclose(fp_corrs);
    fclose(fp_corrtimes);
    fclose(fp_corrdxs);
    fclose(fp_corrl);
    fclose(fp_corrtimel);
    fclose(fp_corrdxl);
    fclose(fp_rx);
    fclose(fp_ry);
    fclose(fp_rz);
    fclose(fp_vaclist);
    fclose(fp_vaclist0);

    /* some test whether the load was succesful or not */
//    FILE *fp_testload;
//    fp_testload = fopen("testload.dat","w");
//    fprintf(fp_testload,"iid=%d treal=%e\n",iid,treal);
//    for (i=0;i<10;i++)
//        fprintf(fp_testload,"corr[%d]=%e\n",i,corr[i]);
//    for (i=0;i<10;i++)
//        fprintf(fp_testload,"corrtime[%d]=%e\n",i,corrtime[i]);
//    for (i=0;i<10;i++)
//        fprintf(fp_testload,"corrdx[%d]=%e\n",i,corrdx[i]);
//    for (i=0;i<10;i++) 
//        fprintf(fp_testload,"rx[%d]=%d\n",i,rx[i]);
//    for (i=0;i<10;i++) 
//        fprintf(fp_testload,"vaclist[%d]=%d\n",i,vaclist[i]);
//    fclose(fp_testload);
    
    /* inherit species and linklist about vacancy from the previous calculation */
    for (i=0;i<_NP;i++)
        if (linklist[i*LINKLEN+0] == 1)
            species[i] = 1;
    for (i=0;i<num_vac;i++)
        species[vaclist[i]] = 3;
    
    
    /* transfer from real position(double) to lattice position(int) */
    for (i=0;i<_NP;i++) {
        *(_RI + 3*i + 0) = (int) round((_R[i].x/latticeconst[0]/latticesize[0][3]+0.5)*NCellX);
        *(_RI + 3*i + 1) = (int) round((_R[i].y/latticeconst[1]/latticesize[1][3]+0.5)*NCellY);
        *(_RI + 3*i + 2) = (int) round((_R[i].z/latticeconst[2]/latticesize[2][3]+0.5)*NCellZ);
    }

    /* make CoulombE bet. vacancies be zero*/
    if (CTL_CoulombE == 0) {
        for (i=0;i<NCellX*NCellY*NCellZ;i++)
	    *(CoulombE + i) = 0.0;
	INFO_Printf("Coulomb interaction bet V-V is turned off\n");
    }

    /* make Energy barrier be same */
    if (CTL_Ebarrier == 0) {
        for (i=0;i<_NP*NUMJUMPDIR;i++)
	    *(Ebarrier + i) = Ebconst;
	INFO_Printf("Energy barrier is set to same numbers;Eb=%e\n",Ebconst);
    }

    /* make Coulomb Potential from Cations be zero*/
    if (CTL_CatPot == 0) {
        for (i=0;i<_NP;i++)
	    *(CatPot + i) = 0.0;
	INFO_Printf("Coulomb interaction bet V-Cat is OFF\n");
    }

//    INFO_Printf("Set EnnY=%e, EsnY=%e\n",EnnY,EsnY);
    
    double kT = KB * _T;
    int unitmove = 2;       /* each vacancy moves by unitmove(int), 2 */
    double unitensmove = (double) unitmove / num_vac; /* average unitmove in view of all vacancies */
    int problength = NUMJUMPDIR*num_vac;

    double rate[problength], cumrate[problength], cumprob[problength];
    double PDInteractE[num_vac][NUMJUMPDIR];

    enpos[0]=enpos[1]=enpos[2]=0;
    for (i=0;i<num_vac;i++) {
        enpos[0] += rx[i];
        enpos[1] += ry[i];
        enpos[2] += rz[i];
    }
    enpos[0] = enpos[0]/num_vac;
    enpos[1] = enpos[1]/num_vac;
    enpos[2] = enpos[2]/num_vac;
    
    /* Initial probability calculation: At each kmc step, only probability difference
       is passed to the next kmc step and updated by multiplying the probability
       difference to the previous probability  */
    INFO_Printf("Continued Initial probability calculation\n");
    for (i=0;i<num_vac;i++) {
        vac_target = vaclist[i];
        vaclist0[i] = vac_target;
        InteractE1 = CalInteract(vac_target);
        for (j=0;j<NUMJUMPDIR;j++) {
            nei = linklist[vac_target*LINKLEN+j+5];
//            printf("nei for %dth vac(%d) in %d direction is %d\n",i,vac_target,j,nei);
            if (species[nei] == 1) {
                vaclist[i] = nei;
                InteractE2 = CalInteract(nei);
                PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT);
//                PDInteractE[i][j] = 1.0;;
                vaclist[i] = vac_target;
            }
            else if (species[nei] == 3) {
                PDInteractE[i][j] = 0.0;
            }
            else {
                INFO_Printf("Wrong Species Error(1)! \n");
                return;
            }
//            if (PDInteractE[i][j] == 0.0)
//                printf("Initial zero:i=%d,j=%d,species=%d,%d(%d,%d,%d)->%d(%d,%d,%d)\n",i,j,species[nei],vac_target,_RI[3*vac_target],_RI[3*vac_target+1],_RI[3*vac_target+2],nei,_RI[3*nei],_RI[3*nei+1],_RI[3*nei+2]);
        }
    }
    

    /* start MC steps */
    iids = iids_c; eids = eids_c;  /*iid and eid is recoverd*/
    iidl = iidl_c; eidl = eidl_c;  /*iid and eid is recoverd*/
    double E_eci1, E_eci2;
    if (curstep_c == num_jump-1) {
        k_c++;
        curstep_c = 0;
    }
    else
        curstep_c++;
    
    for (k=k_c;k<num_loop;k++) {
        for (curstep=curstep_c;curstep<num_jump;curstep++) {
            /* update the probability rate for current kmc step */
            for (i=0;i<num_vac;i++) {
                vac_target = vaclist[i];
                ExtV1 = 0;
                E_eci1 = cal_energy_by_CEM();
                for (j=0;j<NUMJUMPDIR;j++) {
                    nei = linklist[vac_target*LINKLEN+j+5];
                    ExtV2 = 0; 
                    if (species[nei] == 1) {
                        prob = trialfreq*exp(-(Ebarrier[vac_target*NUMJUMPDIR+j]+0.5*(ExtV2-ExtV1))/kT)*exp(-(CatPot[nei]-CatPot[vac_target])/2/kT)*PDInteractE[i][j];
                        prob *= exp(-EnnY*(nnYlist[nei]-nnYlist[vac_target])/2/kT)*exp(-EsnY*(snYlist[nei]-snYlist[vac_target])/2/kT);
                        /*** cluster_expansion_method correction ***/
                        vaclist[i] = nei;
                        species[nei] = 3;
                        species[vac_target] = 1;
                        E_eci2 = cal_energy_by_CEM();
                        vaclist[i] = vac_target;
                        species[nei] = 1;
                        species[vac_target] = 3;
                        /*******************************************/
//                    printf("PDInteractE[%d][%d]=%e\n",i,j,PDInteractE[i][j]);
//                    prob = trialfreq*exp(-Ebarrier[vac_target*NUMJUMPDIR+j]/kT);
//                    if ( PDInteractE[i][j] == 0 ) {
//                        printf("curstep=%d, prob[%d][%d]=0, vac_target=%d(%d,%d,%d),nei=%d(%d,%d,%d),species=%d\n",curstep,i,j,vac_target,_RI[3*vac_target],_RI[3*vac_target+1],_RI[3*vac_target+2],nei,_RI[3*nei],_RI[3*nei+1],_RI[3*nei+2],species[nei]);
//                        return;
//                    }
                    }
                    else if (species[nei] == 3)
                        prob = 0.0;
                    else {
                        INFO_Printf("Wrong Species Error(2)! \n");
                        return;
                    }
                    rate[j+6*i] = prob;
                }
            }

            /* generate cumulated probability list */
            cumrate[0] = rate[0];
            for (i=1;i<problength;i++) {
                cumrate[i] = cumrate[i-1] + rate[i];
            }
            totrate = cumrate[problength-1];
            for (i=0;i<problength;i++) 
                cumprob[i] = cumrate[i] / totrate;
            
            /* event selection */
            num_rand1 = drand48();
            time = -1*log(num_rand1) / totrate;
            num_rand2 = drand48();
            selectnum = -1;
            for (i=0;i<problength;i++)
                if (cumprob[i] > num_rand2) {
                    selectnum = i;
                    break;
                }
            si = (int) floor(selectnum/NUMJUMPDIR);  /* determine vacancy ID */
            sj = selectnum - NUMJUMPDIR*si;          /* determine moving direction */
            treal += time;                           /* update the real time */
            
//        if((curstep%printfreq)==0)
//            INFO_Printf("%dth kmcstep: vac:%d, dir:%d \n",curstep+1,si,sj);
            
            /* update vacancy list, vaclist, and update the ensemble average position
               of total vacancy */
            vac_target = vaclist[si];
            newpos = linklist[vac_target*LINKLEN+sj+5];
            if (species[vac_target] == 1) {
                INFO_Printf("Wrong Species Error(3)!\n");
                return;
            }
            species[vac_target] = 1;
            if (species[newpos] == 3) {
                INFO_Printf("Wrong Species Error(4)!\n");
                return;
            }
            species[newpos] = 3;
            vaclist[si] = newpos;
//        printf("curstep=%d, moved #%d:(%d,%d,%d) -> #%d:(%d,%d,%d)\n",curstep,vac_target,_RI[3*vac_target],_RI[3*vac_target+1],_RI[3*vac_target+2],newpos,_RI[3*newpos],_RI[3*newpos+1],_RI[3*newpos+2]);
            switch (sj) {
            case 0: {rx[si] += unitmove; dx =  unitensmove; dy = 0.0; dz = 0.0; break;}
            case 1: {rx[si] -= unitmove; dx = -unitensmove; dy = 0.0; dz = 0.0; break;}
            case 2: {ry[si] += unitmove; dx = 0.0; dy =  unitensmove; dz = 0.0; break;}
            case 3: {ry[si] -= unitmove; dx = 0.0; dy = -unitensmove; dz = 0.0; break;}
            case 4: {rz[si] += unitmove; dx = 0.0; dy = 0.0; dz =  unitensmove; break;}
            case 5: {rz[si] -= unitmove; dx = 0.0; dy = 0.0; dz = -unitensmove; break;}
            default: dx = dy = dz = 0.0; FATAL("unrecoganized sj"<<sj);
            }
            enpos[0] += dx;
            enpos[1] += dy;
            enpos[2] += dz;
            CorrFcnX(dx);
//        CorrFcnXn(si,sj);
#if 1
            /* update the probability difference */
            for (i=0;i<num_vac;i++) {
                if (i == si) {
                    InteractE1 = CalInteract(newpos);
                    for (j=0;j<NUMJUMPDIR;j++) {
                        nei = linklist[newpos*LINKLEN+5+j];
                        if (species[nei] == 1) {
                            vaclist[i] = nei;
                            InteractE2 = CalInteract(nei);
                            PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                            vaclist[i] = newpos;
                        }
                        else if(species[nei] == 3)
                            PDInteractE[i][j] = 0;                
                        else {
                            INFO_Printf("Wrong Species Error(5)!\n");
                            return;
                        }
                    }
                }
                else {
                    mvac = vaclist[i];
                    for (j=0;j<NUMJUMPDIR;j++) {
                        nei = linklist[mvac*LINKLEN+5+j];
                        if (nei == vac_target) {
                            InteractE1 = CalInteract(mvac);
                            vaclist[i] = nei;
                            InteractE2 = CalInteract(nei);
                            PDInteractE[i][j] = exp(-(InteractE2-InteractE1)/2/kT); 
                            vaclist[i] = mvac;
                        }
                        else                        
                            PDInteractE[i][j] *= PDiffInteract(mvac,j,vac_target,newpos);
                    }
                }
            }
#endif
//        fprintf(fp_enposxyz,"%6d %6.6e %6.6e %6.6e\n",curstep+1,enpos[0],enpos[1],enpos[2]);
//        for (i=0;i<num_vac;i++)
//            fprintf(fp_posx,"%d ",rx[i]);
//        fprintf(fp_posx,"\n");       
//        fprintf(fp_treal,"%6.10e\n",treal);


            if(savecorr)
            {
                if(((curstep+1)%savecorrfreq)==0)
                {
                    cfs.write(this);
                    cfl.write(this);
                }
            }   
            
            if(savecn)
            {
                if(((curstep+1)%savecnfreq)==0&&(curstep!=0))
                    intercn.write(this);            
            }
            
            if(savediff)
            {
                if(((curstep+1)%savedifffreq)==0) {
                    df.write(this);
                }
            }

            if(savehist)
            {
                if(((curstep+1)%savehistfreq)==0) {
                    hf.write(this);
                }
                if(((curstep+1)%savehistfreq)==0) {
                    fp_histEb0 = fopen("histEb0.dat","a");
                    fprintf(fp_histEb0,"%e %e %e\n",histEb0[0],histEb0[1],histEb0[2]);
                }
            }

            if(savecontdata)
            {
                if (((curstep+1)%savecontfreq)==0) {
                    fp_singledata = fopen("singledata.dat","w");
                    fp_corrs = fopen("corrs.dat","w");
                    fp_corrtimes = fopen("corrtimes.dat","w");
                    fp_corrdxs = fopen("corrdxs.dat","w");
                    fp_corrl = fopen("corrl.dat","w");
                    fp_corrtimel = fopen("corrtimel.dat","w");
                    fp_corrdxl = fopen("corrdxl.dat","w");
                    fp_rx = fopen("rx.dat","w");
                    fp_ry = fopen("ry.dat","w");
                    fp_rz = fopen("rz.dat","w");
                    fp_vaclist = fopen("vaclist.dat","w");
                    fp_vaclist0 = fopen("vaclist0.dat","w");
                    
                    fprintf(fp_singledata,"%d %d %d %d %d %d %e\n",k,curstep,iids,eids,iidl,eidl,treal);
                    for (i=0;i<ncorrsteps;i++) {
                        fprintf(fp_corrs,"%e\n",corrs[i]);
                    }
                    for (i=0;i<ncorrbins;i++) {
                        fprintf(fp_corrtimes,"%e\n",corrtimes[i]);
                        fprintf(fp_corrdxs,"%e\n",corrdxs[i]);       
                    }
                    for (i=0;i<ncorrstepl;i++) {
                        fprintf(fp_corrl,"%e\n",corrl[i]);
                    }
                    for (i=0;i<ncorrbinl;i++) {
                        fprintf(fp_corrtimel,"%e\n",corrtimel[i]);
                        fprintf(fp_corrdxl,"%e\n",corrdxl[i]);       
                    }
                    for (i=0;i<num_vac;i++) {
                        fprintf(fp_rx,"%d\n",rx[i]);
                        fprintf(fp_ry,"%d\n",ry[i]);
                        fprintf(fp_rz,"%d\n",rz[i]);
                        fprintf(fp_vaclist,"%d\n",vaclist[i]);
                        fprintf(fp_vaclist0,"%d\n",vaclist0[i]);
                    }
                    fclose(fp_singledata);
                    fclose(fp_corrs);
                    fclose(fp_corrtimes);
                    fclose(fp_corrdxs);
                    fclose(fp_corrl);
                    fclose(fp_corrtimel);
                    fclose(fp_corrdxl);
                    fclose(fp_rx);
                    fclose(fp_ry);
                    fclose(fp_rz);
                    fclose(fp_vaclist);
                    fclose(fp_vaclist0);
                    INFO_Printf("Saved all data\n");
                }
            }
        }
    }

//    INFO_Printf("iid=%d, eid=%d, ncorrstep=%d\n",iid,eid,ncorrstep);
//    fclose(fp_enposxyz);
//    fclose(fp_posx);
    
}

int YSZFrame::Barrier(int vac_target, int nei, int j) {
  return 0;
}

double YSZFrame::CalInteract(int Annum) {
#if 0    
    int k,Vacnum;
    int index1, index2;
    double InteractE;
    int rx, ry, rz;

    InteractE = 0;
    index1 = Annum*3;

    for (k=0;k<num_vac;k++) {
        Vacnum = vaclist[k];
        if (Annum != Vacnum) {
            index2 = Vacnum*3;
            rx = (int) mymod(*(_RI+index2+0) - *(_RI+index1+0) + NCellX/2,NCellX);
            ry = (int) mymod(*(_RI+index2+1) - *(_RI+index1+1) + NCellY/2,NCellY);
            rz = (int) mymod(*(_RI+index2+2) - *(_RI+index1+2) + NCellZ/2,NCellZ);
            InteractE += *(CoulombE + NCellY*NCellZ*rx + NCellZ*ry + rz);
        }
    }
    if (InteractE != InteractE) {
        printf("Total Interaction is NaN\n");
        abort();
    }
    return(InteractE);
#endif
    return 0.0;
}

double YSZFrame::PDiffInteract(int i2, int j2, int old1, int new1) {
#if 1    
    int rx, ry, rz;
    int index1, index2, index3, index4;
    double PDInteractE,/* disto2o1, disto2n1, distn2o1, distn2n1,*/ VIo2o1, VIo2n1, VIn2o1, VIn2n1;
    
    index1 = i2*3;
    index2 = linklist[i2*LINKLEN+5+j2]*3;
    index3 = old1*3;
    index4 = new1*3;
    
    // dist to o2 in view of o1
    rx = (int) mymod(*(_RI+index1+0) - *(_RI+index3+0) + NCellX/2,NCellX);
    ry = (int) mymod(*(_RI+index1+1) - *(_RI+index3+1) + NCellY/2,NCellY);
    rz = (int) mymod(*(_RI+index1+2) - *(_RI+index3+2) + NCellZ/2,NCellZ);
    VIo2o1 = *(CoulombE + NCellY*NCellZ*rx + NCellZ*ry + rz); 
    
    // dist to o2 in view of n1
    rx = (int) mymod(*(_RI+index1+0) - *(_RI+index4+0) + NCellX/2,NCellX);
    ry = (int) mymod(*(_RI+index1+1) - *(_RI+index4+1) + NCellY/2,NCellY);
    rz = (int) mymod(*(_RI+index1+2) - *(_RI+index4+2) + NCellZ/2,NCellZ);
    VIo2n1 = *(CoulombE + NCellY*NCellZ*rx + NCellZ*ry + rz);
    
    // dist to n2 in view of o1
    rx = (int) mymod(*(_RI+index2+0) - *(_RI+index3+0) + NCellX/2,NCellX);
    ry = (int) mymod(*(_RI+index2+1) - *(_RI+index3+1) + NCellY/2,NCellY);
    rz = (int) mymod(*(_RI+index2+2) - *(_RI+index3+2) + NCellZ/2,NCellZ);
    VIn2o1 = *(CoulombE + NCellY*NCellZ*rx + NCellZ*ry + rz);
    if (rx == NCellX/2 && ry == NCellY/2 && rz == NCellZ/2) // maintain same prob
        return(1);
    
    // dist to n2 in view of n1
    rx = (int) mymod(*(_RI+index2+0) - *(_RI+index4+0) + NCellX/2,NCellX);
    ry = (int) mymod(*(_RI+index2+1) - *(_RI+index4+1) + NCellY/2,NCellY);
    rz = (int) mymod(*(_RI+index2+2) - *(_RI+index4+2) + NCellZ/2,NCellZ);
    VIn2n1 = *(CoulombE + NCellY*NCellZ*rx + NCellZ*ry + rz);
    if (rx == NCellX/2 && ry == NCellY/2 && rz == NCellZ/2) { // actually infinite energy
//        printf("curstep=%d,prob0!:%d->%d when %d->%d, distand is 0 between n1:%d(%d,%d,%d) and n2:%d(%d,%d,%d)\n",curstep,i2,linklist[i2*LINKLEN+5+j2],old1,new1,new1,*(_RI+index4+0),*(_RI+index4+1),*(_RI+index4+2),linklist[i2*LINKLEN+5+j2],*(_RI+index2+0),*(_RI+index2+1),*(_RI+index2+2));
        return(0);
    }
    
    PDInteractE = exp(-(VIn2n1 - VIo2n1 - VIn2o1 + VIo2o1)/2/KB/_T);
    
    return(PDInteractE);
#endif
  return 0.0;
}

void YSZFrame::FindMenergy() {
#if 0    
    int i, j, nei;
    int arraynum;
    for (i=0;i<total_anion;i++) {
        for (j=0;j<6;j++) {
            nei = *(YSZAn+i*18+j+7);
        }
    }
#endif
}

double YSZFrame::mymod(int x, int y) {
#if 1
  double temp = x - y*floor((double)x/y);
  return(temp);
#endif
}

double YSZFrame::mymaxabs(int x, int y, int z) {

    int tmp_max;

    tmp_max = abs(x);
    if (abs(y)>tmp_max)
        tmp_max = abs(y);
    if (abs(z)>tmp_max)
        tmp_max = abs(z);

    return 1.0*tmp_max;

}

void YSZFrame::mysort(int n, int *x) {

    int i, j, tmp;
    for (i=1;i<n;i++)
        for (j=0;j<i;j++)
            if (*(x+i) < *(x+j)) {
                tmp = *(x+i);
                *(x+i) = *(x+j);
                *(x+j) = tmp;
                break;
            }

}
/**************************************************************/
/* Calculate Correlation function*/
/**************************************************************/

void YSZFrame::CorrFcnX(double dx)
{
#if 1
    int i, ii, eeid, id, istep, iistep;
    double maxid;
// need publc arrays: corrtime, corr, corrdx
    
    eeid = (int) curstep%ncorrbins;
    if (eeid == 0 && curstep != 0) {
        maxid = floor((treal - corrtimes[0])/corrdts);
        if (maxid < ncorrsteps) {
            ncorrbins*=2;
            Realloc(corrtimes,double,ncorrbins);
            Realloc(corrdxs,double,ncorrbins);
            WARNING("Corr_short info bin Pointers are resized t="<<treal<<",corrtimes[0]="<<corrtimes[0]<<",maxid="<<maxid);
            eeid = (int) curstep%ncorrbins;
        }
//        INFO_Printf("eid(=curstep=%d) is set to zero.\n",curstep);
    }
    corrtimes[eeid] = treal;
    corrdxs[eeid] = dx;
//    INFO_Printf("curstep=%d,iid=%d,eid=%d,eeid=%d,corrdx[%d]=dx=%10.6e, treal=%10.6e\n",
//                curstep,iid,eid,eeid,eeid,corrdx[eeid],treal);
    istep = (int) mymod(iids,ncorrbins);
    iistep = (int) mymod(curstep-iids,ncorrbins);
    if (istep < 0)
      FATAL("Negative istep!, iids="<<iids<<",istep="<<istep);
    for (i=istep;i<=istep+iistep;i++) {
        ii = (int) i%ncorrbins;
        id = (int) floor((treal - corrtimes[ii])/corrdts);
        if ( id < 0)
	  FATAL("Negative correlation bin,,iids="<<iids<<",ii="<<ii<<",treal="<<treal<<",corrtimes="<<corrtimes<<",id="<<id);
        if (id<ncorrsteps) {
            corrs[id]+= corrdxs[eeid]*corrdxs[ii];
//            INFO_Printf("dtime=%6.6e, corr[%5d]=%e,i=%d,id=%d,ii=%d, di*dj=%e*%e,iid=%d,eid=%d\n",
//                        treal-corrtime[ii],id,corr[id],i,id,ii,corrdx[eeid],corrdx[ii],iid,eid);
        }
        else {
//            INFO_Printf("iid++, id >= ncorrstep. time difference is %6.6e\n",treal-corrtime[ii]);
            iids++;
        }
    }
    
/***********************************************************/
    
    eeid = (int) curstep%ncorrbinl;
    if (eeid == 0 && curstep != 0) {
        maxid = floor((treal - corrtimel[0])/corrdtl);
        if (maxid < ncorrstepl) {
            ncorrbinl*=2;
            Realloc(corrtimel,double,ncorrbinl);
            Realloc(corrdxl,double,ncorrbinl);
            WARNING("Corr_long info bin Pointers are resized t="<<treal<<",corrtimel[0]="<<corrtimel[0]<<",maxid="<<maxid);
            eeid = (int) curstep%ncorrbinl;
        }
//        INFO_Printf("eid(=curstep=%d) is set to zero.\n",curstep);
    }
    corrtimel[eeid] = treal;
    corrdxl[eeid] = dx;
//    INFO_Printf("curstep=%d,iid=%d,eid=%d,eeid=%d,corrdx[%d]=dx=%10.6e, treal=%10.6e\n",
//                curstep,iid,eid,eeid,eeid,corrdx[eeid],treal);
    istep = (int) mymod(iidl,ncorrbinl);
    iistep = (int) mymod(curstep-iidl,ncorrbinl);
    if (istep < 0)
      FATAL("Negative istep!, iidl="<<iidl<<",istep="<<istep);
    for (i=istep;i<=istep+iistep;i++) {
        ii = (int) i%ncorrbinl;
        id = (int) floor((treal - corrtimel[ii])/corrdtl);
        if ( id < 0)
	  FATAL("Negative correlation bin,,iidl="<<iidl<<",ii="<<ii<<",treal="<<treal<<",corrtimel="<<corrtimel<<",id="<<id);
        if (id<ncorrstepl) {
            corrl[id]+= corrdxl[eeid]*corrdxl[ii];
//            INFO_Printf("dtime=%6.6e, corr[%5d]=%e,i=%d,id=%d,ii=%d, di*dj=%e*%e,iid=%d,eid=%d\n",
//                        treal-corrtime[ii],id,corr[id],i,id,ii,corrdx[eeid],corrdx[ii],iid,eid);
        }
        else {
//            INFO_Printf("iid++, id >= ncorrstep. time difference is %6.6e\n",treal-corrtime[ii]);
            iidl++;
        }
    }
#endif
}

void YSZFrame::CorrFcnXn(int xi, int xj)
{
#if 0
    int i, ii, eeid, id, istep;
    double dx;
// need publc arrays: corrtime, corr, corrdx 
    if (xj == 0)
      dx = 2.0;
    else if (xj == 1)
      dx = -2.0;
    else
      dx = 0.0;

    eeid = (int) curstep%ncorrbin;
    if (eeid == 0 && curstep != 0) {     
        id = (int) floor((treal - corrtime[0])/corrdt);
        if (id < ncorrstep) {
            ncorrbin++;
            Realloc(corrtime,double,ncorrbin);
            Realloc(corrdx,double,ncorrbin);
            WARNING("Be Careful! Corr info bin Pointers are resized\n");
	    eeid = (int) curstep%ncorrbin;
        }
    }
    corrtime[eeid] = treal;
    for (i=0;i<num_vac;i++)
      corrdx[ncorrbin*i+eeid] = 0.0;
    corrdx[ncorrbin*xi+eeid] = dx;

    istep = iid;
    for (i=istep;i<=curstep;i++) {
        ii = (int) i%ncorrbin;
        id = (int) floor((treal - corrtime[ii])/corrdt);
        if ( id < 0)
            FATAL("Negative correlation bin, id="<<id);
        if (id<ncorrstep) {
            corr[id]+= dx*corrdx[ncorrbin*xi+ii]/num_vac;
        }
        else {
            iid++;
        }
    }
#endif
}

/*********************************************************************/
/* Cluster Expansion Method version 2 - generate clusters            */
/*********************************************************************/
int YSZFrame::ReadClusters()
{
    FILE *fp;
    char string[100];
    char nameforload1[400];
    char nameforload2[400];
    int i, j;

    Realloc(cCC, int,4*nCC);    
    Realloc(cAA, int,4*nAA);
    Realloc(cCA, int,4*nCA);
    Realloc(cCCC,int,3*nCCC);
    Realloc(cCCA,int,3*nCCA);
    Realloc(cCAA,int,3*nCAA);
    Realloc(cAAA,int,3*nAAA);
    Realloc(eci_opt, double,neci);
    memset(cCC, 0,sizeof(int)*4*nCC);
    memset(cAA, 0,sizeof(int)*4*nAA);
    memset(cCA, 0,sizeof(int)*4*nCA);
    memset(cCCC,0,sizeof(int)*3*nCCC);
    memset(cCCA,0,sizeof(int)*3*nCCA);
    memset(cCAA,0,sizeof(int)*3*nCAA);
    memset(cAAA,0,sizeof(int)*3*nAAA);
    memset(eci_opt,0,sizeof(double)*neci);

    Realloc(cluster_spin,int,_NP);
    Realloc(cluster_phi,int,neci);
    Realloc(n_cluster_phi,int,neci);
    Realloc(map_An_clusters,int,NAnion*nonzeroeci);
    memset(cluster_spin,0,sizeof(int)*_NP);
    memset(cluster_phi,0,sizeof(int)*neci);
    memset(n_cluster_phi,0,sizeof(int)*neci);
    memset(map_An_clusters,0,sizeof(int)*NAnion*nonzeroeci);

    Realloc(ebcorrtable, double,nebcorrtable);
    Realloc(ebxcorrtable,double,nebcorrtable);    
    memset(ebcorrtable,0,sizeof(double)*nebcorrtable);
    memset(ebxcorrtable,0,sizeof(double)*nebcorrtable);

    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterCC_file);
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadClusters file, cCC not found!");
    }
    for (i=0;i<nCC;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d %d\n",cCC+4*i+0,cCC+4*i+1,cCC+4*i+2,cCC+4*i+3);
    }
    fclose(fp);

    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterAA_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadClusters file, cAA not found!");
    }
    for (i=0;i<nAA;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d %d\n",cAA+4*i+0,cAA+4*i+1,cAA+4*i+2,cAA+4*i+3);
    }
    fclose(fp);
    
    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterCA_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadClusters file, cCA not found!");
    }
    for (i=0;i<nCA;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d %d\n",cCA+4*i+0,cCA+4*i+1,cCA+4*i+2,cCA+4*i+3);
    }
    fclose(fp);
    
    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterCCC_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadClusters file, cCCC not found!");
    }
    for (i=0;i<nCCC;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d\n",cCCC+3*i+0,cCCC+3*i+1,cCCC+3*i+2);
    }
    fclose(fp);
    
    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterCCA_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadCluster file, cCCA not found!");
    } 
    for (i=0;i<nCCA;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d\n",cCCA+3*i+0,cCCA+3*i+1,cCCA+3*i+2);
    }
    fclose(fp);
    
    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterCAA_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadClusters file, cCAA not found!");
    } 
    for (i=0;i<nCAA;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d\n",cCAA+3*i+0,cCAA+3*i+1,cCAA+3*i+2);
    }
    fclose(fp);

    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,clusterAAA_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ReadClusters file, cAAA not found!");
    } 
    for (i=0;i<nAAA;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%d %d %d\n",cAAA+3*i+0,cAAA+3*i+1,cAAA+3*i+2);
    }
    fclose(fp);
    
    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,eci_opt_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::eci optimized solution file, eci_opt.dat not found!");
    }
    for (i=0;i<neci;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%lf",eci_opt+i);
    }
    fclose(fp);
//    for (i=0;i<4;i++)
//        eci_opt[i] = 0.0;
//    for (i=111;i<neci;i++)
//        eci_opt[i] = 0.0;
    nonzeroeci = 0;
    for (i=0;i<neci;i++)
        if ( fabs(eci_opt[i] ) > 0.0)
            nonzeroeci++;
    INFO_Printf("The number of nonzeroeci = %d\n",nonzeroeci);
    Realloc(clusters_to_use,int,nonzeroeci);
    memset(clusters_to_use,0,sizeof(int)*nonzeroeci);
    j = 0;
    for (i=0;i<neci;i++)
        if ( fabs(eci_opt[i] ) > 0.0) {
            *(clusters_to_use+j) = i;
            j++;
            printf("cluster to use = %d, eci=%e\n",i,eci_opt[i]);
        }
    INFO_Printf("Read clusters is done!\n");
    
    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,ebcorrtable_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ebcorrtable.dat not found!");
    }
    for (i=0;i<nebcorrtable;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%lf",ebcorrtable+i);
    }
    fclose(fp);

    strcpy(nameforload1,cemfilesdir);    
    strcpy(nameforload2,ebxcorrtable_file);    
    fp=fopen(strcat(nameforload1,nameforload2),"r");
    if(fp==NULL)
    {
        FATAL("YSZFrame::ebxcorrtable.dat not found!");
    }
    for (i=0;i<nebcorrtable;i++)
    {
        fgets(string,100,fp);
        sscanf(string,"%lf",ebxcorrtable+i);
    }
    fclose(fp);

    return 0;
}

/**************************************************************************/
/* Cluster Expansion Method version 2 - initially calculate the energy */
/**************************************************************************/
double YSZFrame::cal_energy_by_CEM_initial()
{
    int i, j, k, ii, jj, kk, l, n;
    int Nmap;
    int tmp_map,tmp1,tmp2,tmp3,tmp4;
    int tmp_dat[3];
    int nY=0;
    
    
    //obtain the position info. and update the spin info.
    for (i=0;i<_NP;i++) {
        switch(species[i]) {
        case 0:
            cluster_spin[i] = 1;  //Zr
            break;
        case 1:
            cluster_spin[i] = 1;  //O
            break;
        case 2:
            cluster_spin[i] = -1;  //Y
            nY++;
            break;
        case 3:
            cluster_spin[i] = -1;  //V
            break;
        }
//        INFO_Printf("cluster_spin[%d(%d,%d,%d)]=%d\n",i,_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],cluster_spin[i]);
    }

    for (l=0;l<nonzeroeci;l++) 
        if (clusters_to_use[l] == 0) {
            cluster_phi[0] = 1;
            n_cluster_phi[0] = 1;
            break;
        }
    for (l=0;l<nonzeroeci;l++) 
        if (clusters_to_use[l] == 1) {
            cluster_phi[1] = NAnion-2*num_vac;
            n_cluster_phi[1] = NAnion;
            break;
        }
    for (l=0;l<nonzeroeci;l++) 
        if (clusters_to_use[l] == 2) {
            cluster_phi[2] = NCation-2*nY;
            n_cluster_phi[2] = NCation;
            break;
        }   
        
    Nmap = (int) NAnion*(NAnion-1)/2;
    Realloc(map_AA,int,3*Nmap+1);
    memset(map_AA,0,sizeof(int)*(3*Nmap+1));
    n = 0;
    for (ii=0;ii<NAnion-1;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion;jj++) {
            j = anionlist[jj];            
            tmp_map = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
//            INFO_Printf("An1(%d,%d,%d),An2(%d,%d,%d),mapAA=%d\n",_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],tmp_map);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3) {
                        cluster_phi[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phi[clusters_to_use[l]]++;
                        map_AA[3*n+1] = i;
                        map_AA[3*n+2] = j;
                        map_AA[3*n+3] = clusters_to_use[l];
                        n++;
                        break;
                    }
        }
    }
    map_AA[0] = n;
    Realloc(map_AA,int,3*n+1);    
    INFO_Printf("Found all AA relationship, the sizeof map_AA = %d x 3\n",n);
    
    Nmap = (int) NCation*(NCation-1)/2;
    Realloc(map_CC,int,3*Nmap+1);
    memset(map_CC,0,sizeof(int)*(3*Nmap+1));
    n = 0;
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            tmp_map = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
//            INFO_Printf("Cat1(%d,%d,%d),Cat2(%d,%d,%d),mapCC=%d\n",_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],tmp_map);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA) {
                        cluster_phi[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phi[clusters_to_use[l]]++;
                        map_CC[3*n+1] = i;
                        map_CC[3*n+2] = j;
                        map_CC[3*n+3] = clusters_to_use[l];
                        n++;
                        break;
                    }                
        }
    }
    map_CC[0] = n;
    Realloc(map_CC,int,3*n+1);
    INFO_Printf("Found all CC relationship, the sizeof map_CC = %d x 3\n",n);

    Nmap = (int) NCation*NAnion;
    Realloc(map_CA,int,3*Nmap+1);
    memset(map_CA,0,sizeof(int)*(3*Nmap+1));
    n = 0;
    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion;jj++) {
            j = anionlist[jj];
            tmp_map = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA+nCC) {
                        cluster_phi[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phi[clusters_to_use[l]]++;
                        map_CA[3*n+1] = i;
                        map_CA[3*n+2] = j;
                        map_CA[3*n+3] = clusters_to_use[l];
                        n++;
                        break;
                    }
        }
    }
    map_CA[0] = n;
    Realloc(map_CA,int,3*n+1);
    INFO_Printf("Found all CA relationship, the sizeof map_CA = %d x 3\n",n);

    Nmap = (int) NAnion*(NAnion-1)*(NAnion-2)/6;
    if (Nmap > maxsize_cluster_buffer)
	Nmap = maxsize_cluster_buffer;
    Realloc(map_AAA,int,4*Nmap+1);
    memset(map_AAA,0,sizeof(int)*(4*Nmap+1));
    n = 0;
    for (ii=0;ii<NAnion-2;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_AA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+100; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+100;
                    mysort2(tmp_dat);
                    for (l=0;l<nAAA;l++)
                        if (tmp_dat[0]==cAAA[3*l+0] && tmp_dat[1]==cAAA[3*l+1] && tmp_dat[2]==cAAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
//                INFO_Printf("An1(%d,%d,%d),An2(%d,%d,%d),An3(%d,%d,%d),AA=%d,%d,%d,cAAAtype=%d\n",_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],tmp1+100,tmp2+100,tmp3+100,tmp4);
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA) {
                            cluster_phi[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            map_AAA[4*n+1] = i;
                            map_AAA[4*n+2] = j;
                            map_AAA[4*n+3] = k;
                            map_AAA[4*n+4] = clusters_to_use[l];
                            n++;
//                            INFO_Printf("An1(%d,%d,%d),An2(%d,%d,%d),An3(%d,%d,%d),AA=%d,%d,%d,cAAAtype=%d,prod(spin)=%d\n",_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],tmp1+100,tmp2+100,tmp3+100,clusters_to_use[l],cluster_spin[i]*cluster_spin[j]*cluster_spin[k]);
                            break;
                        }                
            }
        }
    }
    map_AAA[0] = n;
    Realloc(map_AAA,int,4*n+1);
    INFO_Printf("Found all AAA relationship, the sizeof map_AAA = %d x 4\n",n);

    Nmap = (int) NCation*NAnion*(NAnion-1)/2;
    if (Nmap > maxsize_cluster_buffer)
	Nmap = maxsize_cluster_buffer;
    Realloc(map_CAA,int,4*Nmap+1);
    memset(map_CAA,0,sizeof(int)*(4*Nmap+1));
    n = 0;
    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+200; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+200;
                    mysort2(tmp_dat);
                    for (l=0;l<nCAA;l++)
                        if (tmp_dat[0]==cCAA[3*l+0] && tmp_dat[1]==cCAA[3*l+1] && tmp_dat[2]==cCAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA) {
                            cluster_phi[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            map_CAA[4*n+1] = i;
                            map_CAA[4*n+2] = j;
                            map_CAA[4*n+3] = k;
                            map_CAA[4*n+4] = clusters_to_use[l];
                            n++;
                            break;
                        }                
            }
        }
    }
    map_CAA[0] = n;
    Realloc(map_CAA,int,4*n+1);
    INFO_Printf("Found all CAA relationship, the sizeof map_CAA = %d x 4\n",n);
 
    Nmap = (int) NCation*(NCation-1)/2*NAnion;
    if (Nmap > maxsize_cluster_buffer)
	Nmap = maxsize_cluster_buffer;
    Realloc(map_CCA,int,4*Nmap+1);
    memset(map_CCA,0,sizeof(int)*(4*Nmap+1));
    n = 0;
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            for (kk=0;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_CA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1; tmp_dat[1]=tmp2+200; tmp_dat[2]=tmp3+200;
                    mysort2(tmp_dat);
                    for (l=0;l<nCCA;l++)
                        if (tmp_dat[0]==cCCA[3*l+0] && tmp_dat[1]==cCCA[3*l+1] && tmp_dat[2]==cCCA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA+nCAA) {
                            cluster_phi[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            map_CCA[4*n+1] = i;
                            map_CCA[4*n+2] = j;
                            map_CCA[4*n+3] = k;
                            map_CCA[4*n+4] = clusters_to_use[l];
                            n++;
                            break;
                        }                
            }
        }
    }
    map_CCA[0] = n;
    Realloc(map_CCA,int,4*n+1);
    INFO_Printf("Found all CCA relationship, the sizeof map_CCA = %d x 4\n",n);

    Nmap = (int) NCation*(NCation-1)*(NCation-2)/6;
    if (Nmap > maxsize_cluster_buffer)
	Nmap = maxsize_cluster_buffer;
    Realloc(map_CCC,int,4*Nmap+1);
    memset(map_CCC,0,sizeof(int)*(4*Nmap+1));
    n = 0;
    for (ii=0;ii<NCation-2;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation-1;jj++) {
            j = cationlist[jj];
            for (kk=jj+1;kk<NCation;kk++) {
                k = cationlist[kk];
                tmp1 = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_CC(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CC(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1; tmp_dat[1]=tmp2; tmp_dat[2]=tmp3;
                    mysort2(tmp_dat);
                    for (l=0;l<nCCC;l++)
                        if (tmp_dat[0]==cCCC[3*l+0] && tmp_dat[1]==cCCC[3*l+1] && tmp_dat[2]==cCCC[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA+nCAA+nCCA) {
                            cluster_phi[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            map_CCC[4*n+1] = i;
                            map_CCC[4*n+2] = j;
                            map_CCC[4*n+3] = k;
                            map_CCC[4*n+4] = clusters_to_use[l];
                            n++;
                            break;
                        }                
            }
        }
    }
    map_CCC[0] = n;
    Realloc(map_CCC,int,4*n+1);
    INFO_Printf("Found all CCC relationship, the sizeof map_CCC = %d x 4\n",n);

    /***** Build a map from Anion id to cluster *****/
    /* This map requires huge size of memory so use it carefully. */
    /* ii is stored instead of i to reduce the size of arrays   */
    Nmap = (int) NAnion*map_AA[0];
    Realloc(map_An_AA,int,Nmap);
    memset(map_An_AA,-1,sizeof(int)*Nmap);
    Nmap = (int) NAnion*map_CA[0];
    Realloc(map_An_CA,int,Nmap);
    memset(map_An_CA,-1,sizeof(int)*Nmap);
    Nmap = (int) NAnion*map_AAA[0];
    Realloc(map_An_AAA,int,Nmap);
    memset(map_An_AAA,-1,sizeof(int)*Nmap);
    Nmap = (int) NAnion*map_CAA[0];
    Realloc(map_An_CAA,int,Nmap);
    memset(map_An_CAA,-1,sizeof(int)*Nmap);
    Nmap = (int) NAnion*map_CCA[0];
    Realloc(map_An_CCA,int,Nmap);
    memset(map_An_CCA,-1,sizeof(int)*Nmap);    
    for (ii=0;ii<NAnion;ii++) {
        i = anionlist[ii];
        // find map_AA //
        for (j=0;j<map_AA[0];j++)
            if (map_AA[3*j+1]==i || map_AA[3*j+2]==i) {
                map_An_AA[map_AA[0]*ii+j] = 1;
            }
        // find map_CA //
        for (j=0;j<map_CA[0];j++)
            if (map_CA[3*j+1]==i || map_CA[3*j+2]==i) {
                map_An_CA[map_CA[0]*ii+j] = 1;
            }
        // find map_AAA //
        for (j=0;j<map_AAA[0];j++)
            if (map_AAA[4*j+1]==i || map_AAA[4*j+2]==i || map_AAA[4*j+3]==i) {
                map_An_AAA[map_AAA[0]*ii+j] = 1;
            }
        // find map_CAA //
        for (j=0;j<map_CAA[0];j++)
            if (map_CAA[4*j+1]==i || map_CAA[4*j+2]==i || map_CAA[4*j+3]==i) {
                map_An_CAA[map_CAA[0]*ii+j] = 1;
            }
        // find map_CCA //
        for (j=0;j<map_CCA[0];j++)
            if (map_CCA[4*j+1]==i || map_CCA[4*j+2]==i || map_CCA[4*j+3]==i) {
                map_An_CCA[map_CCA[0]*ii+j] = 1;
            }        
    }

//    FILE *fp;
//    fp = fopen("map_An_AA.dat","w");
//    for (ii=0;ii<NAnion;ii++) {
//        i = anionlist[ii];
//        fprintf(fp,"%d %d %d ",_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
//        for (j=0;j<map_AA[0];j++) {
//            fprintf(fp,"%d ",map_AA[3*map_An_AA[map_AA[0]*ii+j]+3]);
//        }
//        fprintf(fp,"\n");
//    }
    
    double energy_cluster = 0.0;
    for (i=0;i<nonzeroeci;i++) {
        INFO_Printf("cluster to use = %d, eci_opt = %e, phi = %d, n_phi = %d\n",clusters_to_use[i],eci_opt[clusters_to_use[i]],cluster_phi[clusters_to_use[i]],n_cluster_phi[clusters_to_use[i]]);
        energy_cluster += eci_opt[clusters_to_use[i]]*cluster_phi[clusters_to_use[i]];
    }
//    INFO_Printf("CEM Calculation is done\n");

    return(energy_cluster);
    
}

/**************************************************************************/
/* Clulster Expansion Method version 2 - calculate energy by CEM          */
/**************************************************************************/
double YSZFrame::cal_energy_by_CEM()
{
    int i, j, k, ii, jj, kk, l;
    int tmp_map,tmp1,tmp2,tmp3,tmp4;
    int tmp_dat[3];
    
    //obtain the position info. and update the spin info.
    for (i=0;i<_NP;i++)
        switch(species[i]) {
        case 0:
            cluster_spin[i] = 1;  //Zr
            break;
        case 1:
            cluster_spin[i] = 1;  //O
            break;
        case 2:
            cluster_spin[i] = -1;  //Y
            break;
        case 3:
            cluster_spin[i] = -1;  //V
            break;
        }
    for (ii=0;ii<NAnion-1;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion;jj++) {
            j = anionlist[jj];            
            tmp_map = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3) {
                        cluster_phi[clusters_to_use[l]]  = cluster_spin[i]*cluster_spin[j];
                        n_cluster_phi[clusters_to_use[l]]++;
                        break;
                    }
        }
    }    
//    INFO_Printf("Found all AA relationship\n");
    
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            tmp_map = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA) {
                        cluster_phi[clusters_to_use[l]]  = cluster_spin[i]*cluster_spin[j];
                        n_cluster_phi[clusters_to_use[l]]++;
                        break;
                    }                
        }
    }
//    INFO_Printf("Found all CC relationship\n");

    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion;jj++) {
            j = anionlist[jj];
            tmp_map = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA+nCC) {
                        cluster_phi[clusters_to_use[l]] = cluster_spin[i]*cluster_spin[j];
                        n_cluster_phi[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
//    INFO_Printf("Found all CA relationship\n");

    for (ii=0;ii<NAnion-2;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_AA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+100; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+100;
                    mysort2(tmp_dat);
                    for (l=0;l<nAAA;l++)
                        if (tmp_dat[0]==cAAA[3*l+0] && tmp_dat[1]==cAAA[3*l+1] && tmp_dat[2]==cAAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA) {
                            cluster_phi[clusters_to_use[l]] = cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
//    INFO_Printf("Found all AAA relationship\n");

    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+200; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+200;
                    mysort2(tmp_dat);
                    for (l=0;l<nCAA;l++)
                        if (tmp_dat[0]==cCAA[3*l+0] && tmp_dat[1]==cCAA[3*l+1] && tmp_dat[2]==cCAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA) {
                            cluster_phi[clusters_to_use[l]] = cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
//    INFO_Printf("Found all CAA relationship\n");
 
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            for (kk=0;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_CA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1; tmp_dat[1]=tmp2+200; tmp_dat[2]=tmp3+200;
                    mysort2(tmp_dat);
                    for (l=0;l<nCCA;l++)
                        if (tmp_dat[0]==cCCA[3*l+0] && tmp_dat[1]==cCCA[3*l+1] && tmp_dat[2]==cCCA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA+nCAA) {
                            cluster_phi[clusters_to_use[l]] = cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
//    INFO_Printf("Found all CCA relationship\n");

    for (ii=0;ii<NCation-2;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation-1;jj++) {
            j = cationlist[jj];
            for (kk=jj+1;kk<NCation;kk++) {
                k = cationlist[kk];
                tmp1 = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_CC(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CC(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1; tmp_dat[1]=tmp2; tmp_dat[2]=tmp3;
                    mysort2(tmp_dat);
                    for (l=0;l<nCCC;l++)
                        if (tmp_dat[0]==cCCC[3*l+0] && tmp_dat[1]==cCCC[3*l+1] && tmp_dat[2]==cCCC[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA+nCAA+nCCA) {
                            cluster_phi[clusters_to_use[l]] = cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phi[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
//    INFO_Printf("Found all CCC relationship\n");

    double energy_cluster = 0.0;
    for (i=0;i<nonzeroeci;i++) {
//        printf("cluster to use = %d, eci_opt = %e, phi = %d\n",clusters_to_use[i],eci_opt[clusters_to_use[i]],cluster_phi[clusters_to_use[i]]);
        energy_cluster += eci_opt[clusters_to_use[i]]*cluster_phi[clusters_to_use[i]];
    }
//    INFO_Printf("CEM Calculation is done\n");

    return(energy_cluster);
    
}

/**************************************************************************/
/* Cluster Expansion Method version 2 - obtain clusters, AA, CC, and CA   */    
/**************************************************************************/
int YSZFrame::obtain_AA(int ix, int iy, int iz, int jx, int jy, int jz)
{
    int i;
    int classAA = -1;
    int kind_t = -1;
    int disp_dat[3];
    int dxij, dyij, dzij;
    int disp_dat2;
    int ixyz=0;
    
    dxij = (int) mymod(jx-ix+NCellX/2,NCellX)-NCellX/2;    //due to PBCs
    dyij = (int) mymod(jy-iy+NCellY/2,NCellY)-NCellY/2;
    dzij = (int) mymod(jz-iz+NCellZ/2,NCellZ)-NCellZ/2;
    disp_dat[0] = abs(dxij);
    disp_dat[1] = abs(dyij);
    disp_dat[2] = abs(dzij);
    mysort2(disp_dat);

    disp_dat2 = disp_dat[0]*disp_dat[0]+disp_dat[1]*disp_dat[1]+disp_dat[2]*disp_dat[2];
    if (disp_dat2*1.0 < cluster_rcut*cluster_rcut) {     // this is a scheme using rcut
//    if (disp_dat[2]*1.0 < cluster_dcut) {
        ixyz = (int) mymod(ix+iy+iz,4);
        if ( ixyz == 1 )
            if ( dxij*dyij*dzij < 0 )
                kind_t = 1;
            else
                kind_t = 0;
        else if ( ixyz == 3 )
            if ( dxij*dyij*dzij > 0 )
                kind_t = 1;
            else
                kind_t = 0;
        else
            ERROR("Anion position is not compatible in obtain_AA");
        
        if ( disp_dat[0]==2 && disp_dat[1]==2 && disp_dat[2]==4 )
            kind_t = 0;
        else if ( disp_dat[0]==4 && disp_dat[1]==4 && disp_dat[2]==4 )
            kind_t = 0;
        else if ( disp_dat[0]==2 && disp_dat[1]==4 && disp_dat[2]==6 )
            kind_t = 0;
        else if ( disp_dat[0]==4 && disp_dat[1]==6 && disp_dat[2]==6 )
            kind_t = 0;
        
        for (i=0;i<nAA;i++)
            if ( disp_dat[0]==cAA[4*i+0] && disp_dat[1]==cAA[4*i+1] && disp_dat[2]==cAA[4*i+2] )
                if (kind_t==cAA[4*i+3]) {
                    classAA = i;
                    return(classAA);
                }
    }
    
    return(classAA);
        
}

int YSZFrame::obtain_CC(int ix, int iy, int iz, int jx, int jy, int jz)
{
    int i;
    int classCC = -1;
    int kind_t = 0;
    int disp_dat[3];
    int dxij, dyij, dzij;
    int disp_dat2;
    
    dxij = (int) mymod(jx-ix+NCellX/2,NCellX)-NCellX/2;    //due to PBCs
    dyij = (int) mymod(jy-iy+NCellY/2,NCellY)-NCellY/2;
    dzij = (int) mymod(jz-iz+NCellZ/2,NCellZ)-NCellZ/2;
    disp_dat[0] = abs(dxij);
    disp_dat[1] = abs(dyij);
    disp_dat[2] = abs(dzij);
    mysort2(disp_dat);

    disp_dat2 = disp_dat[0]*disp_dat[0]+disp_dat[1]*disp_dat[1]+disp_dat[2]*disp_dat[2];
    if (disp_dat2*1.0 < cluster_rcut*cluster_rcut) {     // this is a scheme using rcut
//    if (disp_dat[2]*1.0 < cluster_dcut) {
        for (i=0;i<nCC;i++)
            if ( disp_dat[0]==cCC[4*i+0] && disp_dat[1]==cCC[4*i+1] && disp_dat[2]==cCC[4*i+2] )
                if (kind_t==cCC[4*i+3]) {
                    classCC = i;
                    return(classCC);
                }
    }
    
    return(classCC);

}

int YSZFrame::obtain_CA(int ix, int iy, int iz, int jx, int jy, int jz)
{
    int i;
    int classCA = -1;
    int kind_t = 0;
    int disp_dat[3];
    int dxij, dyij, dzij;
    int disp_dat2;
    
    dxij = (int) mymod(jx-ix+NCellX/2,NCellX)-NCellX/2;    //due to PBCs
    dyij = (int) mymod(jy-iy+NCellY/2,NCellY)-NCellY/2;
    dzij = (int) mymod(jz-iz+NCellZ/2,NCellZ)-NCellZ/2;
    disp_dat[0] = abs(dxij);
    disp_dat[1] = abs(dyij);
    disp_dat[2] = abs(dzij);
    mysort2(disp_dat);

    disp_dat2 = disp_dat[0]*disp_dat[0]+disp_dat[1]*disp_dat[1]+disp_dat[2]*disp_dat[2];
    if (disp_dat2*1.0 < cluster_rcut*cluster_rcut) {     // this is a scheme using rcut
//    if (disp_dat[2]*1.0 < cluster_dcut) {
        for (i=0;i<nCA;i++)
            if ( disp_dat[0]==cCA[4*i+0] && disp_dat[1]==cCA[4*i+1] && disp_dat[2]==cCA[4*i+2] )
                if (kind_t==cCA[4*i+3]) {
                    classCA = i;
                    return(classCA);
                }    
    }

    return(classCA);
    
}

/*****************************************************************************************************/
/* Cluster Expansion Method version 2 - update the energy, requires cal_energy_by_CEM_initial prerun */
/*****************************************************************************************************/
double YSZFrame::cal_energy_by_CEM_update(int vac_i, int vac_j)
{
    int i;
    int Anid1, Anid2;
    int spini, spinj, spink;
    double eci_diff = 0.0;
//    int cluster_phix[neci],n_cluster_phix[neci];

    if ( vac_i<0 || vac_j<0 )
        return 0.0;
    
    // vac_i(spin=-1) is the current vacancy, vac_j(spin=1) is the new vac.
    for (i=0;i<NAnion;i++)
        if (anionlist[i] == vac_i) {
            Anid1 = i;
            break;
        }
    for (i=0;i<NAnion;i++)
        if (anionlist[i] == vac_j) {
            Anid2 = i;
            break;
        }
    // update for Anid1
    for (i=0;i<map_AA[0];i++) {
        if (map_An_AA[map_AA[0]*Anid1+i] == 1) {
            spini = cluster_spin[map_AA[3*i+1]];
            spinj = cluster_spin[map_AA[3*i+2]];
            eci_diff -=  2*eci_opt[map_AA[3*i+3]]*spini*spinj;
        }
//        printf("1.mapAAid=%d, cluster=%d, eci_diff = %e, spini=%d, spinj=%d\n",map_An_AA[map_AA[0]*Anid1+i],map_AA[3*map_An_AA[map_AA[0]*Anid1+i]+3],eci_diff,spini,spinj);
    }
//    for (i=0;i<map_CC[0];i++) {
//        if (map_An_CC[map_CC[0]*Anid1+i] == 1) {
//            spini = cluster_spin[map_CC[3*i+1]];
//            spinj = cluster_spin[map_CC[3*i+2]];
//            eci_diff -=  2*eci_opt[map_CC[3*i+3]]*spini*spinj;
//        }
//    }
    for (i=0;i<map_CA[0];i++) {
        if (map_An_CA[map_CA[0]*Anid1+i] == 1) {
            spini = cluster_spin[map_CA[3*i+1]];
            spinj = cluster_spin[map_CA[3*i+2]];
            eci_diff -=  2*eci_opt[map_CA[3*i+3]]*spini*spinj;
        }
    }
    for (i=0;i<map_AAA[0];i++) {
        if (map_An_AAA[map_AAA[0]*Anid1+i] == 1) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff -=  2*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
        }
    }
    for (i=0;i<map_CAA[0];i++) {
        if (map_An_CAA[map_CAA[0]*Anid1+i] == 1) {
            spini = cluster_spin[map_CAA[4*i+1]];
            spinj = cluster_spin[map_CAA[4*i+2]];
            spink = cluster_spin[map_CAA[4*i+3]];
            eci_diff -=  2*eci_opt[map_CAA[4*i+4]]*spini*spinj*spink;
        }
    }
    for (i=0;i<map_CCA[0];i++) {
        if (map_An_CCA[map_CCA[0]*Anid1+i] == 1) {
            spini = cluster_spin[map_CCA[4*i+1]];
            spinj = cluster_spin[map_CCA[4*i+2]];
            spink = cluster_spin[map_CCA[4*i+3]];
            eci_diff -=  2*eci_opt[map_CCA[4*i+4]]*spini*spinj*spink;
        }
    }
//    for (i=0;i<map_CCC[0];i++) {
//        if (map_An_CCC[map_CCC[0]*Anid1+i] == 1) {
//            spini = cluster_spin[map_CCC[4*i+1]];
//            spinj = cluster_spin[map_CCC[4*i+2]];
//            spink = cluster_spin[map_CCC[4*i+3]];
//            eci_diff -=  2*eci_opt[map_CCC[4*i+4]]*spini*spinj*spink;
//        }
//    }
    
    // update for Anid2
    for (i=0;i<map_AA[0];i++) {
        if (map_An_AA[map_AA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_AA[3*i+1]];
            spinj = cluster_spin[map_AA[3*i+2]];
            eci_diff -=  2*eci_opt[map_AA[3*i+3]]*spini*spinj;
        }
//        printf("2.mapAAid=%d, cluster=%d, eci_diff = %e, spini=%d, spinj=%d\n",map_An_AA[map_AA[0]*Anid1+i],map_AA[3*map_An_AA[map_AA[0]*Anid1+i]+3],eci_diff,spini,spinj);
    }
//    for (i=0;i<map_CC[0];i++) {
//        if (map_An_CC[map_CC[0]*Anid2+i] == 1) {
//            spini = cluster_spin[map_CC[3*i+1]];
//            spinj = cluster_spin[map_CC[3*i+2]];
//            eci_diff -=  2*eci_opt[map_CC[3*i+3]]*spini*spinj;
//        }
//    }
    for (i=0;i<map_CA[0];i++) {
        if (map_An_CA[map_CA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_CA[3*i+1]];
            spinj = cluster_spin[map_CA[3*i+2]];
            eci_diff -=  2*eci_opt[map_CA[3*i+3]]*spini*spinj;
        }
    }
    for (i=0;i<map_AAA[0];i++) {
        if (map_An_AAA[map_AAA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff -=  2*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
        }
    }
    for (i=0;i<map_CAA[0];i++) {
        if (map_An_CAA[map_CAA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_CAA[4*i+1]];
            spinj = cluster_spin[map_CAA[4*i+2]];
            spink = cluster_spin[map_CAA[4*i+3]];
            eci_diff -=  2*eci_opt[map_CAA[4*i+4]]*spini*spinj*spink;
        }
    }
    for (i=0;i<map_CCA[0];i++) {
        if (map_An_CCA[map_CCA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_CCA[4*i+1]];
            spinj = cluster_spin[map_CCA[4*i+2]];
            spink = cluster_spin[map_CCA[4*i+3]];
            eci_diff -=  2*eci_opt[map_CCA[4*i+4]]*spini*spinj*spink;
        }
    }
//    for (i=0;i<map_CCC[0];i++) {
//        if (map_An_CCC[map_CCC[0]*Anid2+i] == 1) {
//            spini = cluster_spin[map_CCC[4*i+1]];
//            spinj = cluster_spin[map_CCC[4*i+2]];
//            spink = cluster_spin[map_CCC[4*i+3]];
//            eci_diff -=  2*eci_opt[map_CCC[4*i+4]]*spini*spinj*spink;
//        }
//    }
    // double-counting Anid1 & Anid2
    for (i=0;i<map_AA[0];i++) {
        if (map_An_AA[map_AA[0]*Anid1+i] == 1 && map_An_AA[map_AA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_AA[3*i+1]];
            spinj = cluster_spin[map_AA[3*i+2]];
            eci_diff +=  4*eci_opt[map_AA[3*i+3]]*spini*spinj;
        }
    }
    for (i=0;i<map_CA[0];i++) {
        if (map_An_CA[map_CA[0]*Anid1+i] == 1 && map_An_CA[map_CA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_CA[3*i+1]];
            spinj = cluster_spin[map_CA[3*i+2]];
            eci_diff +=  4*eci_opt[map_CA[3*i+3]]*spini*spinj;
        }
    }
    for (i=0;i<map_AAA[0];i++) {
        if (map_An_AAA[map_AAA[0]*Anid1+i] == 1 && map_An_AAA[map_AAA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff +=  4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
        }
    }
    for (i=0;i<map_CAA[0];i++) {
        if (map_An_CAA[map_CAA[0]*Anid1+i] == 1 && map_An_CAA[map_CAA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_CAA[4*i+1]];
            spinj = cluster_spin[map_CAA[4*i+2]];
            spink = cluster_spin[map_CAA[4*i+3]];
            eci_diff +=  4*eci_opt[map_CAA[4*i+4]]*spini*spinj*spink;
        }
    }
    for (i=0;i<map_CCA[0];i++) {
        if (map_An_CCA[map_CCA[0]*Anid1+i] == 1 && map_An_CCA[map_CCA[0]*Anid2+i] == 1) {
            spini = cluster_spin[map_CCA[4*i+1]];
            spinj = cluster_spin[map_CCA[4*i+2]];
            spink = cluster_spin[map_CCA[4*i+3]];
            eci_diff +=  4*eci_opt[map_CCA[4*i+4]]*spini*spinj*spink;
        }
    }
    
    

#if 0
    double eci_diff2 = 0.0;
    for (i=0;i<neci;i++)
        cluster_phix[i] = cluster_phi[i];
    for (i=0;i<neci;i++)
        n_cluster_phix[i] = n_cluster_phi[i];

    // vac_i(spin=-1) is the current vacancy, vac_j(spin=1) is the new vac.    
    for (i=0;i<map_AA[0];i++) 
        if (map_AA[3*i+1]==vac_i || map_AA[3*i+2]==vac_i) {
            spini = cluster_spin[map_AA[3*i+1]];
            spinj = cluster_spin[map_AA[3*i+2]];
            eci_diff2 -=  2*eci_opt[map_AA[3*i+3]]*spini*spinj;
            cluster_phix[map_AA[3*i+3]] -= 2*spini*spinj;
        }
    for (i=0;i<map_CC[0];i++) 
        if (map_CC[3*i+1]==vac_i || map_CC[3*i+2]==vac_i) {
            spini = cluster_spin[map_CC[3*i+1]];
            spinj = cluster_spin[map_CC[3*i+2]];
            eci_diff2 -=  2*eci_opt[map_CC[3*i+3]]*spini*spinj;
            cluster_phix[map_CC[3*i+3]] -= 2*spini*spinj;
    }
    for (i=0;i<map_CA[0];i++) 
        if (map_CA[3*i+1]==vac_i || map_CA[3*i+2]==vac_i) {
            spini = cluster_spin[map_CA[3*i+1]];
            spinj = cluster_spin[map_CA[3*i+2]];
            eci_diff2 -=  2*eci_opt[map_CA[3*i+3]]*spini*spinj;
            cluster_phix[map_CA[3*i+3]] -= 2*spini*spinj;
        }
//    INFO_Printf("vac_i=%d,vac_j=%d\n",vac_i,vac_j);
    for (i=0;i<map_AAA[0];i++)
        if (map_AAA[4*i+1]==vac_i || map_AAA[4*i+2]==vac_i || map_AAA[4*i+3]==vac_i) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 -= 2*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] -= 2*spini*spinj*spink;
//            INFO_Printf("map_AAA[%d,:]=[%d %d %d], spin=(%d,%d,%d), cluster=%d\n",i,map_AAA[4*i+1],map_AAA[4*i+2],map_AAA[4*i+3],spini,spinj,spink,map_AAA[4*i+4]);
        }
    for (i=0;i<map_AA[0];i++) 
        if (map_AA[3*i+1]==vac_j || map_AA[3*i+2]==vac_j) {
            spini = cluster_spin[map_AA[3*i+1]];
            spinj = cluster_spin[map_AA[3*i+2]];
            eci_diff2 -=  2*eci_opt[map_AA[3*i+3]]*spini*spinj;
            cluster_phix[map_AA[3*i+3]] -= 2*spini*spinj;
        }
    for (i=0;i<map_CC[0];i++) 
        if (map_CC[3*i+1]==vac_j || map_CC[3*i+2]==vac_j) {
            spini = cluster_spin[map_CC[3*i+1]];
            spinj = cluster_spin[map_CC[3*i+2]];
            eci_diff2 -=  2*eci_opt[map_CC[3*i+3]]*spini*spinj;
            cluster_phix[map_CC[3*i+3]] -= 2*spini*spinj;
    }
    for (i=0;i<map_CA[0];i++) 
        if (map_CA[3*i+1]==vac_j || map_CA[3*i+2]==vac_j) {
            spini = cluster_spin[map_CA[3*i+1]];
            spinj = cluster_spin[map_CA[3*i+2]];
            eci_diff2 -=  2*eci_opt[map_CA[3*i+3]]*spini*spinj;
            cluster_phix[map_CA[3*i+3]] -= 2*spini*spinj;
        }
    for (i=0;i<map_AAA[0];i++)
        if (map_AAA[4*i+1]==vac_j || map_AAA[4*i+2]==vac_j || map_AAA[4*i+3]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 -= 2*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] -= 2*spini*spinj*spink;
        }

    for (i=0;i<map_AAA[0];i++) {
        if (map_AAA[4*i+1]==vac_i && map_AAA[4*i+2]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 += 4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] += 4*spini*spinj*spink;
        }
        if (map_AAA[4*i+2]==vac_i && map_AAA[4*i+1]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 += 4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] += 4*spini*spinj*spink;
        }
        if (map_AAA[4*i+1]==vac_i && map_AAA[4*i+3]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 += 4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] += 4*spini*spinj*spink;
        }
        if (map_AAA[4*i+3]==vac_i && map_AAA[4*i+1]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 += 4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] += 4*spini*spinj*spink;
        }
        if (map_AAA[4*i+2]==vac_i && map_AAA[4*i+3]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 += 4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] += 4*spini*spinj*spink;
        }
        if (map_AAA[4*i+3]==vac_i && map_AAA[4*i+2]==vac_j) {
            spini = cluster_spin[map_AAA[4*i+1]];
            spinj = cluster_spin[map_AAA[4*i+2]];
            spink = cluster_spin[map_AAA[4*i+3]];
            eci_diff2 += 4*eci_opt[map_AAA[4*i+4]]*spini*spinj*spink;
            cluster_phix[map_AAA[4*i+4]] += 4*spini*spinj*spink;
        }
    }

    double energy_update = 0.0;
    for (i=0;i<nonzeroeci;i++) {
        INFO_Printf("cluster to use = %d, eci_opt = %e, phi = %d, n_phi = %d\n",clusters_to_use[i],eci_opt[clusters_to_use[i]],cluster_phix[clusters_to_use[i]],n_cluster_phix[clusters_to_use[i]]);
        energy_update += eci_opt[clusters_to_use[i]]*cluster_phix[clusters_to_use[i]];
    }
    INFO_Printf("eci_diff = %e total energy = %e eci_diff2 = %e\n",eci_diff,energy_update, eci_diff2);
//    INFO_Printf("eci_diff = %e eci_diff2 = %e\n",eci_diff,eci_diff2);
#endif

#if 0
    
    double energy_method1 = 0.0;
    double energy_method2 = 0.0;
    int ii,jj,j,kk,k,l;
    int tmp_map,tmp1,tmp2,tmp3,tmp4;    
    int tmp_dat[3];

    for (i=0;i<neci;i++)
        cluster_phix[i] = 0;
    for (i=0;i<neci;i++)
        n_cluster_phix[i] = 0;
    
    cluster_spin[vac_i] = -cluster_spin[vac_i];
    cluster_spin[vac_j] = -cluster_spin[vac_j];

    for (l=0;l<nonzeroeci;l++) 
        if (clusters_to_use[l] == 0) {
            cluster_phix[0] = 1;
            n_cluster_phix[0] = 1;
            break;
        }
    for (l=0;l<nonzeroeci;l++) 
        if (clusters_to_use[l] == 1) {
            cluster_phix[1] = NAnion-2*num_vac;
            n_cluster_phix[1] = NAnion;
            break;
        }
    for (l=0;l<nonzeroeci;l++) 
        if (clusters_to_use[l] == 2) {
            cluster_phix[2] = NCation-2*2*num_vac;
            n_cluster_phix[2] = NCation;
            break;
        }   
    
    for (ii=0;ii<NAnion-1;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion;jj++) {
            j = anionlist[jj];            
            tmp_map = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3) {
                        cluster_phix[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            tmp_map = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA) {
                        cluster_phix[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }                
        }
    }
    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion;jj++) {
            j = anionlist[jj];
            tmp_map = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA+nCC) {
                        cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
    for (ii=0;ii<NAnion-2;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_AA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+100; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+100;
                    mysort2(tmp_dat);
                    for (l=0;l<nAAA;l++)
                        if (tmp_dat[0]==cAAA[3*l+0] && tmp_dat[1]==cAAA[3*l+1] && tmp_dat[2]==cAAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA) {
                            cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phix[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+200; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+200;
                    mysort2(tmp_dat);
                    for (l=0;l<nCAA;l++)
                        if (tmp_dat[0]==cCAA[3*l+0] && tmp_dat[1]==cCAA[3*l+1] && tmp_dat[2]==cCAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA) {
                            cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phix[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            for (kk=0;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_CA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1; tmp_dat[1]=tmp2+200; tmp_dat[2]=tmp3+200;
                    mysort2(tmp_dat);
                    for (l=0;l<nCCA;l++)
                        if (tmp_dat[0]==cCCA[3*l+0] && tmp_dat[1]==cCCA[3*l+1] && tmp_dat[2]==cCCA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA+nCAA) {
                            cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phix[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
    for (ii=0;ii<NCation-2;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation-1;jj++) {
            j = cationlist[jj];
            for (kk=jj+1;kk<NCation;kk++) {
                k = cationlist[kk];
                tmp1 = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_CC(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_CC(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1; tmp_dat[1]=tmp2; tmp_dat[2]=tmp3;
                    mysort2(tmp_dat);
                    for (l=0;l<nCCC;l++)
                        if (tmp_dat[0]==cCCC[3*l+0] && tmp_dat[1]==cCCC[3*l+1] && tmp_dat[2]==cCCC[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA+nAAA+nCAA+nCCA) {
                            cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phix[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }

    for (i=0;i<nonzeroeci;i++) {
        INFO_Printf("cluster to use = %d, eci_opt = %e, phi = %d, n_phi = %d\n",clusters_to_use[i],eci_opt[clusters_to_use[i]],cluster_phix[clusters_to_use[i]],n_cluster_phix[clusters_to_use[i]]);
        energy_method2 += eci_opt[clusters_to_use[i]]*cluster_phix[clusters_to_use[i]];
    }
    INFO_Printf("eci_diff = %e total energy = %e\n",eci_diff,energy_method2);

    cluster_spin[vac_i] = -cluster_spin[vac_i];
    cluster_spin[vac_j] = -cluster_spin[vac_j];
/*************************************************************************************************************/
/*
    for (i=0;i<neci;i++)
        cluster_phix[i] = 0;
    for (i=0;i<neci;i++)
        n_cluster_phix[i] = 0;
    cluster_spin[vac_i] = -cluster_spin[vac_i];
    
    for (ii=0;ii<NAnion-1;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion;jj++) {
            j = anionlist[jj];            
            tmp_map = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3) {
                        cluster_phix[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            tmp_map = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA) {
                        cluster_phix[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }                
        }
    }
    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion;jj++) {
            j = anionlist[jj];
            tmp_map = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA+nCC) {
                        cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
    for (ii=0;ii<NAnion-2;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_AA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+100; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+100;
                    mysort2(tmp_dat);
                    for (l=0;l<nAAA;l++)
                        if (tmp_dat[0]==cAAA[3*l+0] && tmp_dat[1]==cAAA[3*l+1] && tmp_dat[2]==cAAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA) {
                            cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phix[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }
    
    cluster_spin[vac_j] = -cluster_spin[vac_j];
    for (ii=0;ii<NAnion-1;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion;jj++) {
            j = anionlist[jj];            
            tmp_map = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3) {
                        cluster_phix[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
    for (ii=0;ii<NCation-1;ii++) {
        i = cationlist[ii];
        for (jj=ii+1;jj<NCation;jj++) {
            j = cationlist[jj];
            tmp_map = obtain_CC(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA) {
                        cluster_phix[clusters_to_use[l]]  += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }                
        }
    }
    for (ii=0;ii<NCation;ii++) {
        i = cationlist[ii];
        for (jj=0;jj<NAnion;jj++) {
            j = anionlist[jj];
            tmp_map = obtain_CA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
            if ( tmp_map >=0 )
                for (l=0;l<nonzeroeci;l++) 
                    if (clusters_to_use[l] == tmp_map+3+nAA+nCC) {
                        cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j];
                        n_cluster_phix[clusters_to_use[l]]++;
                        break;
                    }
        }
    }
    for (ii=0;ii<NAnion-2;ii++) {
        i = anionlist[ii];
        for (jj=ii+1;jj<NAnion-1;jj++) {
            j = anionlist[jj];
            for (kk=jj+1;kk<NAnion;kk++) {
                k = anionlist[kk];
                tmp1 = obtain_AA(_RI[3*i+0],_RI[3*i+1],_RI[3*i+2],_RI[3*j+0],_RI[3*j+1],_RI[3*j+2]);
                tmp2 = obtain_AA(_RI[3*j+0],_RI[3*j+1],_RI[3*j+2],_RI[3*k+0],_RI[3*k+1],_RI[3*k+2]);
                tmp3 = obtain_AA(_RI[3*k+0],_RI[3*k+1],_RI[3*k+2],_RI[3*i+0],_RI[3*i+1],_RI[3*i+2]);
                tmp4 = -1;
                if (tmp1>=0 && tmp2>=0 && tmp3>=0) {
                    tmp_dat[0]=tmp1+100; tmp_dat[1]=tmp2+100; tmp_dat[2]=tmp3+100;
                    mysort2(tmp_dat);
                    for (l=0;l<nAAA;l++)
                        if (tmp_dat[0]==cAAA[3*l+0] && tmp_dat[1]==cAAA[3*l+1] && tmp_dat[2]==cAAA[3*l+2]) {
                            tmp4 = l;
                            break;
                        }
                }
                if ( tmp4 >=0 ) 
                    for (l=0;l<nonzeroeci;l++) 
                        if (clusters_to_use[l] == tmp4+3+nAA+nCC+nCA) {
                            cluster_phix[clusters_to_use[l]] += cluster_spin[i]*cluster_spin[j]*cluster_spin[k];
                            n_cluster_phix[clusters_to_use[l]]++;
                            break;
                        }                
            }
        }
    }

    cluster_spin[vac_i] = -cluster_spin[vac_i];
    cluster_spin[vac_j] = -cluster_spin[vac_j];
    
    for (i=0;i<nonzeroeci;i++) {
        INFO_Printf("cluster to use = %d, eci_opt = %e, phi = %d, n_phi = %d\n",clusters_to_use[i],eci_opt[clusters_to_use[i]],cluster_phix[clusters_to_use[i]],n_cluster_phix[clusters_to_use[i]]);
        energy_method1 += eci_opt[clusters_to_use[i]]*cluster_phix[clusters_to_use[i]];
    }
    INFO_Printf("eci_diff = %e total energy = %e\n",eci_diff,energy_method1);
*/
/*****************************************************************************************/
    
#endif
    
    
    return(eci_diff);
    
}

void YSZFrame::mysort2(int *targetarray)
{
    int tmp1, tmp2, tmp3;
    
    if ( *(targetarray+0) <= *(targetarray+1) && *(targetarray+1) <= *(targetarray+2) ) {
        tmp1 = *(targetarray+0);
        tmp2 = *(targetarray+1);
        tmp3 = *(targetarray+2);
    }
    else if ( *(targetarray+0) <= *(targetarray+2) && *(targetarray+2) <= *(targetarray+1) ) {
        tmp1 = *(targetarray+0);
        tmp2 = *(targetarray+2);
        tmp3 = *(targetarray+1);
    }
    else if ( *(targetarray+1) <= *(targetarray+0) && *(targetarray+0) <= *(targetarray+2) ) {
        tmp1 = *(targetarray+1);
        tmp2 = *(targetarray+0);
        tmp3 = *(targetarray+2);
    }
    else if ( *(targetarray+1) <= *(targetarray+2) && *(targetarray+2) <= *(targetarray+0) ) {
        tmp1 = *(targetarray+1);
        tmp2 = *(targetarray+2);
        tmp3 = *(targetarray+0);
    }
    else if ( *(targetarray+2) <= *(targetarray+0) && *(targetarray+0) <= *(targetarray+1) ) {
        tmp1 = *(targetarray+2);
        tmp2 = *(targetarray+0);
        tmp3 = *(targetarray+1);
    }
    else if ( *(targetarray+2) <= *(targetarray+1) && *(targetarray+1) <= *(targetarray+0) ) {
        tmp1 = *(targetarray+2);
        tmp2 = *(targetarray+1);
        tmp3 = *(targetarray+0);
    }

    *(targetarray+0) = tmp1;
    *(targetarray+1) = tmp2;
    *(targetarray+2) = tmp3;

}

/********************************************************************/
/* File input and output */
/********************************************************************/

int YSZFrame::set_corr_filecounter(int i)
{
    cfs.setcount(i);
    INFO_Printf("corrs file counter starts at %d\n",i);
    cfl.setcount(i);
    INFO_Printf("corrl file counter starts at %d\n",i);
    return 0;
}

int YSZFrame::set_diff_filecounter(int i)
{
    df.setcount(i);
    INFO_Printf("diff file counter starts at %d\n",i);
    return 0;
}

int YSZFrame::set_hist_filecounter(int i)
{
    hf.setcount(i);
    INFO_Printf("hist file counter starts at %d\n",i);
    return 0;
}

int YSZFrame::set_atomeye_filecounter(int i)
{
    mf.setcount(i);
    INFO_Printf("atomeye file counter starts at %d\n",i);
    return 0;
}

/* Correlation Function file */
char * CorrFileShort::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","correlations function file");
    return tmp;
}

int CorrFileShort::writeblock(void *p)
{
    int i;
    char tmp[500];
    YSZFrame &d=*((YSZFrame *)p);

/*    sprintf(tmp,"%12d %20.13e %20.13e %20.13e %20.13e\n",
            d.curstep,                  
            d.treal,                    
            d.enpos[0],                 
            d.enpos[1],                 
            d.enpos[2]                  
        );
*/    
    for (i=0;i<d.ncorrsteps;i++) {
//        INFO_Printf("corr[%d]=%e\n",i,d.corrs[i]);
        sprintf(tmp,"%6d %20.13e\n",i,d.corrs[i]);
        *f<<tmp;
    }
//    sprintf(tmp,"curstep = %d corr[0]=%e\n",d.curstep,d.corr[0]);
//    *f<<tmp;
    

    return 0;
}
char * CorrFileLong::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","correlationl function file");
    return tmp;
}

int CorrFileLong::writeblock(void *p)
{
    int i;
    char tmp[500];
    YSZFrame &d=*((YSZFrame *)p);

/*    sprintf(tmp,"%12d %20.13e %20.13e %20.13e %20.13e\n",
            d.curstep,                  
            d.treal,                    
            d.enpos[0],                 
            d.enpos[1],                 
            d.enpos[2]                  
        );
*/    
    for (i=0;i<d.ncorrstepl;i++) {
//        INFO_Printf("corr[%d]=%e\n",i,d.corrl[i]);
        sprintf(tmp,"%6d %20.13e\n",i,d.corrl[i]);
        *f<<tmp;
    }
//    sprintf(tmp,"curstep = %d corr[0]=%e\n",d.curstep,d.corr[0]);
//    *f<<tmp;
    

    return 0;
}

char * ACFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","AC response function file");
    return tmp;
}

int ACFile::writeblock(void *p)
{
    int i;
    char tmp[500];
    YSZFrame &d=*((YSZFrame *)p);

    for (i=0;i<d.nextdt;i++) {
        sprintf(tmp,"%20.13e\n",d.enposwext[i]);
        *f<<tmp;
    }
    

    return 0;
}

char * DiffFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","Diffusion Coefficient file");
    return tmp;
}

int DiffFile::writeblock(void *p)
{
#if 1
    int i;
    char tmp[500];
    YSZFrame &d=*((YSZFrame *)p);
    
    int d_x,d_y,d_z;
    double Dsquare, SelfDiffusivity,SelfDiffusivity_com, Dsquare_com, dx_com,dy_com,dz_com;
    double SelfDiffusivity_x;
    Dsquare = 0.0; 
    dx_com = 0.0; dy_com=0.0; dz_com=0.0;
    for (i=0;i<d.num_vac;i++) {
        d_x = d.rx[i] - d._RI[3*d.vaclist0[i]+0];
        d_y = d.ry[i] - d._RI[3*d.vaclist0[i]+1];
        d_z = d.rz[i] - d._RI[3*d.vaclist0[i]+2];
        Dsquare += d_x*d_x*1.0 + d_y*d_y*1.0 + d_z*d_z*1.0;
        dx_com += d_x;
        dy_com += d_y;
        dz_com += d_z;
    }
    dx_com /= d.num_vac;
    dy_com /= d.num_vac;
    dz_com /= d.num_vac;    
    Dsquare_com = dx_com*dx_com + dy_com*dy_com + dz_com*dz_com;
    SelfDiffusivity = Dsquare / (d.num_vac*6.0*d.treal);
    SelfDiffusivity_com = Dsquare_com / (6.0*d.treal);
    SelfDiffusivity_x = dx_com*dx_com / (2.0*d.treal);
    sprintf(tmp,"%6.4e %6.4e %6.4e %6.10e %6.4e\n",dx_com,dy_com,dz_com,d.treal,Dsquare);
    *f<<tmp;

//    sprintf(tmp,"%6.10e  %6.10e  %6.10e\n",SelfDiffusivity,SelfDiffusivity_com,SelfDiffusivity_x);
//    *f<<tmp;
#endif
    return 0;
}

char * HistFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","Ebarrier Histogram file");
    return tmp;
}

int HistFile::writeblock(void *p)
{
    int i;
    char tmp[500];
    YSZFrame &d=*((YSZFrame *)p);

    for (i=0;i<1000;i++) {
        sprintf(tmp,"%20.13e\n",d.histEb[i]);
        *f<<tmp;
    }    
    return 0;
}

char * AtomEyeFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","Atomeye data file ");
    return tmp;
}

int AtomEyeFile::writeblock(void *p)
{    
    int i;
    char tmp[500];
    YSZFrame &d=*((YSZFrame *)p);

    sprintf(tmp,"Number of particles = %d\n",d.num_vac);
    *f<<tmp;
    sprintf(tmp,"A = 1.0 Angstrom (basic length-scale)\n");
    *f<<tmp;
    sprintf(tmp,"H0(1,1) = 25.7 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(1,2) = 0 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(1,3) = 0 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(2,1) = 0 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(2,2) = 25.7 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(2,3) = 0 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(3,1) = 0 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(3,2) = 0 A\n");
    *f<<tmp;
    sprintf(tmp,"H0(3,3) = 25.7 A\n");
    *f<<tmp;
    sprintf(tmp,".NO_VELOCITY.\n");
    *f<<tmp;
    sprintf(tmp,"entry_count = 1\n");
    *f<<tmp;
    sprintf(tmp,"16.0\n");
    *f<<tmp;
    sprintf(tmp,"O\n");
    *f<<tmp;
    
    for (i=0;i<d.num_vac;i++) {
        sprintf(tmp,"%f %f %f\n",fmod(d.rx[i],20)*(5.14/4),fmod(d.ry[i],20)*(5.14/4),fmod(d.rz[i],20)*(5.14/4));
        *f<<tmp;
    }
    return 0;
}

#ifdef _TEST

/* Main Program Begins */
class YSZFrame sim;

#include "main.cpp"
#endif//_TEST

