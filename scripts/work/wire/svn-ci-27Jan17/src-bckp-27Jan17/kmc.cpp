/*
  kmc.h
  by Wei Cai  caiwei@stanford.edu, Eunseok Lee  euniv@stanford.edu
  Last Modified : Wed Jan 20 16:31:59 2010

  FUNCTION  :  Easy-to-use KMC simulation Framefork

  Featuring :  1. Scripting input
               2. X-window display
               3. automatic simulation log
               4. convenient configuration and output file handle
               
  This is a single CPU code.
*/

#include "kmc.h"


/********************************************************************/
/* Begining of class KMCFrame definition */
/********************************************************************/


/********************************************************************/
/* Initialize script file parser */
/********************************************************************/
void KMCFrame::initparser()
{
    char s[100];

    /* Parser */
    bindvar("myname",myname,STRING);
    bindvar("command",command,STRING);
    bindvar("input",input,DOUBLE);
    bindvar("output",output,STRING);
    bindvar("output_fmt",output_fmt,STRING);


    /* property variables */
    bindvar("NP",&_NP,INT);
    bindvar("H",&(_H[0][0]),DOUBLE);
    bindvar("H0",&(_H0[0][0]),DOUBLE);
    bindvar("OMEGA",&_OMEGA,DOUBLE);
    bindvar("EPOT",&_EPOT,DOUBLE);
    bindvar("T",&_T,DOUBLE);
    bindvar("RLIST",&_RLIST,DOUBLE);
    bindvar("curstep",&curstep,INT);
    
    bindvar("H_11",&(_H[0][0]),DOUBLE);
    bindvar("H_12",&(_H[0][1]),DOUBLE);
    bindvar("H_13",&(_H[0][2]),DOUBLE);
    bindvar("H_21",&(_H[1][0]),DOUBLE);
    bindvar("H_22",&(_H[1][1]),DOUBLE);
    bindvar("H_23",&(_H[1][2]),DOUBLE);
    bindvar("H_31",&(_H[2][0]),DOUBLE);
    bindvar("H_32",&(_H[2][1]),DOUBLE);
    bindvar("H_33",&(_H[2][2]),DOUBLE);

    bindvar("H",&(_H[0][0]),DOUBLE);  
    bindvar("H0",&(_H0[0][0]),DOUBLE);  

    /* Kinetic Monte Carlo */
    bindvar("nspecies",&nspecies,INT);
    bindvar("atommass",_ATOMMASS,DOUBLE);
    bindvar("totalsteps",&totalsteps,INT);
    bindvar("equilsteps",&equilsteps,INT);
    bindvar("timestep",&_TIMESTEP,DOUBLE);
    bindvar("randseed",&_RANDSEED,INT);
    bindvar("select_atom",&select_atom,INT);
    bindvar("select_dir",&select_dir,INT);
    bindvar("NIC",&_NIC,INT);
    bindvar("NNM",&_NNM,INT);
    
    /* Configuration manipulation */
    bindvar("crystalstructure",crystalstructure,STRING);
    bindvar("latticeconst",latticeconst,DOUBLE);
    bindvar("latticesize",latticesize,DOUBLE);
    
    /* File input and output */
    bindvar("savecn",&savecn,INT);
    bindvar("saveprop",&saveprop,INT);
    bindvar("savecnfreq",&savecnfreq,INT);
    bindvar("savepropfreq",&savepropfreq,INT);
    bindvar("printfreq",&printfreq,INT);
    bindvar("filecounter",&filecounter,INT);
    bindvar("incnfile",incnfile,STRING);
    bindvar("finalcnfile",finalcnfile,STRING);
    bindvar("outpropfile",outpropfile,STRING);
    bindvar("intercnfile",intercnfile,STRING);
    
    /* Interface to atomeye and other viewers */
    bindvar("atomeyepath",atomeyepath,STRING);
    bindvar("atomeyerepeat",atomeyerepeat,INT);
    bindvar("atomeyeexe",atomeyeexe,STRING);

    /* Visualization */
    bindvar("win_width",&win_width,INT);
    bindvar("win_height",&win_height,INT);
    bindvar("plotfreq",&plotfreq,INT);
    bindvar("atomradius",atomradius,DOUBLE);
    bindvar("bondradius",&bondradius,DOUBLE);
    bindvar("bondlength",&bondlength,DOUBLE);
    bindvar("bondcolor",bondcolor,STRING);
    bindvar("highlightcolor",highlightcolor,STRING);
    bindvar("fixatomcolor",fixatomcolor,STRING);
    bindvar("backgroundcolor",backgroundcolor,STRING);
    bindvar("rotateangles",rotateangles,DOUBLE);

    bindvar("plot_limits",plot_limits,DOUBLE);
    bindvar("plot_atom_info",&plot_atom_info,INT);
    bindvar("plot_map_pbc",&plot_map_pbc,INT);
    bindvar("plot_color_windows",plot_color_windows,DOUBLE);
    bindvar("plot_color_bar",plot_color_bar,DOUBLE);
    bindvar("energycolorbar",plot_color_bar,DOUBLE); /* for backward compatibility */
    
    bindvar("plot_color_axis",&plot_color_axis,INT);
    bindvar("autowritegiffreq",&autowritegiffreq,INT);

    for(int i=0;i<MAXSPECIES;i++)
    {
        sprintf(s,"element%d",i);
        bindvar(s,element[i],STRING);
        sprintf(s,"atomcolor%d",i);
        bindvar(s,atomcolor[i],STRING);
    }
    bindvar("atomcolor",atomcolor[0],STRING);
    for(int i=0;i<MAXCOLORS;i++)
    {
        sprintf(s,"color%02d",i);
        bindvar(s,colornames[i],STRING);
    }
}    

int KMCFrame::exec(virtual char *name)
{
    if(Organizer::exec(name)==0) return 0;

    /* Parser */
    bindcommand(name,"opendir",welcome());
    bindcommand(name,"print_randseed",INFO_Printf("randseed = %d\n",_RANDSEED));
    bindcommand(name,"runcommand",runcommand());

    /* Kinetic Monte Carlo simulation */
    bindcommand(name,"kmc_step",kmc_step());
    bindcommand(name,"runkmc",runkmc());    
    bindcommand(name,"eval",{refreshneighborlist();eval();});
    bindcommand(name,"multieval",multieval());
    bindcommand(name,"refreshnnlist",refreshneighborlist());
    bindcommand(name,"randomposition",randomposition());
    
    /* Random number generators */
    bindcommand(name,"srand",srand((unsigned int)_RANDSEED));
    bindcommand(name,"srand48",srand48((unsigned int)_RANDSEED));
    bindcommand(name,"srandbytime",srand((unsigned int)time(NULL)));
    bindcommand(name,"srand48bytime",srand48((unsigned int)time(NULL)));
    bindcommand(name,"print_drand48",INFO_Printf("rand#1=%20.14e\n",drand48()));
    
    /* Configuration manipulation */
    bindcommand(name,"makecrystal",makecrystal());

    bindcommand(name,"scaleH",scaleH());
    bindcommand(name,"setH",setH());
    bindcommand(name,"saveH",saveH());
    bindcommand(name,"restoreH",restoreH());
    bindcommand(name,"reorientH",_H=_H.reorient());
    bindcommand(name,"changeH_keepR",shiftbox());
    bindcommand(name,"changeH_keepS",redefinepbc());
    bindcommand(name,"shiftbox",shiftbox());      /* for backward compatibility */
    bindcommand(name,"redefinepbc",redefinepbc());/* for backward compatibility */
    bindcommand(name,"extendbox",extendbox());
    
    bindcommand(name,"moveatom",moveatom());
    bindcommand(name,"movegroup",movegroup());
    bindcommand(name,"pbcshiftatom",pbcshiftatom());
    
    bindcommand(name,"fixatoms_by_ID",fixatoms_by_ID());
    bindcommand(name,"fixatoms_by_position",fixatoms_by_position());
    bindcommand(name,"fixallatoms",fixallatoms());
    bindcommand(name,"freeallatoms",freeallatoms());
    bindcommand(name,"reversefixedatoms",reversefixedatoms());

    bindcommand(name,"setfixedatomsspecies",setfixedatomsspecies());
    bindcommand(name,"setfixedatomsgroup",setfixedatomsgroup());
    bindcommand(name,"reversespecies",reversespecies());
    bindcommand(name,"movefixedatoms",movefixedatoms());
    bindcommand(name,"removefixedatoms",removefixedatoms());
    bindcommand(name,"markremovefixedatoms",markremovefixedatoms());
    bindcommand(name,"removeellipsoid",removeellipsoid());
    bindcommand(name,"removerectbox",removerectbox());
        
    /* File input and output */
    bindcommand(name,"writecn",writefinalcnfile());
    bindcommand(name,"readcn",readcn());
    bindcommand(name,"openintercnfile",openintercnfile());
    bindcommand(name,"writeintercn",writeintercnfile());
    bindcommand(name,"setfilecounter",setfilecounter());
    bindcommand(name,"openpropfile",openpropfile());

    /* Interface with AtomEye and other viewers */
    bindcommand(name,"atomeye",atomeye());
    bindcommand(name,"writeatomeyecfg",writeatomeyecfgfile(finalcnfile));
    bindcommand(name,"convertCNtoCFG",convertCNtoCFG());

    /* Print out energies of current configuration */
    bindcommand(name,"writeENERGY",writeENERGY(finalcnfile));
    bindcommand(name,"writePOSITION",writePOSITION(finalcnfile));
    bindcommand(name,"GnuPlotHistogram",GnuPlotHistogram());    

    /* Visualization */
    bindcommand(name,"openwin",openwin());
    bindcommand(name,"plot",plot());
    bindcommand(name,"alloccolors",alloccolors());
    bindcommand(name,"alloccolorsX",alloccolorsX());
    bindcommand(name,"testcolor",win->testcolor());
    bindcommand(name,"reversergb",win->reversergb());
    bindcommand(name,"rotate",rotate());
    bindcommand(name,"saverot",saverot());
    bindcommand(name,"wintogglepause",wintogglepause());

    return -1;
}

void KMCFrame::runcommand()
{
    char extcomm[400];
    
    LFile::SubHomeDir(command,command);
    INFO("run:  "<<command);
    if(nolog)
        system(command);
    else
    {
        if (ncom==0) system("rm -f B*.log");
        strcpy(extcomm,command);
        strcpy(extcomm+strlen(extcomm)," > B");
        sprintf(extcomm+strlen(extcomm),"%d.log",ncom);
        system(extcomm); 
        sprintf(extcomm,"echo -- '\n'B%d.log: \"%s\"'\n' >> B.log",ncom,command);
        system(extcomm);
        sprintf(extcomm,"cat B%d.log >> B.log",ncom);
        system(extcomm);
        ncom++;
    }
}

void KMCFrame::initvars()
{
    initparser();
    strcpy(command,"echo Hello World");
    
    strcpy(incnfile,"con.cn");
    strcpy(intercnfile,"inter.cn");
    strcpy(finalcnfile,"final.cn");
    strcpy(outpropfile,"prop.out");

    strcpy(bondcolor,"red");
    strcpy(highlightcolor,"purple");

    for(int i=0;i<MAXSPECIES;i++)
    {
        atomradius[i]=0.1;
        strcpy(atomcolor[i],"orange");
        strcpy(element[i],"Mo");
    }
    for(int i=0;i<MAXCOLORS;i++)
    {
        strcpy(colornames[i],"gray50");
    }

    strcpy(element[0],"O");
    strcpy(element[1],"V_O");
    strcpy(element[2],"Y");
    strcpy(element[3],"Zr");

    atomeyerepeat[0]=0;

}


#ifdef _USETCL
/********************************************************************/
/* Initialize Tcl parser */
/********************************************************************/
int KMCFrame::Tcl_AppInit(Tcl_Interp *interp)
{
    /*
     * Tcl_Init reads init.tcl from the Tcl script library.
     */
    if (Tcl_Init(interp) == TCL_ERROR) {
        return TCL_ERROR;
    }
    /*
     * Register application-specific commands.
     */
    Tcl_CreateCommand(interp, "MD++", MD_Cmd,
                      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "MD++_Set", MD_SetVar,
                      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "MD++_Get", MD_GetVar,
                      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    
    /*
      Random_Init(interp);
      Blob_Init(interp);
    */
    return TCL_OK;
}

int Tcl_Main_Wrapper(int argc, char *argv[])
{
    int i, pid;
    char oldargv[1000], newargv[1000], cmd[1000], c, *p, *filename;
    FILE *oldfile, *newfile;

    pid=getpid();
    for(i=1;i<argc;i++)
    {
        strcpy(oldargv,argv[i]);
        /* identify suffix .script */
        p = strrchr(oldargv,'.');
        if(p!=NULL) p = strstr(p,".script");
        if(p!=NULL)
        {
            /* identify file name */
            filename = strrchr(oldargv,'/');
            if(filename==NULL)
                filename = oldargv;
            else
                filename ++;
            
            p[0] = 0; sprintf(newargv,"/tmp/%s.tmp%d.%s",filename,pid,"tcl");
            p[0] = '.';
            
            INFO_Printf("replace: %s -> %s\n",oldargv,newargv);            

            strcpy(argv[i],newargv);

            /* create new file */
            oldfile = fopen(oldargv,"r");
            newfile = fopen(newargv,"w");

            if(oldfile==NULL) ERROR("old input file open failure");
            if(newfile==NULL) ERROR("new input file open failure");
            fprintf(newfile,"# -*-shell-script-*-\n");
            fprintf(newfile,"#TCL script file created from %s\n\n",oldargv);
            fprintf(newfile,"MD++ {\n");
            for(;;)
            {
                c=getc(oldfile);
                if(c==EOF) break;
                putc(c,newfile);
            }
            fprintf(newfile,"\n}\n");
            fclose(oldfile);
            fclose(newfile);
        }
    }
    INFO_Printf("\n");
    
    Tcl_Main(argc, argv, Tcl_AppInit);

    /* remove tmp files */
    sprintf(cmd,"rm -f /tmp/*.tmp%d.tcl\n",pid);
    //INFO_Printf(cmd);
    system(cmd);

    return TCL_OK;
}
#endif





/********************************************************************/
/* Coordinate transformation */
/********************************************************************/

void KMCFrame::SHtoR()
{ for(int i=0;i<_NP;i++) _R[i]=_H*_SR[i]; }
void KMCFrame::RHtoS()
{ Matrix33 hinv=_H.inv();
 for(int i=0;i<_NP;i++) _SR[i]=hinv*_R[i]; }
void KMCFrame::RtoR0()
{ int i; SHtoR(); for(i=0;i<_NP;i++) _R0[i]=_R[i]; }
void KMCFrame::R0toR()
{ int i; for(i=0;i<_NP;i++) _R[i]=_R0[i]; }

bool KMCFrame::Bond(int I, int J) const
{/* Returns true if bond(i,j) is stored in I */
    return (I<=J)? ((I^J)&1)  : !((I^J)&1);
}

/* Memory Allocation */
void KMCFrame::Alloc()
{
    int size;
    size=_NP*allocmultiple;
    
    Realloc(_SR,Vector3,size);
    Realloc(fixed,int,size);
    Realloc(_R,Vector3,size);
    Realloc(_R0,Vector3,size);
    Realloc(species,int,size);
    Realloc(group,int,size);
    Realloc(_EPOT_IND,double,size);

    memset(fixed,0,sizeof(int)*size);
    memset(species,0,sizeof(int)*size);
    memset(group,0,sizeof(int)*size);
    memset(_EPOT_IND,0,sizeof(double)*size);

    bindvar("SR",_SR,DOUBLE);
    bindvar("R",_R,DOUBLE);
    bindvar("R0",_R0,DOUBLE);
    bindvar("fixed",fixed,INT);
    bindvar("species",species,INT);
    bindvar("group",group,INT);
    bindvar("EPOT_IND",_EPOT_IND,DOUBLE);
    
}




/********************************************************************/
/* Kinetic Monte Carlo */
/********************************************************************/

void KMCFrame::potential()
{
    int_potential();
    ext_potential();
}

void KMCFrame::ext_potential()
{
    INFO("KMCFrame::ext_potential() not implemented yet.");
}

double KMCFrame::ext_potential(int iatom)
{
    INFO("KMCFrame::ext_potential() not implemented yet.");
    return _EPOT_IND[iatom];
}

void KMCFrame::int_potential()
{
    INFO("KMCFrame::int_potential() not implemented yet.");
}


double KMCFrame::int_potential(int iatom)
{
    INFO("KMCFrame::int_potential() not implemented yet.");
    return _EPOT;
}

double KMCFrame::int_potential(int iatom, int jatom)
{
    INFO("KMCFrame::int_potential() not implemented yet.");
    return _EPOT;
}

void KMCFrame::runkmc()
{
    /* Multi Process will execute this function simultaneously */

    for(curstep=0;curstep<totalsteps;curstep++)
    {
        while(win!=NULL)
        {
            if(win->IsPaused()) sleep(1);
            else break;
        }

        if(curstep%savepropfreq==0)INFO("curstep="<<curstep);
        winplot();
        
        if(curstep>=equilsteps)
        {
            if(saveprop)
            {
                if((curstep%savepropfreq)==0)
                {
                    int_potential();/* added by Keonwook Kang, Aug 14, 2006 */
                    calcprop(); /* added by Keonwook Kang, Aug 14, 2006 */
                    calcoutput();
                    pf.write(this);
                }
            }

            if(savecn)
                if((curstep%savecnfreq)==0&&(curstep!=0))
                    intercn.write(this);
        }

        kmc_step();
        
        if(win!=NULL)
            if(win->IsAlive())
            if(autowritegiffreq>0)
                if((curstep%autowritegiffreq)==0)
                {
                    win->LockWritegif(); //lock and wait for writegif
                }
    }
    if(saveprop)
    {
        int_potential();/* added by Keonwook Kang, Aug 14, 2006 */
        calcprop(); /* added by Keonwook Kang, Aug 14, 2006 */
        pf.write(this);
    }
}

void KMCFrame::calcprop()  /*calculate physical properties*/
{

}

void KMCFrame::calcoutput()
{
    char names[10000], *p, *q;
    /* separate output_fmt into short strings between white spaces */

    if(strlen(output_fmt)==0)
    {
        output[0]=0;
        return;
    }
    
    strcpy(names,output_fmt);
    p=names;

    output[0]=0;
    while(strlen(p)>0)
    {
        switch(*p)
        {
        case ' ':
        case '\t': /* white spaces */
            p++; break;
        default:
            q=strchr(p,' ');
            if(q!=NULL) *q=0;
            identify(p);
            if(curn>=0)
                switch(vartype[curn])                                    
                {                                                            
                case(INT): sprintf(output+strlen(output),"%10d ",*((int *)varptr[curn]+shift)); break; 
                case(LONG): sprintf(output+strlen(output),"%10ld ",*((long *)varptr[curn]+shift)); break;
                case(DOUBLE): sprintf(output+strlen(output),"%23.16e ",*((double *)varptr[curn]+shift)); break;
                case(STRING): sprintf(output+strlen(output),"%s ",(char *)varptr[curn]+shift); break;
                default: FATAL("unknown vartype ("<<vartype[curn]<<")");
                }
            if(q!=NULL) p=q+1;
            else goto it;
            break;
        }
    }
 it:
    sprintf(output+strlen(output),"%s","\n");
    //INFO_Printf("%s",output); return;
}

    
void KMCFrame::eval()
{
    /* multi process functions */
    potential();
    calcprop();
    calcoutput();
    INFO_Printf("%s\n",output_fmt);
    INFO_Printf("%s\n",output);
}
    
void KMCFrame::multieval()
{   /* for timing purposes */
    int i;
    for(i=0;i<totalsteps;i++)
        potential();
}










/********************************************************************/
/* Kinetic Monte Carlo Simulation */
/********************************************************************/

int KMCFrame::kmc_step()
{
    Vector3 dr, ds, dr0, ds0;
    Matrix33 hinv;
    double Eatom0,Eatom1;
    
    hinv=_H.inv();

    /* randomly select trial atom */
    select_atom = (int)floor(drand48()*_NP);
    while(fixed[select_atom]==-1) {select_atom = (int)floor(drand48()*_NP);}
        
    Eatom0=int_potential(select_atom);

    /* displace the atom, select_dir = +1,-1,+2,-2,+3,-3 */
    /* _SR[select_atom].[abs(select_dir)-1]+=(select_dir>0?1:-1); */

    Eatom1=int_potential(select_atom);

    /* update event list */

    /* select event from list */


    /* execute event */
    
    return select_atom;
}





































/********************************************************************/
/* Configuration manipulations */
/********************************************************************/

/* basis data for different crystal structures */
#include "lattice.h"

void KMCFrame::makecrystal()
{   /* create perfect crystal structure */
    UnitCell my_unitcell;
    Matrix33 myorient, lattice;
    Vector3 tmpv;
    int i,j,k,m,n;

    INFO("makecrystal");
    myorient.set(latticesize[0][0],latticesize[0][1],latticesize[0][2],
                 latticesize[1][0],latticesize[1][1],latticesize[1][2],
                 latticesize[2][0],latticesize[2][1],latticesize[2][2]);

    if(strcmp(crystalstructure,"simple-cubic")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(1,sc_basis);   /* simple-cubic */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"body-centered-cubic")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(2,bcc_basis); /* body-centered-cubic */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"face-centered-cubic")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,fcc_basis); /* face-centered-cubic */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"NaCl")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,nacl_basis,nacl_species); /* NaCl (B1) */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"CsCl")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(2,bcc_basis,cscl_species); /* NaCl (B1) */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"L1_2")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,fcc_basis,l12_species);   /* L1_2 */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"L1_0")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,fcc_basis,l10_species);   /* L1_0 */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"diamond-cubic")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,dc_basis);   /* diamond-cubic */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"beta-tin")==0)
    {
        if(latticeconst[1]<=0) latticeconst[1]=latticeconst[0];
        if(latticeconst[2]<=0) latticeconst[2]=latticeconst[0];
        //latticeconst[2]=latticeconst[0]*sqrt(2); /* when c=a*sqrt(2) is equivalent to diamond-cubic */
        //latticeconst[2]=latticeconst[0]*0.6;     /* when c=a*0.6 for beta-tin phase of Silicon */        
        my_unitcell.set(4,betatin_basis);   /* beta-tin */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"beta-tin-2")==0)
    {
        if(latticeconst[1]<=0) latticeconst[1]=latticeconst[0];
        if(latticeconst[2]<=0) latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,betatin_basis,betatin_2_species);   /* beta-tin-2 */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"zinc-blende")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,dc_basis,zb_species);   /* zinc-blende */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"CaF2")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(12,caf2_basis,caf2_species);   /* zinc-blende */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"hexagonal-ortho")==0)
    {
        latticeconst[1]=latticeconst[0]*1.73205080756887729352;
        my_unitcell.set(4,hex_ortho_basis);   /* hexagonal-ortho */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"rhombohedral-Si-Ge")==0)
    {
        latticeconst[1]=latticeconst[0]*1.73205080756887729352;
        latticeconst[2]=latticeconst[0]*4.89897948556636;
        my_unitcell.set(24,rhsige_basis,rhsige_species);   /* rhomohedral Si_Ge */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"ice-Ih-O")==0)
    {
        latticeconst[1]=latticeconst[0]*1.73205080756887729352;
        my_unitcell.set(8,ice_Ih_O_basis);   /* ice-Ih-Oxygen */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"ice-Ih")==0)
    {
        latticeconst[1]=latticeconst[0]*1.73205080756887729352;
        my_unitcell.set(24,ice_Ih_basis,ice_Ih_species);   /* ice-Ih */
        my_unitcell=my_unitcell*myorient.tran();
        nspecies=2;
    }
    else
    {
        ERROR("unknown lattice structure :"<<crystalstructure
              <<" for makecn");
    }
    INFO("my_unitcell\n"<<my_unitcell);

    n=0;
    _NP=(int)(latticesize[0][3]*latticesize[1][3]
              *latticesize[2][3]*my_unitcell.n);
    INFO("makecn: _NP="<<_NP);
    
    /* Multi process function */
    Alloc();
    
    for(i=0;i<(int)latticesize[0][3];i++)
        for(j=0;j<(int)latticesize[1][3];j++)
            for(k=0;k<(int)latticesize[2][3];k++)
            {
                for(m=0;m<my_unitcell.n;m++)
                {
                    tmpv.set(i,j,k);
                    _SR[n+m]=my_unitcell.basis[m]+tmpv;
                    species[n+m]=my_unitcell.species[m];
                }
                n+=my_unitcell.n;
            }
    for(i=0;i<_NP;i++)
    {
        _SR[i].x/=latticesize[0][3];
        _SR[i].y/=latticesize[1][3];
        _SR[i].z/=latticesize[2][3];
    }
    for(i=0;i<_NP;i++)
    {
        _SR[i]-=0.5;
        fixed[i]=0;
    }

    lattice.set(latticeconst[0],0.,0., 0.,latticeconst[1],0., 0.,0.,latticeconst[2]);
    _H=lattice*myorient.tran();
    lattice.set(latticesize[0][3],0.,0., 0.,latticesize[1][3],0., 0.,0.,latticesize[2][3]);
    _H=_H*lattice;
    
    _H=_H.reorient();

    for(i=0;i<_NP;i++)
    {
      _R[i] = _H*_SR[i];
    }
    
}



void KMCFrame::scaleH()
{
    /* multiply H matrix by input[0] */
    INFO(_H);
     _H*=input[0];
    INFO(_H);
}
void KMCFrame::saveH()
{
    /* copy H to H0 */
    INFO(_H);
    _H0=_H;
    INFO(_H);
}
void KMCFrame::setH()
{
    /* input = [ i j c ]
       H[i][j] = c
    */
    int i,j;
    double h;
    i=(int)input[0];
    j=(int)input[1];
    h=input[2];
    _H[i][j]=h;
}
void KMCFrame::restoreH()
{
    /* copy H0 to H */
    INFO(_H);
    _H=_H0;
    INFO(_H);
}






void KMCFrame::shiftbox()
{ /*
   * input = [ i j delta ]
   *
   *   H[:][i--] += H[:][j--]*delta
   *   while keeping s fixed (affine transformation)
   */   
    int xind, yind;
    double frac;
    if(input[0]==0)
    {
        ERROR("input(): no shift geometry set");
        return;
    }
    xind=(int)input[0]; //1x/2y/3z
    yind=(int)input[1];
    xind--;yind--;
    frac=input[2];
    
    INFO("The old H=\n"<<_H);
    
    _H[0][xind]+=_H[0][yind]*frac;
    _H[1][xind]+=_H[1][yind]*frac;
    _H[2][xind]+=_H[2][yind]*frac;
    
    INFO("The new H=\n"<<_H);
    SHtoR();
}

void KMCFrame::redefinepbc()
{ /*
   * input = [ i j delta ]
   *
   *   H[:][i--] += H[:][j--]*delta
   *   while keeping r fixed (introduce displacement at box boundary)
   */   
    int xind, yind;
    double frac;
    if(input[0]==0)
    {
        ERROR("input(): no shift geometry set");
        return;
    }
    xind=(int)input[0]; //1x/2y/3z
    yind=(int)input[1];
    xind--;yind--;
    frac=input[2];
    
    INFO("The old H=\n"<<_H);

    SHtoR();
    
    _H[0][xind]+=_H[0][yind]*frac;
    _H[1][xind]+=_H[1][yind]*frac;
    _H[2][xind]+=_H[2][yind]*frac;
    
    INFO("The new H=\n"<<_H);
    RHtoS();
}


void KMCFrame::extendbox()
{ /*
   * input = [ i n ]
   * 
   * The i-th direction (i=1,2,3) of the box will be multiplied by 
   *  n times (by inserting more atoms)
   * Usually allocmultiple = n should be set before
   *  makecrystal or readcn is called before this function
   */
    int dir, n;
    int i,j;

    dir=(int)input[0]; n=(int)input[1];
    if((dir<=0)||(dir>=4)) return;
    if(n<=1) return;
    dir--;

    for(i=0;i<_NP;i++)
    {
        for(j=1;j<n;j++)
        {
            _SR[i+_NP*j]=_SR[i];
            _SR[i+_NP*j][dir]+=j;
            fixed[i+_NP*j]=fixed[i];
            _EPOT_IND[i+_NP*j]=_EPOT_IND[i];
        }
    }
    for(i=0;i<_NP*n;i++)
    {
        _SR[i][dir]=(_SR[i][dir]+0.5)/n-0.5;
    }
    _H[0][dir]*=n;
    _H[1][dir]*=n;
    _H[2][dir]*=n;
    _NP*=n;
    SHtoR();
    INFO("extenbox: new NP="<<_NP);
}




void KMCFrame::moveatom()
{ /*
   * Move specified atoms
   *
   * input = [ n dx dy dz i1 i2 ... in ]
   *
   * n:            number of atoms to be moved (n>0)
   * (dx,dy,dz):   displacement vector in real space
   * i1,i2,...,in: indices of atoms to be moved
   *
   */
    int i, n, id, fix;
    double mx, my, mz;
    double x0, y0, z0, r, r2, dmin, dmax;
    Vector3 s0, ds, dr;

    n=(int)input[0];

    INFO_Printf("n=%d\n",n);

    if(n>0)
    {
        mx=input[1];
        my=input[2];
        mz=input[3];
        
        SHtoR();
        for(i=0;i<n;i++)
        {
            id=(int)input[4+i];
            _R[id].x+=mx;
            _R[id].y+=my;
            _R[id].z+=mz;
        }
        RHtoS();
    }
    else if(n==-2)
    {
        /* see above */
    }
    if(n==0)
    {
        INFO("move atoms in circular plate");
        x0=input[1];
        y0=input[2];
        z0=input[3];
        r =input[4];
        dmin=input[5];
        dmax=input[6];
        mx=input[7];
        my=input[8];
        mz=input[9];
        fix=(int)input[10];
        
        r2=r*r;
        s0.set(x0,y0,z0);
        for(i=0;i<_NP;i++)
        {            
            if((_SR[i].y>=dmax)||(_SR[i].y<dmin)) continue;
            if(fix)
            {
                fixed[i]=1;
                INFO_Printf("fix atom %d in plane\n",i);
            }
            ds=_SR[i]-s0;
            ds.subint();
            
            dr=_H*ds;
            if(dr.norm2()>r2) continue;
            
            _SR[i].x+=mx;
            _SR[i].y+=my;
            _SR[i].z+=mz;
            INFO_Printf("move atom %d\n",i);
        }        
        /* shift box */
        SHtoR();        
    }
    return;
}

void KMCFrame::movegroup()
{ /*
   * Move specified groups of atoms
   *
   * input = [ n dx dy dz grp1 grp2 ... grpn ]
   *
   * n:            number of groups to be moved (n>0)
   * (dx,dy,dz):   displacement vector in real space
   * grp1,grp2,...,grpn: indices of group to be moved
   *
   */
    int i, n, id, ip;
    double mx, my, mz;

    n=(int)input[0];

    INFO_Printf("n=%d\n",n);

    if(n>0)
    {
        SHtoR();
        mx=input[1];
        my=input[2];
        mz=input[3];
        
        SHtoR();
        for(i=0;i<n;i++)
        {
            id=(int)input[4+i];
            for(ip=0;ip<_NP;ip++)
            {
                if(group[ip]==id)
                {
                    _R[ip].x+=mx;
                    _R[ip].y+=my;
                    _R[ip].z+=mz;
                }
            }
        }
        RHtoS();
    }
    return;
}


void KMCFrame::pbcshiftatom()
{
    int i;
    double dx,dy,dz;
    dx=input[0];
    dy=input[1];
    dz=input[2];

    for(i=0;i<_NP;i++)
    {
        _SR[i].x+=dx;_SR[i].y+=dy;_SR[i].z+=dz;
        _SR[i].subint();
    }
    SHtoR();
}




void KMCFrame::fixatoms_by_ID()
{ /* fix a set of atoms whose ID are specified
   *
   * input = [ n i1 i2 ... in ]
   *
   * n:             number of atoms to be fixed
   * i1,i2,...,in:  indices of atoms to be fixed
   */
    int i, n;
    n=(int)input[0];
//    for(i=0;i<_NP;i++)
//    {
//        fixed[i]=0;
//    }
    for(i=1;i<=n;i++)
    {
        fixed[(int)input[i]]=1;
    }
}

void KMCFrame::fixatoms_by_position()
{ /* fix all atoms whose position falls within a specified regime
   *
   * input = [ enable x0 x1 y0 y1 z0 z1 ]
   *
   * enable: set to 1 to activate
   * atoms whose reduced coordinate (x,y,z) satisfie
   *   x0 <= x < x1
   *   y0 <= y < y1
   *   z0 <= z < z1
   * will be fixed
   */
    int i, n, nfixed;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    n=(int)input[0];
    xmin=input[1];
    xmax=input[2];
    ymin=input[3];
    ymax=input[4];
    zmin=input[5];
    zmax=input[6];

    INFO_Printf("setfixedatoms = [\n");
    nfixed=0;
    for(i=0;i<_NP;i++)
    {
        if((_SR[i].x>=xmin)&&(_SR[i].x<xmax)
           &&(_SR[i].y>=ymin)&&(_SR[i].y<ymax)
           &&(_SR[i].z>=zmin)&&(_SR[i].z<zmax))
        {
            fixed[i]=1;
            nfixed++;
        }
    }
    INFO_Printf("%d\n",nfixed);
    for(i=0;i<_NP;i++)
    {
        if((_SR[i].x>=xmin)&&(_SR[i].x<xmax)
           &&(_SR[i].y>=ymin)&&(_SR[i].y<ymax)
           &&(_SR[i].z>=zmin)&&(_SR[i].z<zmax))
        {
            INFO_Printf("%d ",i);
        }
    }
    INFO_Printf("\n]\n");
    INFO_Printf("%d atoms fixed\n",nfixed);
}


void KMCFrame::fixallatoms()
{
    int i;

    for(i=0;i<_NP;i++)
    {
        fixed[i]=1;
    }
}

void KMCFrame::freeallatoms()
{
    int i;

    for(i=0;i<_NP;i++)
    {
        fixed[i]=0;
    }
}

void KMCFrame::reversefixedatoms()
{
    int i;

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==0) fixed[i]=1;
        else if(fixed[i]==1) fixed[i]=0;
    }
}


void KMCFrame::setfixedatomsspecies()
{ /* set the species value of all fixed atoms to input[0]
   */
    int i;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
            species[i]=(int)input[0];
    }
}
       
void KMCFrame::setfixedatomsgroup()
{ /* set the group ID of all fixed atoms to input[0]
   */
    int i;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
            group[i]=(int)input[0];
    }
}
       
void KMCFrame::reversespecies()
{ /* set the species value of all fixed atoms to input[0]
   */
    int i;
    for(i=0;i<_NP;i++)
    {
        if(species[i]==0)
            species[i]=1;
        else if(species[i]==1)
            species[i]=0;        
    }
}
       
void KMCFrame::removefixedatoms()
{ /* remove all atom i if fixed[i] is not equal to zero */
    int i, j;
    
    j=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==0)
        {
            _SR[j]=_SR[i];
            fixed[j]=fixed[i];
            species[j]=species[i];
            _EPOT_IND[j]=_EPOT_IND[i];
            j++;
        }
    }
    _NP=j;
    SHtoR();
}

void KMCFrame::markremovefixedatoms()
{ /* if fixed[i]==1, change it to fixed[i]=-1 */
    int i;
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
        {
            fixed[i]=-1;
        }
    }
}

void KMCFrame::movefixedatoms()
{ /* remove all atom i if fixed[i] is not equal to zero */
    int i;
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
        {
            _SR[i].x+=input[0];
            _SR[i].y+=input[1];
            _SR[i].z+=input[2];
        }
    }
}

void KMCFrame::removeellipsoid()
{
    /* remove all atoms within an ellipsoid
     *
     * input = [ enable x0 y0 z0 a b c ]
     *
     * enable:     if 1 then removeellipsoid is enabled
     * (x0,y0,z0): center of the ellipsoid
     * a,b,c:      semi-axis of the ellipsoid
     */
    
    int i, n, removenum;
    double x0, y0, z0, arem, brem, crem;
    Vector3 r0, dr;
    
    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("removeellipsoid(): no geometry set");
        return;
    }

    /* read in parameters */
    x0=input[1];
    y0=input[2];
    z0=input[3];
    arem=input[4];
    brem=input[5];
    crem=input[6];
    
    r0.set(x0,y0,z0);
    /* remove all atoms outside rrem */
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        dr = _R[i]-r0;
        dr.x/=arem;
        dr.y/=brem;
        dr.z/=crem;
        if(dr.norm2()<=1)
            fixed[i]=-1;
    }
    
    /* commit removing atoms */
    n=0; SHtoR();
    removenum=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]!=-1)
        {
            _R0[n]=_R[i];
            fixed[n]=fixed[i];
            n++;
        }
        else
        {
            removenum++;
        }
    }
    INFO("originally "<<_NP<<" atoms");
    INFO("remove "<<removenum<<" atoms");
    _NP-=removenum;
    INFO("now "<<_NP<<" atoms");
    R0toR();  RHtoS();
    INFO("NP="<<_NP<<"  n="<<n);

}

void KMCFrame::removerectbox()
{ /* Remove a rectangular box of atoms
   *
   * input = [ ax ay az La bx by bz Lb cx cy cz Lc x0 y0 z0 plotonly outside ]
   *
   */
    int i, j, assist, plotonly, outside;
    Vector3 ds, dss, sa, sb, sc, s0;
    Vector3 ha, hb, hc;
    Matrix33 sh, shinv, latt, lattinv;
    

    ha.set(latticesize[0][0],latticesize[0][1],latticesize[0][2]);
    hb.set(latticesize[1][0],latticesize[1][1],latticesize[1][2]);
    hc.set(latticesize[2][0],latticesize[2][1],latticesize[2][2]);
    ha*=latticesize[0][3];
    hb*=latticesize[1][3];
    hc*=latticesize[2][3];
    latt.setcol(ha,hb,hc);
    lattinv=latt.inv();
        
    sa.set(input[0],input[1],input[2]);
    sb.set(input[4],input[5],input[6]);
    sc.set(input[8],input[9],input[10]);
    sa*=input[3];
    sb*=input[7];
    sc*=input[11];
    s0.set(input[12],input[13],input[14]);
    assist=(int)input[15];
    plotonly=(int)input[16];
    outside=(int)input[17];

    sa=lattinv*sa;
    sb=lattinv*sb;
    sc=lattinv*sc;
    s0=lattinv*s0;
    
    INFO("assist="<<assist);
    INFO("plotonly="<<plotonly);
    
    sh.setcol(sa,sb,sc);
    shinv=sh.inv();
    INFO("vol="<<sh.det());
    INFO("s0="<<s0);
    
    for(i=0;i<_NP;i++)
    {
        ds=_SR[i]-s0;
        /*ds.subint();*/
        dss=shinv*ds;
        if(assist&&(!plotonly))
        {
            if((dss.x>=-0.5)&&(dss.x<0.5)&&
               (dss.y>=-0.5)&&(dss.y<0.5))
            {
                if(fabs(dss.z)<2)
                {
                    if((dss.z>=-0.5)&&(dss.z<0.5))
                    {
                        fixed[i]=-1;/* mark for deletion */
                    }
                    else if(dss.z>0.5)
                    {
                        dss.z=(dss.z-0.5)/(2-0.5)*(2-0.1)+0.1-dss.z;
                        dss.x=0; dss.y=0;
                        INFO("dss="<<dss);
                        _SR[i]+=sh*dss;
                    }
                    else
                    {
                        dss.z=(dss.z+0.5)/(-2+0.5)*(-2+0.1)-0.1-dss.z;
                        dss.x=0; dss.y=0;
                        INFO("dss="<<dss);
                        _SR[i]+=sh*dss;
                    }
                }
            }
        }
        else
        {
            if((dss.x>=-0.5)&&(dss.x<0.5)&&
               (dss.y>=-0.5)&&(dss.y<0.5)&&
               (dss.z>=-0.5)&&(dss.z<0.5))
            {
                fixed[i]=-1;/* mark for deletion */
            }
        }
    }
    if(outside)
    {
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]==-1)
            {
                fixed[i]=0;
            }
            else if(fixed[i]==0)
            {
                fixed[i]=-1;
            }
        }
    }
    if(!plotonly)
    {
        j=0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]!=-1)
            {
                _SR[j]=_SR[i];
                fixed[j]=0;
                _EPOT_IND[j]=_EPOT_IND[i];
                j++;
            }
        }
        INFO("removerectbox: "<<_NP-j<<" of atoms removed");
        _NP=j;
        INFO("new number of atoms "<<_NP);
        SHtoR();
    }
}

void KMCFrame::randomposition()
{
    int i;
    for(i=0;i<_NP;i++)
    {
        _SR[i].x=drand48()-0.5;
        _SR[i].y=drand48()-0.5;
        _SR[i].z=drand48()-0.5;
    }
}




/********************************************************************/
/* Neighbor list */
/********************************************************************/
void KMCFrame::NbrList_reconstruct()
{
    /* use cell list combined with Verlet list */
    int ncx,ncy,ncz;
    int ncx0,ncx1,ncy0,ncy1,ncz0,ncz1,np0,np1;
    int i,j,k,n,i1,j1,k1,n1,i2,j2,k2;
    int ipt,jpt, itmp;
    double rsmax2;
    Vector3 s,sij,rij,h;
    int maxnic, maxnnm, m_nic, m_nnm;

    NbrList_init(_NP,_NNM);
    
    rsmax2=_RLIST*_RLIST;
    h=_H.height();

    INFO("reconstruct neighborlist");
    DUMP("reconstruct_cell()  ");
    INFO("H="<<_H<<"  _RLIST="<<_RLIST);
    
    ncx=(int)floor(h[0]/((_RLIST)*1.05));if(ncx==0)ncx=1;
    ncy=(int)floor(h[1]/((_RLIST)*1.05));if(ncy==0)ncy=1;
    ncz=(int)floor(h[2]/((_RLIST)*1.05));if(ncz==0)ncz=1;
    
    NbrList_initcell(ncx,ncy,ncz,_NIC);
    
    /* Load distribution */
    if((ncx>=ncy)&&(ncx>=ncz))
    {
        ncx0=0;
        ncx1=ncx;
        ncy0=0; ncy1=ncy; 
        ncz0=0; ncz1=ncz;
    }
    else if((ncy>=ncx)&&(ncy>=ncz))
    {
        ncy0=0;
        ncy1=ncy;
        ncx0=0; ncx1=ncx; 
        ncz0=0; ncz1=ncz;
    }
    else 
    {
        ncz0=0;
        ncz1=ncz;
        ncx0=0; ncx1=ncx; 
        ncy0=0; ncy1=ncy;
    }

    np0=0;
    np1=_NP;

    for(i=ncx0;i<ncx1;i++)
        for(j=ncy0;j<ncy1;j++)
            for(k=ncz0;k<ncz1;k++)
                celllist[i][j][k][0]=0;
    
    m_nic = 0;
    for(ipt=0;ipt<_NP;ipt++)
    {
        s=_SR[ipt];
        s.subint();
        i=(int)floor((s.x+0.5)*ncx);i=(i+2*ncx)%ncx;
        j=(int)floor((s.y+0.5)*ncy);j=(j+2*ncy)%ncy;
        k=(int)floor((s.z+0.5)*ncz);k=(k+2*ncz)%ncz;
        if((i>=ncx0)&&(i<ncx1)&&(j>=ncy0)&&(j<ncy1)&&(k>=ncz0)&&(k<ncz1))
        {
            //INFO_Printf("NbrList_reconstruct: ipt=%d i=%d j=%d k=%d\n",ipt,i,j,k);
            celllist[i][j][k][0]++;
            if(celllist[i][j][k][0]>_NIC-1)
                FATAL("reconstruct(): too many atoms per cell "
                      <<celllist[i][j][k][0]<<" > limit("<<(_NIC-1)<<")"
                      <<"increase NIC in script");
            celllist[i][j][k][celllist[i][j][k][0]]=ipt;
            if(m_nic<celllist[i][j][k][0]) m_nic=celllist[i][j][k][0];
        }
    }
    maxnic=m_nic;
    
    DUMP("celllist constructed");

    for(ipt=np0;ipt<np1;ipt++)
    {
        nn[ipt]=0;
        for(itmp=0;itmp<_NNM;itmp++)
            nindex[ipt][itmp]=0;
    }
    DUMP("list cleared");

    m_nnm=0;
    for(i=ncx0;i<ncx1;i++)
        for(j=ncy0;j<ncy1;j++)
            for(k=ncz0;k<ncz1;k++)
                for(n=1;n<=celllist[i][j][k][0];n++)
                {
                    ipt=celllist[i][j][k][n];
                    //neighboring cell
                    for(i1=((ncx==1)?i:(i-1));i1<=((ncx<=2)?i:(i+1));i1++)
                        for(j1=((ncy==1)?j:(j-1));j1<=((ncy<=2)?j:(j+1));j1++)
                            for(k1=((ncz==1)?k:(k-1));k1<=((ncz<=2)?k:(k+1));k1++)
                            {
                                i2=(i1+2*ncx)%ncx;
                                j2=(j1+2*ncy)%ncy;
                                k2=(k1+2*ncz)%ncz;
                                //INFO_Printf("i2=%d j2=%d k2=%d\n",i2,j2,k2);
                                for(n1=1;n1<=celllist[i2][j2][k2][0];n1++)
                                {
                                    jpt=celllist[i2][j2][k2][n1];
                                    if(ipt==jpt) continue;
                                    sij=_SR[ipt]-_SR[jpt];
                                    sij.subint(); // subtract the nearest integer
                                    rij=_H*sij;
                                    if(rij.norm2()<rsmax2)
                                    {
                                        nindex[ipt][nn[ipt]]=jpt;
                                        nn[ipt]++;
                                        if(nn[ipt]>=_NNM)
                                            FATAL("reconstruct: index["<<ipt<<"] ("
                                                  <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                                  "increase NNM in script");
                                        if(m_nnm<nn[ipt]) m_nnm=nn[ipt];
                                    }
                                }
                            }
                }
    
    maxnnm=m_nnm;
    
    /* combine local lists to shm list */
    for(i=np0;i<np1;i++)
    {
        _R[i]=_H*_SR[i];
        _R0[i]=_R[i];
    }

    INFO_Printf("reconstruct() finished. maxnic=%d maxnnm=%d\n",maxnic,maxnnm);
}

void KMCFrame::NbrList_print() 
{   /* debug use */
    int i;
    fprintf(stderr,"rlist=%18.10f skin=%18.10f\n", _RLIST, _SKIN);
    for(i=0;i<10;i++)
    {
        fprintf(stderr,"nn(%d)=%d\n",i,nn[i]);
    }
    for(i=_NP-10;i<_NP;i++)
    {
        fprintf(stderr,"nn(%d)=%d\n",i,nn[i]);
    }
}

bool KMCFrame::NbrList_needrefresh()
{
    /* all processes follow same procedure */
    double maxd;
    int i, need;
    if (firsttime)
    {
        firsttime=false;
        return true;
    }
    maxd=0;
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        maxd=max(fabs(_R[i].x-_R0[i].x),maxd);
        maxd=max(fabs(_R[i].y-_R0[i].y),maxd);
        maxd=max(fabs(_R[i].z-_R0[i].z),maxd);
    }
    need=(maxd>_SKIN);

    return (need!=0);
}

void KMCFrame::NbrList_init(int mx,int mz)
{
    /* Multi process function */
    /* init two dimensional array */
    int shft1,shft2;//, mx, mz;
    int i;

    DUMP("initlist("<<mx<<","<<mz<<")");
    shft1=mx*mz*sizeof(int);
    shft2=mx*sizeof(int *);
    if(shft1+shft2==0) return;
    
    Realloc(nindex_mem,char,(shft1+shft2));

    nindex=(int **)(nindex_mem+shft1);

    memset(nindex_mem,0,shft1+shft2);
    for(i=0;i<mx;i++)
    {
        nindex[i]=(int *)(nindex_mem+i*mz*sizeof(int));
    }

    Realloc(nn,int,mx);
}

void KMCFrame::NbrList_initcell(int mx, int my, int mz, int mc)
{
    /* init two dimensional array */
    int shft1,shft2,shft3,shft4;
    int i, j, k;

//    INFO_Printf("MDFrame::NbrList_initcell(%d,%d,%d,%d)\n",mx,my,mz,mc);
//    freecelllist();
    shft1=mx*my*mz*mc*sizeof(int);
    shft2=mx*my*mz*sizeof(int *);
    shft3=mx*my*sizeof(int **);
    shft4=mx*sizeof(int ***);
    if(shft1+shft2+shft3+shft4==0) return;
    
    INFO("initcell("<<mx<<","<<my<<","<<mz<<","<<mc<<") cell list size "<<shft1+shft2+shft3+shft4);
    Realloc(cell_mem,char,(shft1+shft2+shft3+shft4));
    celllist=(int ****)(cell_mem+shft1+shft2+shft3);

    memset(cell_mem,0,shft1+shft2+shft3+shft4);
    for(i=0;i<mx;i++)
    {
        celllist[i]=(int ***)(cell_mem+shft1+shft2+i*my*sizeof(int **));
        for(j=0;j<my;j++)
        {
            celllist[i][j]=(int **)(cell_mem+shft1+(i*my*mz+j*mz)*sizeof(int *));
            for(k=0;k<mz;k++)
            {
                celllist[i][j][k]=
                    (int *)(cell_mem+(i*my*mz*mc+j*mz*mc+k*mc)*sizeof(int));
            }
        }
    }
//    INFO_Printf("MDFrame::NbrList_initcell done cell=%x\n",celllist);
//    INFO_Printf(" celllist[0][0][0][0] = %d\n",celllist[0][0][0][0]);
}

void KMCFrame::NbrList_free()
{
    nindex_mem=NULL;
    nn=NULL;
    nindex=NULL;
}

void KMCFrame::NbrList_freecell()
{
    free(cell_mem);
    cell_mem=NULL;
    celllist=NULL;
}

void KMCFrame::NbrList_refresh()
{
    if(NbrList_needrefresh()) NbrList_reconstruct();
}






/********************************************************************/
/* File input and output */
/********************************************************************/
    
int KMCFrame::readcn()
{
    initcn.open(incnfile,LFile::O_Read);
    return initcn.read(this);
}

int KMCFrame::setfilecounter()
{
    intercn.setcount(filecounter);
    return 0;
}

int KMCFrame::writefinalcnfile(int zip, bool bg)
{
    finalcn.open(finalcnfile);
    return finalcn.write(this,zip,bg);
}

void KMCFrame::saveintercn(int step)
{
    if(savecn)
        if((step%savecnfreq)==0&&(step!=0))
            writeintercnfile();
}

int KMCFrame::writeintercnfile(int zip,bool bg)
{
    return intercn.write(this,zip,bg);
}

int KMCFrame::openintercnfile()
{
    return intercn.open(intercnfile);
}

int KMCFrame::openpropfile()
{
    return pf.open(outpropfile);
}

/* Configuration Files */
char * CNFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","Configuration File for Atom Positions, Format:\n"
            "sx1 sy1 sz1; ... ; sxn syn szn; H(3x3)");
    return tmp;
}

int CNFile::writeblock(void *p)
{
    int i;
    KMCFrame &d=*((KMCFrame *)p);
    f->Printf("%d\n",d._NP);
    for(i=0;i<d._NP;i++)
    {
        f->Printf("%25.17E %25.17E %25.17E %2d %2d %2d\n",
                  d._SR[i].x,d._SR[i].y,d._SR[i].z,
                  d.fixed[i],d.species[i],d.group[i]);
    }
    f->Printf("%25.17E %25.17E %25.17E\n"
              "%25.17E %25.17E %25.17E\n"
              "%25.17E %25.17E %25.17E\n",
              d._H[0][0],d._H[0][1],d._H[0][2],
              d._H[1][0],d._H[1][1],d._H[1][2],
              d._H[2][0],d._H[2][1],d._H[2][2]);
    f->Printf("%d ",d.nspecies);
    for(i=0;i<d.nspecies;i++) f->Printf("%s ",d.element[i]);
    f->Printf("\n");
    return 0;
}

int CNFile::readblock(void *p)
{
    int i;
    char *buffer, *pp, *q;
    KMCFrame &d=*((KMCFrame *)p);
    
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    sscanf(pp, "%d", &d._NP);
    INFO("readblock: NP="<<d._NP);    
    
    d.Alloc();
    
    pp=strchr(pp, '\n');
    pp++;
    for(i=0;i<d._NP;i++)
    {
        q=pp;
        pp=strchr(pp,'\n');
        if(pp) *(char *)pp++=0;
        d.fixed[i]=0;
        sscanf(q, "%lf %lf %lf %d %d %d",
               &(d._SR[i].x),&(d._SR[i].y),&(d._SR[i].z),
               &(d.fixed[i]),&(d.species[i]),&(d.group[i]));
    }

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%lf %lf %lf",&d._H[0][0],&d._H[0][1],&d._H[0][2]);
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%lf %lf %lf",&d._H[1][0],&d._H[1][1],&d._H[1][2]);
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%lf %lf %lf",&d._H[2][0],&d._H[2][1],&d._H[2][2]);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d %s %s %s %s %s %s %s %s %s %s",
           &(d.nspecies),
           d.element[0],d.element[1],d.element[2],d.element[3],d.element[4],
           d.element[5],d.element[6],d.element[7],d.element[8],d.element[9]);
    Free(buffer);
    DUMP("readblock finished");
    return 0;
}


/* Output Property Files */
char * PropFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","property file at each time step");
    return tmp;
}

int PropFile::writeentry(void *p)
{
    char tmp[500];
    KMCFrame &d=*((KMCFrame *)p);
    if(strlen(d.output)>0)
    {
        *f<<d.output;
    }
    else
    {   /* default output format */
        sprintf(tmp,"%12d %20.13e %6d %6d\n",
                d.curstep,                  /* 1  */
                d._EPOT,                    /* 2  */
                d.select_atom,              /* 3 */
                d.select_dir                /* 4 */
            );
        *f<<tmp;
    }
    return 0;
}

















/* interface with Li Ju's AtomEye viewer */
void KMCFrame::writeatomeyecfgfile(char *fname)
{
    FILE *fp;
    int i, nplot, nr, n1, n2, n3, nn, ix, iy, iz, k;
    Matrix33 h;  Vector3 s;
    int nw, cind; double Emin, Emax;
    
    INFO("KMCFrame::writeatomeyecfgfile "<<fname);
    
    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writeatomeyecfgfile: open file failure");
    }

    nplot = 0;
    nr = atomeyerepeat[0];
    if(nr!=0)
    {
        n1=atomeyerepeat[1];
        n2=atomeyerepeat[2];
        n3=atomeyerepeat[3];
        if(n1<=0) n1=1;
        if(n2<=0) n2=1;
        if(n3<=0) n3=1;
        nn = n1*n2*n3;
        INFO_Printf("atomeyerepeat ( %d %d %d )\n",
                    n1, n2, n3);
    }
    else
    {
        nn = n1 = n2 = n3 = 1;
    }

    nw = (int)plot_color_windows[0];
    for(i=0;i<_NP;i++)
    {
        if(nw>0)
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    continue;
                }
                for(k=0;k<nw;k++)
                {
                    Emin = plot_color_windows[k*3+1];
                    Emax = plot_color_windows[k*3+2];
                    cind = (int) plot_color_windows[k*3+3];

                    if((color_ind[i]>=Emin)&&
                       (color_ind[i]<=Emax))
                    {
                        nplot ++;
                        continue;    
                    }
                }
        }
        else
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    
                    continue;
                }        
            nplot ++;
        }
    }
    h=_H.tran();
    fprintf(fp,"Number of particles = %d\n",nplot*nn);
    fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
    fprintf(fp,"H0(1,1) = %f A\n",h[0][0]*n1);
    fprintf(fp,"H0(1,2) = %f A\n",h[0][1]*n1);
    fprintf(fp,"H0(1,3) = %f A\n",h[0][2]*n1);
    fprintf(fp,"H0(2,1) = %f A\n",h[1][0]*n2);
    fprintf(fp,"H0(2,2) = %f A\n",h[1][1]*n2);
    fprintf(fp,"H0(2,3) = %f A\n",h[1][2]*n2);
    fprintf(fp,"H0(3,1) = %f A\n",h[2][0]*n3);
    fprintf(fp,"H0(3,2) = %f A\n",h[2][1]*n3);
    fprintf(fp,"H0(3,3) = %f A\n",h[2][2]*n3);
#define ATOMEYE_EXTCFG
#ifdef ATOMEYE_EXTCFG
    fprintf(fp,".NO_VELOCITY.\n");
    fprintf(fp,"entry_count = 4\n");
    fprintf(fp,"auxiliary[0] = pote [eV]\n");
//    fprintf(fp,"1.00000\n");
//    fprintf(fp,"%s\n",element[0]);
    for(i=0;i<_NP;i++)
    {
        s=_SR[i];
        if(plot_map_pbc==1) s.subint();
        if(plot_limits[0]==1)
            if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
               ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
               ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
            {
                continue;
            }
        
        if(nw>0)
        {
            for(k=0;k<nw;k++)
            {
                Emin = plot_color_windows[k*3+1];
                Emax = plot_color_windows[k*3+2];
                cind = (int) plot_color_windows[k*3+3];
                
                if((color_ind[i]>=Emin)&&
                   (color_ind[i]<=Emax))
                {
                    for(ix=0;ix<n1;ix++)
                        for(iy=0;iy<n2;iy++)
                            for(iz=0;iz<n3;iz++)
                            {
                                fprintf(fp,"%f\n%s\n %f %f %f %f  \n",
                                        _ATOMMASS[species[i]],element[species[i]],
                                        (s.x+0.5+ix)/n1,
                                        (s.y+0.5+iy)/n2,
                                        (s.z+0.5+iz)/n3,
                                        color_ind[i]);
                            }

                    continue;    
                }
            }
        }
        else
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    continue;
                }        
            for(ix=0;ix<n1;ix++)
                for(iy=0;iy<n2;iy++)
                    for(iz=0;iz<n3;iz++)
                    {
                        fprintf(fp,"%f\n%s\n%f %f %f %f  \n",
                                _ATOMMASS[species[i]],element[species[i]],                                
                                (s.x+0.5+ix)/n1,
                                (s.y+0.5+iy)/n2,
                                (s.z+0.5+iz)/n3,
                                color_ind[i]);
                    }
        }
    }        
#else
    for(i=0;i<_NP;i++)
    {
        for(ix=0;ix<n1;ix++)
            for(iy=0;iy<n2;iy++)
                for(iz=0;iz<n3;iz++)
                {
                    fprintf(fp,"%f %s %f %f %f %f %f %f\n",
                            _ATOMMASS[species[i]],element[species[i]],
                            (s.x+0.5+ix)/n1,
                            (s.y+0.5+iy)/n2,
                            (s.z+0.5+iz)/n3,
                            _VSR[i].x/n1,
                            _VSR[i].y/n2,
                            _VSR[i].z/n3);
                }
    }
#endif
    fclose(fp);
}

void KMCFrame::convertCNtoCFG()
{
    int i, istart, iend;
    char cfgfname[100];
    
    istart = (int) input [0];
    iend   = (int) input [1];

    for(i=istart;i<=iend;i++)
    {
        if(strlen(incnfile)>=12)
        {
            sprintf(incnfile + strlen(incnfile)-12,"inter%04d.cn",i);
            strcpy(cfgfname,incnfile);
            sprintf(cfgfname + strlen(cfgfname)-12,"inter%04d.cfg",i);
        }
        else
        {
            sprintf(incnfile,"inter%04d.cn",i);
            sprintf(cfgfname,"inter%04d.cfg",i);
        }
        readcn();
        writeatomeyecfgfile(cfgfname);
    }
        
}



void KMCFrame::atomeye()
{
    char extpath[200], comsav[200];

    LFile::SubHomeDir(atomeyepath,extpath);
    
    writeatomeyecfgfile("atomeye.cfg");
    strcpy(comsav,command);

#if 1 /* compress */
    sprintf(command,"gzip -f atomeye.cfg");
    runcommand();

    sprintf(command,"\"%s/%s\" atomeye.cfg.gz &",extpath,atomeyeexe);
    runcommand();
#else /* no compress */
    sprintf(command,"\"%s/%s\" atomeye.cfg &",extpath,atomeyeexe);
    runcommand();
#endif
    strcpy(command,comsav);
}




void KMCFrame::writeENERGY(char *fname)
{
    FILE *fp;
    int i;
    
    INFO("KMCFrame::writeENERGY "<<fname);

    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    default: ERROR("writeENERGY() unknown coloraxis "<<plot_color_axis); return;
    }
    
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writeENERGY: open file failure");
    }

    for(i=0;i<_NP;i++)
    {
        fprintf(fp,"  %20.16e\n",color_ind[i]);
    }

    fclose(fp);
}

void KMCFrame::writePOSITION(char *fname)
{
    FILE *fp;
    int i;
    
    INFO("KMCFrame::writePOSITION "<<fname);

    SHtoR();
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writePOSITION: open file failure");
    }

    if(input[0]==1) INFO("only write positions of fixed atoms");
    for(i=0;i<_NP;i++)
    {
        if(input[0]==1)
        {
            if(fixed[i]!=0)
                fprintf(fp,"%6d  %20.12e %20.12e %20.12e %4d\n",
                        i,_R[i].x,_R[i].y,_R[i].z,fixed[i]);
        }
        else
        {
            fprintf(fp,"%6d  %20.12e %20.12e %20.12e %4d\n",
                    i,_R[i].x,_R[i].y,_R[i].z,fixed[i]);
        }
    }
    fclose(fp);
}


void KMCFrame::GnuPlotHistogram()
{ /* plot histogram of atom properties (e.g. energy, topology)
   * by calling Gnuplot
   *
   * input = [ E0 E1 n ]
   *
   * histogram from E0 to E1 with n bins
   */
    
    FILE *fp;
    int i, nbin, *counts, n;
    double Emin, Emax, *Ebin, dE;
    char comsav[200];
    
    INFO("KMCFrame::GnuPlotEnergyHistogram");

    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    /* determining plot specification */
    Emin = input[0];
    Emax = input[1];
    nbin = (int)input[2];

    if(Emin==Emax)
    {
        Emin=Emax=color_ind[0];
        for(i=1;i<_NP;i++)
        {
            if(Emin>color_ind[i]) Emin=color_ind[i];
            if(Emax<color_ind[i]) Emax=color_ind[i];
        }
    }
    if(fabs(Emin-Emax)<2e-8)
    {
        Emin-=1e-8;
        Emax+=1e-8;
    }
    if(nbin<=0) nbin=20;
    dE=(Emax-Emin)/nbin;
    
    /* constructing energy histogram */
    counts=(int    *)malloc(sizeof(double)*nbin);
    Ebin  =(double *)malloc(sizeof(double)*nbin);
    for(i=0;i<nbin;i++)
    {
        Ebin[i]=dE*(i+0.5)+Emin;
        counts[i]=0;
    }
    
    for(i=0;i<_NP;i++)
    {
        n=(int)floor((color_ind[i]-Emin)/dE);
        if((n>=0)&&(n<nbin))
            counts[n]++;
    }
    /* write energy into file */
    writeENERGY("eng.out");
    /* write energy histogram data into file */
    fp=fopen("enghist.out","w");
    if(fp==NULL)
    {
        FATAL("GnuPlotEnergyHistogram: open file failure");
    }

    for(i=0;i<nbin;i++)
    {
        fprintf(fp,"  %20.16e %d\n",Ebin[i],counts[i]);
    }
    fclose(fp);

    free(counts);
    free(Ebin);

    /* launch GnuPlot */
    /* write energy histogram data into file */
    fp=fopen("plotenghist.gp","w");
    if(fp==NULL)
    {
        FATAL("GnuPlotEnergyHistogram: open file failure");
    }
    fprintf(fp,"plot 'enghist.out' with histeps\n"); 
    fclose(fp);
    
    fp=fopen("gpw","w");
    if(fp==NULL)
    {
        FATAL("GnuPlotEnergyHistogram: open file failure");
    }
    fprintf(fp,
            "#!/bin/tcsh\n"
            /* These are for seaborg.nersc.gov */
            //"setenv LD_LIBRARY_PATH ~/Tools/Lib:\n"
            //"module add gnuplot\n"
            "gnuplot plotenghist.gp -\n");
    fclose(fp);
    
    strcpy(comsav,command);
    sprintf(command,"chmod u+x gpw"); runcommand();
    sprintf(command,"xterm -fn 7x13 -sb -e ./gpw &\n");
    runcommand();
    strcpy(command,comsav);
}

















/********************************************************************/
/* Visualization */
/********************************************************************/

void KMCFrame::openwin()
{
    openwindow(win_width,win_height,dirname);
}

int KMCFrame::openwindow(int w,int h,const char *n)
{
    if(win!=NULL)
        if(win->alive) return -1;

    win=new YWindow(w,h,n,true,true,false);
    if(win==NULL) return -1;
    else
    {
        win->setinterval(100);
        win->Run();
        return 0;
    }
}

void KMCFrame::closewindow() {delete(win);}

void KMCFrame::winplot()
{
    if(win!=NULL)
        if(win->alive)
            if((curstep%plotfreq)==0)
            {
                plot();
            }
}

void KMCFrame::winplot(int step)
{
    if(win!=NULL)
        if(win->alive)
            if((step%plotfreq)==0)
            {
                plot();
            }
}        

void KMCFrame::wintogglepause()
{
    if(win!=NULL)
        if(win->IsAlive())
            win->TogglePause();
}

#include "namecolor.c"
void KMCFrame::alloccolors()
{
    int r,g,b;
    if(win!=NULL)
    {
        if(win->alive)
        {
            for(int i=0;i<MAXCOLORS;i++)
            {
                Str2RGB(colornames[i],&r,&g,&b);
                colors[i]=win->AllocShortRGBColor(r,g,b);
            }
            
            Str2RGB(atomcolor[0],&r,&g,&b);
            colors[MAXCOLORS+0]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(bondcolor,&r,&g,&b);
            colors[MAXCOLORS+1]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(highlightcolor,&r,&g,&b);
            colors[MAXCOLORS+2]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(fixatomcolor,&r,&g,&b);
            colors[MAXCOLORS+3]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(backgroundcolor,&r,&g,&b);
            win->bgcolor=win->AllocShortRGBColor(r,g,b);

            for(int i=0;i<MAXSPECIES;i++)
            {
                Str2RGB(atomcolor[i],&r,&g,&b);
                colors[MAXCOLORS+5+i]=win->AllocShortRGBColor(r,g,b);
            }
        }
        else
        {
            WARNING("No window to allocate color for!");
        }
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}
void KMCFrame::alloccolorsX()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            colors[0]=win->AllocNamedColor(atomcolor[0]);
            colors[1]=win->AllocNamedColor(bondcolor);
            colors[2]=win->AllocNamedColor(highlightcolor);
            colors[3]=win->AllocNamedColor(fixatomcolor);
            win->bgcolor=win->AllocNamedColor(backgroundcolor);
        }
        else
        {
            WARNING("No window to allocate color for!");
        }
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}

void KMCFrame::rotate()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            win->horizontalRot(rotateangles[0]*M_PI/180);
            win->verticalRot(rotateangles[1]*M_PI/180);
            win->spinRot(rotateangles[2]*M_PI/180);
            if(rotateangles[3]!=0) win->zoom(rotateangles[3]);
            if(rotateangles[4]!=0) win->project(rotateangles[4]);
        }
        else
        {
            WARNING("No window to rotate for!");
        }
    }
    else
    {
        WARNING("No window to rotate for!");
    }
}

void KMCFrame::saverot()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            //win->saveRot(); win->saveScale();
            win->saveView();
        }
        else
        {
            WARNING("No window to rotate for!");
        }
    }
    else
    {
        WARNING("No window to rotate for!");
    }
}

#include "colormap.h"
  
void KMCFrame::plot()
{
    int i,j,k, r,g,b, nw, cind;
    unsigned ce; bool high; char s[100];
    double Emin, Emax, L, alpha;
    Vector3 s1, s2, vr1, vr2, sri, srj, ri, rj;

#if 0 /* draw bonds */
    int jp, show;
    double r2, x1,y1,z1, x2,y2,z2, dx,dy,dz, dr;
#endif

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=max(_H[0][0],_H[1][1]);
    L=max(L,_H[2][2])*.5;

    switch(plot_color_axis){
     case(0): color_ind=_EPOT_IND; break;
     case(1): color_ind=_EPOT_IND; break;
     case(2): color_ind=_TOPOL; break;
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
        
        /* if plot_limits[0] equals to
           0: then this option is ignored
           1: then only plot atoms falling within the limits
              given by plot_limits[1~6]
           2: then atoms falling within the limits will be plotted
              with the same color as the fixed atoms
        */
        if(plot_limits[0]==1)
            if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6]))
            {
                continue;
            }
        if(plot_limits[0]==2)
            if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6]))
            {
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,colors[MAXCOLORS+3],s,2);
                continue;
            }
        
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
            else
                s[0]=0;
        }
        
        nw = (int)plot_color_windows[0];
        if(nw>0)
        {
            for(k=0;k<nw;k++)
            {
                Emin = plot_color_windows[k*3+1];
                Emax = plot_color_windows[k*3+2];
                cind = (int) plot_color_windows[k*3+3];

                if((color_ind[i]>=Emin)&&(color_ind[i]<=Emax))
                {
                    
                    if(plot_color_bar[0]==0)
                    {
                        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                                       atomradius[species[i]]/L,colors[cind],s,2);
                    }
                    else
                    { /* color atoms according to their properties given by color_ind
                         which is specified by coloraxis
                         there are two different color maps as specified by
                         plot_color_bar[0]: 1 or 2
                      */
                        alpha=(color_ind[i]-plot_color_bar[1])
                            /(plot_color_bar[2]-plot_color_bar[1]);
                        if(plot_color_bar[0]==1)
                            colormap1(alpha,&r,&g,&b);
                        else colormap2(alpha,&r,&g,&b);
                        
                        ce=win->AllocShortRGBColor(r,g,b);
                        if(plot_atom_info==5) /* print atom color (rgb) if clicked */
                            sprintf(s,"%d,%d,%d,%d,%x",i,r,g,b,ce);
                        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                                       atomradius[species[i]]/L,ce,s,2);
                    }                
                }
            }
            if(plot_color_windows[nw*3+1]) /* draw fixed atoms */
                if(fixed[i])   /* plot fixed atoms */
                {
                    //INFO_Printf("draw fixed atom [%d]\n",i);
                    win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                         atomradius[species[i]]/L,colors[MAXCOLORS+3],s,2);
                }
            
        }
        else
        {
            high=false;
            for(j=1;j<=plot_highlight_atoms[0];j++)
                if(i==plot_highlight_atoms[j])
                {
                    high=true;
                    break;
                }
            
            if(fixed[i])   /* plot fixed atoms */
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,colors[MAXCOLORS+3],s,2);
            else if (high) /* plot highlight atoms */
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,colors[MAXCOLORS+2],s,2);
            
            else if(plot_color_bar[0]!=0)
            { /* color atoms according to their properties given by color_ind
                 which is specified by coloraxis
                 there are two different color maps as specified by
                 plot_color_bar[0]: 1 or 2
              */
                alpha=(color_ind[i]-plot_color_bar[1])
                        /(plot_color_bar[2]-plot_color_bar[1]);
                if(plot_color_bar[0]==1)
                    colormap1(alpha,&r,&g,&b);
                else colormap2(alpha,&r,&g,&b);
                
                ce=win->AllocShortRGBColor(r,g,b);
                if(plot_atom_info==5) /* print atom color (rgb) if clicked */
                    sprintf(s,"%d,%d,%d,%d,%x",i,r,g,b,ce);
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,ce,s,2);
            }
            else
            {
                /* otherwise color atoms according to their species */
                ce=colors[MAXCOLORS+5+species[i]];
                if(plot_atom_info==5) /* print atom species and color code if clicked */
                    sprintf(s,"%d,%d,%x",i,species[i],ce);
                
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,ce,s,2);
            }
        }
    }
#if 0    
    /* Draw Bonds */
    if(bondlength>1e-3)
    {
        for(i=0;i<_NP;i++)
        {
            sri=_SR[i];
            if(plot_map_pbc==1) sri.subint();
            ri = _H*sri;
            
            show = 0;
            nw = (int)plot_color_windows[0];
            if(nw>0)
            {
                for(k=0;k<nw;k++)
                {
                    Emin = plot_color_windows[k*3+1];
                    Emax = plot_color_windows[k*3+2];
                    cind = (int) plot_color_windows[k*3+3];
                    
                    if((color_ind[i]>=Emin)&&
                       (color_ind[i]<=Emax))
                    {
                        show = 1;
                        break;
                    }
                }
            }
            else
            {
                show = 1;
            }
            if(!show) continue;
            
            for(jp=0;jp<nn[i];jp++)
            {
                j=nindex[i][jp];
                srj=_SR[j];
                if(plot_map_pbc==1) srj.subint();
                rj = _H*srj;

                show = 0;
                nw = (int)plot_color_windows[0];
                if(nw>0)
                {
                    for(k=0;k<nw;k++)
                    {
                        Emin = plot_color_windows[k*3+1];
                        Emax = plot_color_windows[k*3+2];
                        cind = (int) plot_color_windows[k*3+3];
                        
                        if((color_ind[j]>=Emin)&&
                           (color_ind[j]<=Emax))
                        {
                            show = 1;
                            break;
                        }
                    }
                }
                else
                {
                    show = 1;
                }
                if(!show) continue;
                
                if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                        continue;
                r2=(ri-rj).norm2();
                //INFO_Printf(" %d - %d : %f %f\n",i,j,r2,bondlength*bondlength);
                if(r2<bondlength*bondlength)
                {
                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 
                    x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
                    x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
                    dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
                    dx/=dr;dy/=dr;dz/=dr;
                    win->DrawLine(x1+dx*atomradius[species[i]]/L,
                                  y1+dy*atomradius[species[i]]/L,
                                  z1+dz*atomradius[species[i]]/L,
                                  x2-dx*atomradius[species[j]]/L,
                                  y2-dy*atomradius[species[j]]/L,
                                  z2-dz*atomradius[species[j]]/L,
                                  colors[MAXCOLORS+1],bondradius/L,1);
                }
            }
        }
    }
#endif
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

/* set variables from arg list */
void KMCFrame::set_dirname(char *s)
{
    strcpy(dirname,s);
}

void KMCFrame::set_randseed(char *s)
{
    sscanf(s,"%d",&_RANDSEED);
}


/********************************************************************/
/* End of class KMCFrame definition */
/********************************************************************/




/********************************************************************/
/* Main program starts below */
/********************************************************************/

#ifdef _TEST

class TESTKMC : public KMCFrame
{
public:
    virtual void int_potential()
    {
        /* Multi process function */
        /* a void potential */
        INFO("Empty potential");
    }
    virtual void initvars()
    {
        KMCFrame::initvars();
    }
            
};

class TESTKMC sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
