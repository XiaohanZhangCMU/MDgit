/*
  vaspbox.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Jan  1 21:29:35 2007

  FUNCTION  : relaxation wrapper for VASP
*/

#include "mdparallel.h"


class VaspBox : public MDPARALLELFrame
{
public:
    char outcarfile[500];
    char zxcgrfile[500];
    double vaspstress[6], dstress[6], dstrain[6];
    double Scompliance[6][6];
    Matrix33 gh;
    double vaspboxstepsize;

    VaspBox():vaspboxstepsize(1){};
    
    void potential()
    {
        /* Multi process function */
        /* a void potential */
        INFO("Empty potential");
    }
    virtual void initvars()
    {
        MDPARALLELFrame::initvars();
        _RLIST=3.8;
    }
    virtual void initparser();
    virtual int exec(char *name);
    
    void readoutcarstress();
    int searchstringfromfile(char *inLine, const char *match,
                              int len, FILE *fp);
    void calgh();
    void movebox();
    void scalebox();
    void writeZXCGR();
};

void VaspBox::initparser()
{
    MDPARALLELFrame::initparser();
    
    bindvar("outcarfile",outcarfile,STRING);
    bindvar("zxcgrfile",zxcgrfile,STRING);
    bindvar("Scompliance",Scompliance,DOUBLE);
    bindvar("vaspboxstepsize",&vaspboxstepsize,DOUBLE);
}

int VaspBox::exec(char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;
    bindcommand(name,"readoutcarstress",readoutcarstress());
    bindcommand(name,"movebox",movebox());
    bindcommand(name,"scalebox",scalebox());
    
    return -1;
}

int VaspBox::searchstringfromfile(char *inLine, const char *match,
                                   int len, FILE *fp)
{
    while(fgets(inLine, len, fp)!=NULL)
    {
        if(strncmp(inLine,match,strlen(match))==0)
        {
            INFO(inLine);
            return 0;
            break;
        }
    }
    return -1;
}

void VaspBox::readoutcarstress()
{
    FILE *fp;
    char extname[500], inLine[500], *p;
    int i, succ;
    
    INFO("readoutcarstress: "<<outcarfile);
    LFile::SubHomeDir(outcarfile,extname);
    INFO("readoutcarstress: "<<extname);

    fp=fopen(extname,"r");
    if(fp==NULL)
    {
        FATAL("readoutcarstress: open file failure");
    }

    inLine[0]=0; succ=0;
    while(1)
    {
        if(searchstringfromfile(inLine,
           "  energy  without entropy=", 500, fp)!=-1)
        {
            p=strrchr(inLine, '=');
            sscanf(p+1,"%lf",&_EPOT);
        }
        if(searchstringfromfile(inLine, "  FORCE on cell", 500, fp)==-1)
            break ;
        searchstringfromfile(inLine, "  Direction", 500, fp);
        searchstringfromfile(inLine, "  in kB", 500, fp);
        sscanf(inLine+strlen("  in kB"),"%lf %lf %lf %lf %lf %lf",
               vaspstress+0,vaspstress+1,vaspstress+2,
               vaspstress+3,vaspstress+4,vaspstress+5);
        for(i=0;i<6;i++)
            vaspstress[i]*=100;
//        for(i=0;i<3;i++)
        for(i=0;i<6;i++)
            vaspstress[i]*=-1;
        
        INFO_Printf("vaspstress (in MPa) ="
                    " (%3.3f %3.3f %3.3f %3.3f %3.3f %3.3f)\n",
                    vaspstress[0],vaspstress[1],vaspstress[2],
                    vaspstress[3],vaspstress[4],vaspstress[5]);
        succ=1;
    }
    if(!succ)
        FATAL("readoutcarstress: didn't find stress data");

    INFO("Energy = "<<_EPOT);
    
    fclose(fp);
}

void VaspBox::calgh()
{
    int i, j;
    double sum;
    
    dstress[0]=vaspstress[0]-_EXTSTRESS[0][0]; /* xx */
    dstress[1]=vaspstress[1]-_EXTSTRESS[1][1]; /* yy */
    dstress[2]=vaspstress[2]-_EXTSTRESS[2][2]; /* zz */
    dstress[3]=vaspstress[3]-_EXTSTRESS[0][1]; /* xy */
    dstress[4]=vaspstress[4]-_EXTSTRESS[1][2]; /* yz */
    dstress[5]=vaspstress[5]-_EXTSTRESS[2][0]; /* zx */

    for(i=0;i<6;i++)
    {
        sum = 0;
        for(j=0;j<6;j++)
            sum += Scompliance[i][j]*dstress[j];
        dstrain[i]=-sum;
    }
    gh.clear();

    gh[0][0] = _H[0][0] * dstrain[0] * vaspboxstepsize;
    gh[1][1] = _H[1][1] * dstrain[1] * vaspboxstepsize;
    gh[2][2] = _H[2][2] * dstrain[2] * vaspboxstepsize;
    gh[0][1] = _H[1][1] * dstrain[3] * vaspboxstepsize;
    gh[2][0] = _H[0][0] * dstrain[5] * vaspboxstepsize;
    gh[2][1] = _H[0][0] * dstrain[4] * vaspboxstepsize;
}

void VaspBox::writeZXCGR()
{
    FILE *fp;
    char extname[500], inLine[500]; //, *p;
    int n; //int i, j, succ;
    
    INFO("writeZXCGR: "<<zxcgrfile);
    LFile::SubHomeDir(zxcgrfile,extname);
    INFO("writeZXCGR: "<<extname);

    fp=fopen(extname,"r");
    conj_step = 0;
    if(fp==NULL)
    {
        fp=fopen(extname,"w");
    }
    else
    {
        while(searchstringfromfile(inLine,
                                   "step =", 500, fp)!=-1)
        {
            sscanf(inLine+strlen("step ="),"%d",&conj_step);
        }
        fclose(fp);
        INFO("last step = "<<conj_step);
        fp=fopen(extname,"a");
    }
    conj_step ++;
    
    fprintf(fp,"step = %d\n",conj_step);

    n = 9;

    fprintf(fp,"%d %e\n",n,_EPOT);
    
    fprintf(fp,"%e %e %e\n",_H[0][0],gh[0][0],vaspstress[0]);
    fprintf(fp,"%e %e %e\n",_H[1][0],gh[1][0],vaspstress[3]);
    fprintf(fp,"%e %e %e\n",_H[2][0],gh[2][0],vaspstress[5]);
    fprintf(fp,"%e %e %e\n",_H[0][1],gh[0][1],vaspstress[3]);
    fprintf(fp,"%e %e %e\n",_H[1][1],gh[1][1],vaspstress[1]);
    fprintf(fp,"%e %e %e\n",_H[2][1],gh[2][1],vaspstress[4]);
    fprintf(fp,"%e %e %e\n",_H[0][2],gh[0][2],vaspstress[5]);
    fprintf(fp,"%e %e %e\n",_H[1][2],gh[1][2],vaspstress[4]);
    fprintf(fp,"%e %e %e\n",_H[2][2],gh[2][2],vaspstress[2]);
    
//    for(j=0;j<3;j++)
//        for(i=0;i<3;i++)
//            fprintf(fp,"%e %e\n",_H[i][j],gh[i][j]);
    fclose(fp);
}

void VaspBox::movebox()
{
    int i, j;
    
    calgh();
    writeZXCGR();
    
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            _H[i][j]+=gh[i][j];
}

void VaspBox::scalebox()
{
    int i, j;
    double pres, invB;

    /* compute pressure */
    pres = ( (vaspstress[0]-_EXTSTRESS[0][0])
             +(vaspstress[1]-_EXTSTRESS[1][1])
             +(vaspstress[2]-_EXTSTRESS[2][2]) ) / 3;
    /* inverse Bulk modulus */
    invB = Scompliance[0][0] + Scompliance[0][1] + Scompliance[0][2];

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
        {
            gh[i][j] = -_H[i][j] * pres * invB * vaspboxstepsize;
        }
    
    writeZXCGR();
    
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            _H[i][j]+=gh[i][j];
}
/* main program begins */
class VaspBox sim;

/* The main program is defined here */
#include "main.cpp"

