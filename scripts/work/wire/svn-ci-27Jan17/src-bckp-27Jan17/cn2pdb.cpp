/*
  cn2pdb.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Fri May 25 23:00:43 2001

  FUNCTION  :  Translating atomic configuration files CNFile (.cn)
               to Protein Data Bank (.pdb) files for Rasmol viewer
               and vpdb
               
  cn2pdb -Si -r1.2 *.cn.bz2
  cn2pdb -Mo -r1.1 *.cn.bz2
*/

#include <stdio.h>
#include "general.h"
#include "linalg3.h"
#include "filecls.h"

#define NPM 1000000
/*length in Angstrom, 1au = 0.5291772083 \AA */
/*
  static char *AtomName[]=
  {
  "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
  "P","S","Cl","Ar","K","Ca",
  };
*/
class Converter{
public:
    int _NP;
    char matname[100];
    int reflash; double ex,ey,ez,sx0,sy0,sz0;
    double ratio;
    class Vector3 _SR[NPM];//, _F[NPM], _VSR[NPM];
//positions of and forces on atoms
    int _LATOMTYPE[NPM];
    class Matrix33 _H, _VIRIAL;

    Converter():reflash(1),ex(1),ey(1),ez(1),sx0(0),sy0(0),sz0(0),ratio(1){};
    int readcn(char *name)
    {
        int i;
        char *buffer; char *pp,fname[100];
        class LFile *f;
        
        strcpy(fname,name);
        f=LFile::Open(fname,LFile::O_Read);
        if(f==NULL)
        {
            sprintf(fname,"%s.bz2",name);
            f=LFile::Open(fname,LFile::O_Read);
            if(f==NULL)
            {
                sprintf(fname,"%s.gz",name);
                f=LFile::Open(fname,LFile::O_Read);
                if(f==NULL) return -1;
            }
        }
        
        LFile::LoadToString(fname,buffer,0);
        
        pp=buffer;
        sscanf(pp, "%d", &_NP); //INFO("NP="<<_NP);
        
        pp=strchr(pp, '\n');
        pp++;
        for(i=0;i<_NP;i++)
        {
            char *q;
            q=pp;
            pp=strchr(pp,'\n');
            if(pp) *(char *)pp++=0;
            sscanf(q, "%lf %lf %lf",
                   &(_SR[i].x),&(_SR[i].y),&(_SR[i].z));
        }

        sscanf(pp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &(_H[0][0]),&(_H[0][1]),&(_H[0][2]),
               &(_H[1][0]),&(_H[1][1]),&(_H[1][2]),
               &(_H[2][0]),&(_H[2][1]),&(_H[2][2]));

        free(buffer);
        LFile::Close(f);
        return 0;
    }
    void reflashS()
        { for(int i=0;i<_NP;i++) _SR[i].subint(); }
    int writepdb(char *fname)
    {
        double a, b, c, alpha, beta, gamma;
        double b1,b2,c1,c2,c3;
        Vector3 va, vb, vc;   Matrix33 h0;
        FILE *f;

        f=fopen(fname,"w");
        
        va.set(_H[0][0],_H[1][0],_H[2][0]); va*=ratio;
        vb.set(_H[0][1],_H[1][1],_H[2][1]); vb*=ratio;
        vc.set(_H[0][2],_H[1][2],_H[2][2]); vc*=ratio;
        a=va.norm();
        b=vb.norm();
        c=vc.norm();
        alpha=acos(vb*vc/b/c);
        beta=acos(vc*va/c/a);
        gamma=acos(va*vb/a/b);
        fprintf(f,"HEADER    cn2pdb NP=%d\n",_NP);
        fprintf(f,"HEADER  500 MUST_PBC\n");
        fprintf(f,"CRYST1%9.3lf%9.3lf%9.3lf"
                "%7.2lf%7.2lf%7.2lf P 1           1\n",
                a,b,c,alpha*180/M_PI,beta*180/M_PI,gamma*180/M_PI);
//                a,b,c,90.0,90.0,90.0);
        b1=b*cos(gamma);
        b2=b*sin(gamma);
        c1=c*cos(beta);
        c2=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
        c3=sqrt(c*c-c1*c1-c2*c2);
        h0.set(a,b1,c1,
               0,b2,c2,
               0,0,c3);
//                fprintf(stdout,"%f %f %f\n%f %f %f\n%f %f %f\n",
//                        a,b1,c1,0.0,b2,c2,0.0,0.0,c3);
        for(int i=0;i<_NP;i++)
        {
            Vector3 r,s;
            s=_SR[i];  

            s.subint();
            s.x+=0.5; s.y+=0.5; s.z+=0.5;
//            if(s.x>0.99999)s.x=0.0001;
//            if(s.y>0.99999)s.y=0.0001;
//            if(s.z>0.99999)s.z=0.0001;
//            if(s.x<0.0001)s.x=0.0001;
//            if(s.y<0.0001)s.y=0.0001;
//            if(s.z<0.0001)s.z=0.0001;
            
            r=h0*s;
            
            fprintf(f,"ATOM  %5d %2s"
                    "           1    %8.3f%8.3f%8.3f\n",
                    i+1,
                    ((_LATOMTYPE[i]==0)?matname:
                     ((_LATOMTYPE[i]==16)?"O":"Si")),
                    r.x,r.y,r.z);
        }
        fclose(f);
        LFile::BZipFile(fname,true);
        return 0;
    }        
} con;

int cn2pdb(int argv, char *argc[])
{
    char fname[100], *ext;
    int i,n;
    if(argv==1)
    {
        fprintf(stderr,"Usage: %s [-<NAME>] [-f] [-e[x,y,z]<n>]"
                "[-s(x,y,z)<n>] [-r<n>] filename.cn\n",argc[0]);
        return -1;
    }
    sprintf(con.matname,"%s","Mo"); //Mo as default material name
    n=argv;
    for(i=1;i<n;i++)
    {
        if(argc[i][0]=='-')
        {
            if((argc[i][1]=='f')&&(argc[i][2]==0))
            {
                con.reflash=0; continue;
            }
            if((argc[i][1]=='e')&&(argc[i][2]!=0))
            {
                if((argc[i][2]!='x')&&(argc[i][2]!='y')&&(argc[i][2]!='z'))
                {
                    sscanf(argc[i]+2,"%lf",&con.ex);
                    con.ey=con.ez=con.ex;continue;
                }
                else
                {
                    switch(argc[i][2]) {
                    case('x'): sscanf(argc[i]+3,"%lf",&con.ex);continue; 
                    case('y'): sscanf(argc[i]+3,"%lf",&con.ey);continue; 
                    case('z'): sscanf(argc[i]+3,"%lf",&con.ez);continue; 
                    }
                }
            }
            if((argc[i][1]=='s')&&(argc[i][2]!=0))
            {
                if((argc[i][2]!='x')&&(argc[i][2]!='y')&&(argc[i][2]!='z'))
                {
                    continue;
                }
                else
                {
                    switch(argc[i][2]) {
                    case('x'): sscanf(argc[i]+3,"%lf",&con.sx0);continue; 
                    case('y'): sscanf(argc[i]+3,"%lf",&con.sy0);continue; 
                    case('z'): sscanf(argc[i]+3,"%lf",&con.sz0);continue;
                    }
                }
            }
            if((argc[i][1]=='r')&&(argc[i][2]!=0))
            {
                sscanf(argc[i]+2,"%lf",&con.ratio);
                continue;
            }
            sscanf(argc[i]+1,"%s",con.matname);
//            fprintf(stderr,"%s",con.matname);
            continue;
        }

        
        sprintf(fname,"%s",argc[i]);

        ext=strstr(fname, ".cn");
        if(ext==NULL)
        {
            fprintf(stderr,HIR"invalid file"NOR"(not .cn) %s\n",fname);
            continue;
        }
            
        con.readcn(fname);

        if(con.reflash) con.reflashS();

        fprintf(stderr,"%s -> ",fname);
        ext[1]='p';
        ext[2]='d';
        ext[3]='b';
        ext[4]=0;
        fprintf(stderr,"%s\n",fname);
        con.writepdb(fname);
    }
    return 0;
}

int main(int argv, char *argc[])
{
//#define _CN2PDB
//#define _MCN2CN
//#define _VCN2CN
    
#ifdef _CN2PDB
    cn2pdb(argv,argc);
#elif defined(_MCN2CN)
    maurice2cn(argv,argc);
#elif defined(_VCN2CN)
    vasily2cn(argv,argc);
#elif defined(_LCN2PDB)
    lcn2pdb(argv,argc);
#elif defined(_MCN2PDB)
    mcn2pdb(argv,argc);
#endif
    return 0;
}

