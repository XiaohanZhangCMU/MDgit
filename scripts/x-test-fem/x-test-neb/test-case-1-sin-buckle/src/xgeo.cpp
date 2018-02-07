/*
  xgeo.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Jan 11 11:58:45 2007

  FUNCTION  :
*/

#include "xgeo.h"

void XGeo::initparser()
{
    bindvar("geo_box",geo_box,DOUBLE);
    bindvar("geo_line",geo_line,DOUBLE);
    bindvar("geo_plane",geo_plane,DOUBLE);
    bindvar("geo_xpoint",geo_xpoint,INT);
    bindvar("geo_xline",geo_xline,INT);
    bindvar("geo_vector",geo_vector,DOUBLE);

    bindvar("win_width",&win_width,INT);
    bindvar("win_height",&win_height,INT);
    bindvar("pointradius",&pointradius,DOUBLE);
    bindvar("geo_box_color",geo_box_color,STRING);
    bindvar("geo_line_color",geo_line_color,STRING);
    bindvar("geo_plane_color",geo_plane_color,STRING);
    bindvar("geo_xpoint_color",geo_xpoint_color,STRING);
    bindvar("geo_xline_color",geo_xline_color,STRING);
    bindvar("geo_vector_color",geo_vector_color,STRING);
    bindvar("backgroundcolor",backgroundcolor,STRING);
}

void XGeo::initvars()
{
    initparser();
}

int XGeo::exec(const char *name)
{
    if(Organizer::exec(name)==0) return 0;
    bindcommand(name,"openwin",openwin());
    bindcommand(name,"plot",plot());
    bindcommand(name,"alloccolors",alloccolors());
    return -1;
}

#ifdef _USETCL
/********************************************************************/
/* Initialize Tcl parser */
/********************************************************************/
int XGeo::Tcl_AppInit(Tcl_Interp *interp)
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

    /*
      Random_Init(interp);
      Blob_Init(interp);
    */
    return TCL_OK;
}

int Tcl_Main_Wrapper(int argc, char *argv[])
{
#define MAXNARGS  20
#define MAXARGLEN 500    
//    char newargv[MAXNARGS][MAXARGLEN];

    if(argc>MAXNARGS)
        ERROR("maximum number of arguments exceeded!");

    Tcl_Main(argc, argv, Tcl_AppInit);
//    Tcl_Main(argc, (char **)newargv, Tcl_AppInit);

    return TCL_OK;
}
#endif


#include "namecolor.c"
void XGeo::alloccolors()
{
    int r,g,b;
    if ((win!=NULL)&&(win->alive))
    {
        Str2RGB(geo_box_color,&r,&g,&b);
        colors[0]=win->AllocShortRGBColor(r,g,b);
        Str2RGB(geo_line_color,&r,&g,&b);
        colors[1]=win->AllocShortRGBColor(r,g,b);
        Str2RGB(geo_plane_color,&r,&g,&b);
        colors[2]=win->AllocShortRGBColor(r,g,b);
        Str2RGB(geo_xpoint_color,&r,&g,&b);
        colors[3]=win->AllocShortRGBColor(r,g,b);
        Str2RGB(geo_xline_color,&r,&g,&b);
        colors[4]=win->AllocShortRGBColor(r,g,b);
        Str2RGB(geo_vector_color,&r,&g,&b);
        colors[5]=win->AllocShortRGBColor(r,g,b);
        Str2RGB(backgroundcolor,&r,&g,&b);
        win->bgcolor=win->AllocShortRGBColor(r,g,b);
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}

void XGeo::openwin()
{
    openwindow(win_width,win_height,dirname);
    alloccolors();
}

int XGeo::openwindow(int w,int h,const char *n)
{
    if(win!=NULL)
        if(win->alive) return -1;
    win=new YWindow(w,h,n);
    if(win==NULL) return -1;
    else
    {
        win->setinterval(100);
        win->Run();
        return 0;
    }
}

void XGeo::closewindow() {delete(win);}

void XGeo::plot()
{
    if((win!=NULL)&&(win->alive))
    {
        win->Lock();
        win->Clear();
        drawBox();
        drawLines();
        drawVectors();
        drawPlanes();
        drawXpoints();
        drawXlines();
    }
    win->Unlock();
    win->Refresh();
}

void XGeo::drawBox()
{
    Vector3 a,b,c,e1,e2,e3,r1,r2;
    
    a.set(geo_box[0][0],geo_box[0][1],geo_box[0][2]);
    b.set(geo_box[1][0],geo_box[1][1],geo_box[1][2]);
    c.set(geo_box[2][0],geo_box[2][1],geo_box[2][2]);
    a*=geo_box[0][3];
    b*=geo_box[1][3];
    c*=geo_box[2][3];
    L=a.norm();
    if(L<b.norm()) L=b.norm();
    if(L<c.norm()) L=c.norm();

    L/=2;
//    a/=L; b/=L; c/=L;
    _H.setcol(a,b,c);

    e1=a/a.norm();
    e2=b/b.norm(); e2=schmidt(e2,e1);
    e3=cross(e1,e2);
    _RV.setcol(e1,e2,e3); _RV=_RV.inv(); _RV/=L;
    
#define _Geo_DoDraw(_1,_2,_3,_4,_5,_6) \
    r1=(a*(_1 1)+b*(_2 1)+c*(_3 1))*0.5; \
    r2=(a*(_4 1)+b*(_5 1)+c*(_6 1))*0.5; \
    r1=_RV*r1; r2=_RV*r2; \
    win->DrawLine(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z,colors[0]);

    _Geo_DoDraw(+,-,+, +,+,+);
    _Geo_DoDraw(-,-,+, +,-,+);
    _Geo_DoDraw(-,+,+, -,-,+);
    _Geo_DoDraw(+,+,+, -,+,+);
    _Geo_DoDraw(+,-,-, +,+,-);
    _Geo_DoDraw(-,-,-, +,-,-);
    _Geo_DoDraw(-,+,-, -,-,-);
    _Geo_DoDraw(+,+,-, -,+,-);
    _Geo_DoDraw(+,+,+, +,+,-);
    _Geo_DoDraw(+,-,+, +,-,-);
    _Geo_DoDraw(-,-,+, -,-,-);
    _Geo_DoDraw(-,+,+, -,+,-);
}

void XGeo::drawLines()
{
    int i, nlines;
    Vector3 dl, sr0, r1, r2;
    
    nlines=(int)geo_line[0];
    if(nlines<=0) return;

    for(i=0;i<nlines;i++)
    {
        dl.set(geo_line[i*7+2],geo_line[i*7+3],geo_line[i*7+4]);
        sr0.set(geo_line[i*7+5],geo_line[i*7+6],geo_line[i*7+7]);
        getlineendpoints(dl,sr0,r1,r2);
        r1=_RV*r1; r2=_RV*r2;
        win->DrawLine(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z,colors[1],0.015);
    }
}

void XGeo::drawVectors()
{
    int i, nvectors;
    Vector3 dl, sr0, r1, r2;
    double scale;
    
    nvectors=(int)geo_vector[0];
    if(nvectors<=0) return;

    for(i=0;i<nvectors;i++)
    {
        dl.set(geo_vector[i*8+2],geo_vector[i*8+3],geo_vector[i*8+4]);
        sr0.set(geo_vector[i*8+5],geo_vector[i*8+6],geo_vector[i*8+7]);
        scale=geo_vector[i*8+8];
        r1=_H*sr0;
        r2=r1+dl*scale/L/dl.norm();
        r1=_RV*r1; r2=_RV*r2;
        win->DrawLine(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z,colors[5],0.02);
    }
}

bool XGeo::getlineendpoints(Vector3 &dl, Vector3 &sr0,
                            Vector3 &r1, Vector3 &r2)
{
    Vector3 sl, sr[6], dr; Matrix33 hinv;
    double x, y, z, a;
    int i, n; bool succ;

//    INFO("dl ="<<dl);
//    INFO("sr0="<<sr0);
    hinv=_H.inv();
    sl=hinv*dl;
    x=y=z=0;

#define testendpoint(_1,_2,_3)                             \
    if(sl._1!=0)                                           \
    {                                                      \
        a=(0.5-sr0._1)/sl._1;                              \
        _1=0.5;                                            \
        _2=a*sl._2+sr0._2;                                 \
        _3=a*sl._3+sr0._3;                                 \
        if((fabs(_2) <= 0.5001)&&(fabs(_3) <= 0.5001))     \
        {                                                  \
            sr[n].set(x,y,z); n++;                         \
        }                                                  \
        a=(-0.5-sr0._1)/sl._1;                             \
        _1=-0.5;                                           \
        _2=a*sl._2+sr0._2;                                 \
        _3=a*sl._3+sr0._3;                                 \
        if((fabs(_2) <= 0.5001)&&(fabs(_3) <= 0.5001))     \
        {                                                  \
            sr[n].set(x,y,z); n++;                         \
        }                                                  \
    }                                                      

    n=0;
    /* x=+-0.5 */
    testendpoint(x,y,z);
    
    /* y=+-0.5 */
    testendpoint(y,x,z);

    /* z=+-0.5 */
    testendpoint(z,x,y);

    if(n<2)
    {
        WARNING("getlineendpoints n="<<n<<" !=2");
        r1.clear(); r2.clear();
        return false;
    }
    
    succ=false;
    for(i=1;i<n;i++)
    {
        dr=sr[i]-sr[0];
        if(dr.norm2()>1e-6)
        {
            succ=true;
            break;
        }
    }
    sr[1]=sr[i];
    
    if(!succ)
    {
        WARNING("getlineendpoints failed");
        r1.clear(); r2.clear();
        return false;
    }
    else
    {
        r1=_H*sr[0];
        r2=_H*sr[1];
//        INFO("sr[0]="<<sr[0]);
//        INFO("sr[1]="<<sr[1]);
//        INFO("r1="<<r1);
//        INFO("r2="<<r2);
        return true;
    }
}

bool XGeo::getxline(Vector3 &dn1, Vector3 &sr1,
                    Vector3 &dn2, Vector3 &sr2,
                    Vector3 &dl,  Vector3 &sr)
{
    Vector3 r1, r2, r;
    Vector3 e1, e2;
    Vector3 x, b;
    Matrix33 A, Ainv, hinv;
    
    hinv=_H.inv();
    r1=_H*sr1;
    r2=_H*sr2;            
    dl=cross(dn1, dn2);
//    INFO("dl="<<dl);               
    if(dl.norm2()>1e-6)
    {                             
        e2=cross(dn2,dl); e1=cross(dn1,dl);               
        A.setcol(e1,e2,dl); Ainv=A.inv();                
        b=r2-r1;                                         
        x=Ainv*b;                                        
        r=r1+e1*x[0]; sr=hinv*r;
//        INFO("getxline: sr="<<sr);
        return true;
    }
    else
    {
        sr.clear();
        return false;
    }
}

void XGeo::drawPlanes()
{
    int i, nplanes;
    Vector3 dn, sr1, ra, rb;
    Vector3 dn2, sr2;
    Vector3 dl, srl;
    unsigned c; double alpha;
    
    nplanes=(int)geo_plane[0];
    if(nplanes<=0) return;

    for(i=0;i<nplanes;i++)
    {
        dn.set(geo_plane[i*7+2],geo_plane[i*7+3],geo_plane[i*7+4]);
        sr1.set(geo_plane[i*7+5],geo_plane[i*7+6],geo_plane[i*7+7]);

        alpha=1.0*i/nplanes/2;
        c=(unsigned)(colors[2]*(1-alpha)+colors[5]*alpha);
        
#define testendline(_1,_2)                                \
        dn2.clear(); dn2._1=1;                            \
        sr2.clear(); sr2._1=_2 0.5;                       \
        if(getxline(dn,sr1,dn2,sr2,dl,srl)){              \
        if(getlineendpoints(dl,srl,ra,rb)){               \
        ra=_RV*ra; rb=_RV*rb;                             \
        win->DrawLine(ra.x,ra.y,ra.z,rb.x,rb.y,rb.z,c);}}

        /* x=+-0.5 */
        testendline(x,+);
        testendline(x,-);
        /* y=+-0.5 */
        testendline(y,+);
        testendline(y,-);
        /* z=+-0.5 */
        testendline(z,+);
        testendline(z,-);
    }
}

bool XGeo::getxpoint(Vector3 &dl1, Vector3 &sr1,
                     Vector3 &dn2, Vector3 &sr2, Vector3 &r)
{
    Vector3 r1, r2, dr;
    double a, t;
    
    r1=_H*sr1;
    r2=_H*sr2;

//    INFO("getxpoint: dl1="<<dl1<<"  r1="<<r1);
//    INFO("getxpoint: dn2="<<dn2<<"  r2="<<r2);

    dr=r2-r1;
    t=dl1*dn2;
    if(fabs(t)<1e-6)
    {
        r.clear();
        return false;
    }
    else
    {
        a=(dr*dn2)/t;
        r=r1+dl1*a;
        return true;
    }
}

void XGeo::drawXpoints()
{
    int i, nxpoints, l1, p2;
    Vector3 dl1, sr1;
    Vector3 dn2, sr2;
    Vector3 r;
    
    nxpoints=(int)geo_xpoint[0];
//    INFO("nxpoints="<<nxpoints);
    if(nxpoints<=0) return;

    for(i=0;i<nxpoints;i++)
    {
        l1=geo_xpoint[i*3+2]; l1--;
        p2=geo_xpoint[i*3+3]; p2--;

//        INFO("drawXlines: l1="<<l1<<"  p2="<<p2);
        
        dl1.set(geo_line[l1*7+2],geo_line[l1*7+3],geo_line[l1*7+4]);
        sr1.set(geo_line[l1*7+5],geo_line[l1*7+6],geo_line[l1*7+7]);
        dn2.set(geo_plane[p2*7+2],geo_plane[p2*7+3],geo_plane[p2*7+4]);
        sr2.set(geo_plane[p2*7+5],geo_plane[p2*7+6],geo_plane[p2*7+7]);

        if(getxpoint(dl1,sr1,dn2,sr2,r))
        {
//            INFO("xpoint: r="<<r);
            r=_RV*r;
            win->DrawPoint(r.x,r.y,r.z,pointradius/L,colors[3]);
        }
    }
}

void XGeo::drawXlines()
{
    int i, nxlines, p1, p2;
    Vector3 dn1, sr1, ra, rb;
    Vector3 dn2, sr2;
    Vector3 dl, sr;
    
    nxlines=(int)geo_xline[0];
//    INFO("nxlines="<<nxlines);
    if(nxlines<=0) return;

    for(i=0;i<nxlines;i++)
    {
        p1=geo_xline[i*3+2]; p1--;
        p2=geo_xline[i*3+3]; p2--;

//        INFO("drawXlines: p1="<<p1<<"  p2="<<p2);
        
        dn1.set(geo_plane[p1*7+2],geo_plane[p1*7+3],geo_plane[p1*7+4]);
        sr1.set(geo_plane[p1*7+5],geo_plane[p1*7+6],geo_plane[p1*7+7]);
        dn2.set(geo_plane[p2*7+2],geo_plane[p2*7+3],geo_plane[p2*7+4]);
        sr2.set(geo_plane[p2*7+5],geo_plane[p2*7+6],geo_plane[p2*7+7]);

        if(getxline(dn1,sr1,dn2,sr2,dl,sr))
        {
            if(getlineendpoints(dl,sr,ra,rb))
            {
                ra=_RV*ra; rb=_RV*rb;                             
                win->DrawLine(ra.x,ra.y,ra.z,rb.x,rb.y,rb.z,colors[4]);
            }
        }
    }
}

class XGeo sim;

/* The main program is defined here */
#include "main.cpp"

