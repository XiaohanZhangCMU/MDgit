/*
  xgeo.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Dec 30 21:20:52 2006

  FUNCTION  :
*/

#ifndef _XGEO_H
#define _XGEO_H

#include "mdparallel.h"

#define MAXGEOOBJ 1000
class XGeo : public Organizer
{
public:
    double geo_box[3][4];
    double geo_line[7*MAXGEOOBJ+5];
    double geo_plane[7*MAXGEOOBJ+5];
    int geo_xpoint[3*MAXGEOOBJ+3];
    int geo_xline[3*MAXGEOOBJ+3];
    double geo_vector[8*MAXGEOOBJ+3];

    class YWindow *win;
    int win_width, win_height;

    double pointradius;
    char geo_box_color[30], geo_line_color[30], geo_plane_color[30];
    char geo_vector_color[30];
    char geo_xpoint_color[30], geo_xline_color[30], backgroundcolor[30];
    unsigned colors[10];

private:
    double L; /* size of box */
    Matrix33 _H, _RV;

    
public:    
    XGeo():win_width(500),win_height(300),pointradius(0.1),L(1){};
    virtual ~XGeo(){};
    
    virtual void initvars();
    virtual void initparser();
    virtual int exec(const char *name);            
#ifdef _TCL    
    virtual int Tcl_AppInit(Tcl_Interp *interp);
#endif

    int openwindow(int w,int h,const char *n);
    void closewindow();
    virtual void openwin();

    virtual void plot();
    void wintogglepause();
    void alloccolors();
    void rotate();
    void saverot();

private:
    void drawBox();
    void drawLines();
    void drawPlanes();
    void drawXpoints();
    void drawXlines();
    void drawVectors();

    bool getlineendpoints(Vector3 &, Vector3 &, Vector3 &, Vector3 &);
    bool getxline(Vector3 &, Vector3 &, Vector3 &,
                  Vector3 &, Vector3 &, Vector3 &);
    bool getxpoint(Vector3 &, Vector3 &, Vector3 &,
                   Vector3 &, Vector3 &);
};

#endif // _XGEO_H

