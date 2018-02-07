/*
  display_win.h (Windows version of display.h)
  by David Park and Wei Cai  caiwei@stanford.edu
  Last Modified : Sat Aug 13 16:00:01 2005

  FUNCTION  :  X-window display package
*/

#ifndef _DISPLAY_WIN_H
#define _DISPLAY_WIN_H

#define nUSE_X  

#define _DISPLAY_VERSION 1.05

#ifdef NO_GENERAL
#define bool int
#define true 1
#define false 0
#endif

#include <stdio.h>
#include <stdlib.h>


#ifdef USE_X
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#ifndef NO_XPM
#include <X11/xpm.h>
#endif

#else

#define GLUT_STATIC
#include "glut.h"

/*------------------------------------------------------------*
 *              ascii codes for keys                          *
 *------------------------------------------------------------*/

#define KEY_ESCAPE     27
#define KEY_RETURN     13
#define KEY_BACKSPACE   8

#define KEY_ARROW_LEFT   100
#define KEY_ARROW_UP     101
#define KEY_ARROW_RIGHT  102
#define KEY_ARROW_DOWN   103

typedef enum {
  WINDOW_ACTION_NONE,
  WINDOW_ACTION_SCALE,
  WINDOW_ACTION_SCALE_X,
  WINDOW_ACTION_SCALE_Y,
  WINDOW_ACTION_SCALE_Z,
  WINDOW_ACTION_ROTATE_XY,
  WINDOW_ACTION_ROTATE_X,
  WINDOW_ACTION_ROTATE_Y,
  WINDOW_ACTION_ROTATE_Z,
  WINDOW_ACTION_TRANSLATE_XY,
  WINDOW_ACTION_TRANSLATE_X,
  WINDOW_ACTION_TRANSLATE_Y,
  WINDOW_ACTION_TRANSLATE_Z,
  WINDOW_ACTION_LOOKAT,
  WINDOW_ACTION_PICK,
  WINDOW_ACTION_KEY_PRESS,
  WINDOW_ACTION_MAX
  } WindowAction;

typedef enum {
  BUTTON_UNKNOWN,
  BUTTON_LEFT,
  BUTTON_MIDDLE,
  BUTTON_RIGHT,
  BUTTON_SHIFT_LEFT,
  BUTTON_SHIFT_MIDDLE,
  BUTTON_SHIFT_RIGHT,
  BUTTON_CTRL_LEFT,
  BUTTON_CTRL_MIDDLE,
  BUTTON_CTRL_RIGHT,
  BUTTON_MAX
  } EventButtonId;

typedef enum {
  EVENT_BUTTON_STATE_NONE,
  EVENT_BUTTON_STATE_SHIFT,
  EVENT_BUTTON_STATE_CONTROL,
  EVENT_BUTTON_STATE_CONTROL_SHIFT
  } EventButtonState;

typedef enum {
  EVENT_NONE,
  EVENT_BUTTON_PRESS,
  EVENT_BUTTON_RELEASE,
  EVENT_MOUSE_MOTION,
  EVENT_WINDOW_EXPOSE,
  EVENT_KEY_PRESS,
  EVENT_KEY_RELEASE
  } EventType;

typedef  enum {
  SPHERE_REP_DISC,
  SPHERE_REP_DOT,
  SPHERE_REP_POINT,
  SPHERE_REP_LINE,
  SPHERE_REP_SOLID,
  } SphereRep;

typedef struct GlutXmotion {
  int x, y;
  } GlutXmotion;

typedef struct XEvent {
  int mx, my;
  int mouse_pos_set;
  EventType type;
  int button_id;
  EventButtonState button_state;
  WindowAction action;
  void (*action_func)();
  GlutXmotion xmotion;
  } XEvent;

typedef char KeySym; 

void display (void);

void event_Butp (int button, int state, int x, int y);

void event_Butr (int mx, int my, EventButtonId button_id);

void event_IdleProc ();

void event_KeyboardProc (unsigned char ch, int x, int y);

void event_Motion (int x, int y);

void event_SpecialProc (int key, int x, int y);

void event_WinObjGet (int win, int *p);

void event_WinReshape (int width, int height);

void opengl_SurfGenSphere (SphereRep dtype, float radius, int slices, int stacks);

void viewport_Set (int win, int w, int h);

#endif



#ifndef NO_THREAD
#include <pthread.h>
#endif

#include <sys/ipc.h>
#include <sys/sem.h>
#if !defined(_SEM_SEMUN_UNDEFINED)
/* union semun is defined by including <sys/sem.h> */
#else
/* according to X/OPEN we have to define it ourselves */

union semun
{
    int val;                    /* value for SETVAL */
    struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
    unsigned short int *array;  /* array for GETALL, SETALL */
    struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif

#ifndef NO_GENERAL
#include "general.h"
#endif
//#include "filecls.h"

#include <math.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

#ifndef DEPTH
//#define DEPTH 16
#define DEPTH 32
#endif

#define DEPTH_REAL 8
#if DEPTH == 16
#define RGB RGB16
#define GETRED RED16
#define GETGREEN GREEN16
#define GETBLUE BLUE16
#elif DEPTH == 32
#define RGB RGB32
#define GETRED RED32
#define GETGREEN GREEN32
#define GETBLUE BLUE32
#endif

#define RGB32(R, G, B) ((unsigned)(R)<<16)+((unsigned)(G)<<8)+((unsigned)(B))
#define RGB16(R, G, B) ((unsigned)((R)&0xF8)<<8)+((unsigned)((G)&0xF8)<<3)+((unsigned)((B)>>3))

#define RGBany(R, G, B) ((((unsigned)(R)&CCT[0][0])<<CCT[0][1])>>CCT[0][2])\
                       +((((unsigned)(G)&CCT[1][0])<<CCT[1][1])>>CCT[1][2])\
                       +((((unsigned)(B)&CCT[2][0])<<CCT[2][1])>>CCT[2][2])

#define RED32(C)   (((C&0xFF0000)>>16)*1.0/(0x00FF))
#define GREEN32(C) (((C&0xFF00)>>8)*1.0/(0x00FF))
#define BLUE32(C)  (((C&0x00FF)*1.0/(0x00FF)))

#define RED16(C)   (((C&0xF800)>>11)*1.0/(0x001F))
#define GREEN16(C) (((C&0x07C0)>>6)*1.0/(0x001F))
#define BLUE16(C)  (((C&0x001F)*1.0/(0x001F)))

#define REDany(C)   (((C&CCT[0][3])>>CCT[0][4])*1.0/(CCT[0][5]))
#define GREENany(C) (((C&CCT[1][3])>>CCT[1][4])*1.0/(CCT[1][5]))
#define BLUEany(C)  (((C&CCT[2][3])>>CCT[2][4])*1.0/(CCT[2][5]))

/* smaller set */
//enum { MaxPoints=1000, MaxLines=1000, rInterval=50, DSCLEN=10};
//enum { MaxPoints=10000, MaxLines=5000, rInterval=50, DSCLEN=30};
enum { MaxPoints=40000, MaxLines=5000, rInterval=50, DSCLEN=30};
/* larger set */
//enum { MaxPoints=40000, MaxLines=80000, rInterval=50, DSCLEN=30};

inline void delay(long sec, long usec)
{
    //Define to use timeval without including <sys/time.h>
    struct {time_t tv_sec, tv_usec;} tv; 
    tv.tv_sec=sec;
    tv.tv_usec=usec;
    select(0,NULL,NULL,NULL,(timeval *)&tv);
}

extern unsigned long colorconvert_16[3][6], colorconvert_32[3][6];

class SYWindow//simple-display
{
public:
    enum{ ROTATION=0,TRANSLATION=1,SCALING=2,PROJECTION=3,
             PBCTRANSLATION=4,ASPECTRATIO=5,PBCX=6,PBCY=7,PBCZ=8,PBCGLIDE=9};
    struct YPoint
    {
        double x, y, z, r; //coordinates always lie in [-1, 1]
        unsigned long c;
        unsigned long attr;//0:disc 1:circle
        char dscrpt[DSCLEN];
    }Points[MaxPoints];
    struct YLine
    {
        double x0, y0, z0, x1, y1, z1, r;
        unsigned long c;
        unsigned long attr;//0:disc 1:circle
    }Lines[MaxLines];
    struct ZInd
    {
        double Z;
        int index;
    }Z_ind[MaxLines+MaxPoints];
    
    static int cmp(const void *p1, const void *p2);

    unsigned long CCT[3][6];
    
    int nP, nL; //Number of points/lines
    int semID; //for Lock/Unlock mechanism
    int semID2; //for Lock/Unlock mechanism
    bool alive,pause,drawframe;
    int loop;
    int rinterval;
    int msrsp;//mouse response

    KeySym ksfeedback; // key stroke feedback to main program

#ifdef USE_X
    
    Display *theDisplay;
    int theScreen;
    XVisualInfo visinfo;
    Visual *vis;
    XSetWindowAttributes attr;
    int PixDepth;
    GC gc;
    Pixmap pixmap;
    Window theWindow, Root;

#else

    int theWindow, Root;

    struct EventInfo {
      int mx, my;
      int mouse_pos_set;
      EventType type;
      int button_id;
      EventButtonState button_state;
      WindowAction action;
      void (*action_func)(int, int, int, int, int );
      } Event_info;

    struct Extent {
      float xmin, ymin, zmin;
      float xmax, ymax, zmax;
      } VpExtent;

    struct XformMap {
      float rot;
      float translate[3];
      float scale;
      } SceneXformMap;

    struct Xform {
      float rot[3], trans[3], scale[3], center[3];
      } SceneXform;

    struct Dlists {
      int disc_dlist, point_dlist, surf_dlist, line_dlist;
      } SceneDlist;

    int sphere_render; //different ways to render sphere

#endif

    unsigned int width, height;
    unsigned int depth;

#ifdef USE_X
    Colormap cmap;
#endif

    unsigned long cblack, cwhite, bgcolor, framecolor;

    double A11,A12,A13,A21,A22,A23,A31,A32,A33;//Rotation matrix.
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;//saved Rotation matrix.
    double B11,B12,B13,B21,B22,B23,B31,B32,B33;//Incremental rot. matrix.

    int R, B; //Current radius, box size(in pixels)
//    double xp,yp;//Last mouse position
    int Xp,Yp;
    int gifcount,pscount;
    unsigned long lastDrag;
    
    bool dirty; //true if need update
    int X0,Y0,X00,Y00;//offset
    double Scale,Scale0,Aspr,Aspr0;//scale and aspect ratio
    double D,D0; //projection
    bool square, sort;

    bool autowritegif;

    int scalepoints;
    int enable_pbc;
    double pbcshift[3];
    double maxpointradius, maxlinewidth;
    
    double CXs(double x, double z)
        { return (width/2+(B*x/(D==0?1:(1-z/D)))*Scale+X0); }
    double CYs(double y, double z)
        { return (height/2-(B*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
    double CXr(double x, double z)
        { return (width/2+(width/4*x/(D==0?1:(1-z/D)))*Scale+X0); }
    double CYr(double y, double z)
        { return (height/2-(height/4*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
    double CX(double x, double z) {return (square)?CXs(x,z):CXr(x,z); }
    double CY(double y, double z) {return (square)?CYs(y,z):CYr(y,z); }
    double CR(double r, double z) {if(scalepoints) return CRs(r,z); else return CRf(r,z); }
    double CRs(double r, double z) { return (B*r*Scale/(D==0?1:(1-z/D))); }
    double CRf(double r, double z) { return r; }
//    int CXs(double x, double z)
//        { return (int)(width/2+(B*x/(D==0?1:(1-z/D)))*Scale+X0); }
//    int CYs(double y, double z)
//        { return (int)(height/2-(B*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
//    int CXr(double x, double z)
//        { return (int)(width/2+(width/4*x/(D==0?1:(1-z/D)))*Scale+X0); }
//    int CYr(double y, double z)
//        { return (int)(height/2-(height/4*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
//    int CX(double x, double z) {return (square)?CXs(x,z):CXr(x,z); }
//    int CY(double y, double z) {return (square)?CYs(y,z):CYr(y,z); }
//    int CR(double r) { return (int)(B*r*Scale); }
    double DEG(double a) { return (M_PI*a/180); }

    void initRot()
    {
        a11=a22=a33=B11=B22=B33=1;
        a12=a13=a21=a23=a31=a32=B12=B13=B21=B23=B31=B32=0;
    }
    void copyRot()
    { A11=a11;A12=a12;A13=a13;
      A21=a21;A22=a22;A23=a23;
      A31=a31;A32=a32;A33=a33; }
    void saveRot()
    { a11=A11;a12=A12;a13=A13;
      a21=A21;a22=A22;a23=A23;
      a31=A31;a32=A32;a33=A33; }
    void saveScale()
    { Scale0=Scale;}
    void saveView()
    { saveRot();Scale0=Scale;X00=X0;Y00=Y0;D0=D;Aspr0=Aspr;}
        
    SYWindow(const SYWindow &){}
    const SYWindow &operator =(const SYWindow &yw){ return yw;}
public:
    SYWindow(int width_hint, int height_hint,
             const char *winname, bool s=true, bool so=false,bool fm=false);
    virtual ~SYWindow();
    void setinterval(int i) {if (i>0) rinterval=i; }
    bool IsAlive() { return alive; }
    bool TogglePause() { pause=!pause; return pause;}
    bool IsPaused() { return pause; }
#ifdef _NOLOCK
    void Lock()
    {
    }
    void Unlock()
    {
    }
    void InitSem()
    {
    }
    void InitSem2()
    {
    }
    void LockWritegif()
    {
    }
    void UnlockWritegif()
    {
    }
#else
#ifndef _MYOWNLOCK
    void Lock()
    {
        static sembuf sbl={0, -1, 0};
        semop(semID,&sbl,1);
//        printf("Lock semID=%d\n",semID);
    }
    void Unlock()
    {
        static sembuf sbu={0, 1, 0};
        semop(semID,&sbu,1);
//        printf("Unlock semID=%d\n",semID);
    }
    void InitSem()
    {
        semID=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
//        printf("Init semID=%d\n",semID);
        if(semID==-1) printf("semget failure! use -D_MYOWNLOCK\n");
    }
    void InitSem2()
    {
        semID2=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
//        printf("Init semID2=%d\n",semID2);
        if(semID2==-1) printf("semget failure! use -D_MYOWNLOCK\n");
    }
    void LockWritegif()
    {
        static sembuf sbl={0, -1, 0};
        autowritegif=true;
        semop(semID2,&sbl,1);
//        printf("Lock semID2=%d\n",semID2);
    }
    void UnlockWritegif()
    {
        static sembuf sbu={0, 1, 0};
        autowritegif=false;
        semop(semID2,&sbu,1);
    }
#else
    void Lock()
    {
//        printf("Lock semID=%d\n",semID);
        semID--;
        while(semID<0) sleep(1);
//        printf("Lock semID=%d\n",semID);
    }
    void Unlock()
    {
//        printf("Unlock semID=%d\n",semID);
        semID++;
//        printf("Unlock semID=%d\n",semID);
    }
    void InitSem()
    {
        printf("-D_MYOWNLOCK defined, use semaphore mimic\n");
        semID=0;
    }
    void InitSem2()
    {
        printf("-D_MYOWNLOCK defined, use semaphore mimic\n");
        semID2=0;
    }
    void LockWritegif()
    {
//        printf("Lock semID2=%d\n",semID2);
        semID2--;
        while(semID2<0) sleep(1);
    }
    void UnlockWritegif()
    {
//        printf("Unlock semID2=%d\n",semID2);
        semID2++;
    }
#endif
#endif
    void DrawPoint(double x, double y, double z, double r, unsigned long c,
                   unsigned long attr=0);
    void DrawPoint(double x, double y, double z, double r, unsigned long c,
                   char *ds, unsigned long attr=0);
    void DrawLine(double x0, double y0, double z0,
                  double x1, double y1, double z1, unsigned long c,
                  double r=0.01, unsigned long attr=0);
    void Draw3DLine(YLine);
    void Draw3DPixel(YPoint);
    void Draw3DLinetoPS(FILE *,YLine);
    void Draw3DPixeltoPS(FILE *,YPoint);
    unsigned long AllocNamedColor(char *name);
    unsigned long AllocRGBColor(unsigned r, unsigned g, unsigned b);
    unsigned long AllocShortRGBColor(unsigned r, unsigned g, unsigned b);
    unsigned long StdColor(int c)
    {
        return AllocShortRGBColor(c&1?0xff:0, c&2?0xff:0, c&4?0xff:0);
    }
    
    void Clear() { nP=nL=0; }
    void Refresh() { dirty=true; }

    virtual void paint();
    virtual void update();
    virtual void newGraph();
//    virtual void newGraph(bool init);
    virtual void Evolve();
    virtual int identify(int px, int py);

    virtual void Routine()
    {

#ifdef USE_X
        while(alive)
        {
            if(rinterval==0) delay(0,1000*rInterval);
            else delay(0,1000*rinterval);

            Evolve();
        }

#endif
    }

    virtual int ExtKeyHandler(KeySym ks)
    {
        ksfeedback=ks;
        if((int)ks)return 0;else return 1;
    };
    virtual void FreeResource()
    {
        union semun su={0};
        if(alive)
        {

#ifdef USE_X
            XFreeGC(theDisplay, gc);
            XFreePixmap(theDisplay, pixmap);
            XFreeColormap(theDisplay, cmap);
            XDestroyWindow(theDisplay, theWindow);
            XCloseDisplay(theDisplay);
#endif
            
            semctl(semID, 0, IPC_RMID, su);
            semctl(semID2, 0, IPC_RMID, su);
        }
        alive=false;
        pause=false;
        
    }    
#ifndef NO_THREAD
    static void *thread_routine(void *p) ;
    void Run() ;
#endif
    
#ifndef NO_XPM
    void writeXpm(char *name)
    {
        XpmWriteFileFromPixmap(theDisplay,name,pixmap,0,NULL);
    }
    void writegif();
#endif
    void importgif();
    void writeps();
    void testcolor();
    void reversergb();
};

class YWindow: public SYWindow
{
    //all about rotation
public:
    bool enableRot;
    void update();
    void applyRot();
    void horizontalRot(double);
    void verticalRot(double);
    void spinRot(double);
    void zoom(double);
    void rotateTo(XEvent);
    void translateTo(XEvent);
    void pbcshiftTo(XEvent,int);
    void pbcglideTo(XEvent);
    void pbcglide(double,double);
    void scaleTo(XEvent);
    void projectTo(XEvent);
    void applyRotate();
    void setWinSpec(int x0,int y0,double s,double d,double a[3][3]);
    YWindow(int w,int h,const char *n,bool s=true,bool so=true,bool fm=true):
        SYWindow(w,h,n,s,so,fm),enableRot(false){};
    void printWinSpec();
    virtual void Evolve(); //overload Evolve to handle more exceptions
    virtual void help();
    virtual void drawBoxFrame();
};



#endif // _DISPLAY_WIN_H

