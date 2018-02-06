/*
  display_win.cpp (Windows version of display.cpp)
  by David Park and Wei Cai caiwei@stanford.edu
  Last Modified : Sat Aug 13 16:35:58 2005

  FUNCTION  :  X-window display package

  Problem:  rotation matrix A can become non-unitary
*/

#include "display_win.h"

#ifdef TMCC
extern "C"
#endif

static int num_win = 0;
static int win_table[100][2];

#ifndef NO_THREAD

#ifdef TMCC
extern "C"
#endif

#define PBCSHIFT(a,shift,origin) \
    if(enable_pbc){a-=origin;a/=2;a+=shift;a-=rint(a);a*=2;a+=origin;}

void * SYWindow::thread_routine(void *p)
{
    int win;

    static float diffuse[] = {0.6, 0.6, 0.6, 1.0};
    //static float diffuse[] = {0.4, 0.4, 0.4, 1.0};
    
    //static float ambient[] = {0.2, 0.2, 0.2, 1.0};
    static float ambient[] = {0.4, 0.4, 0.4, 1.0};
    
    static float position[] = {0.0, 0.0, 1.0, 0.0};  /* infinite light location. */

    /**************
     ***  body  ***
     **************/
    
    //fprintf (stderr, "\n---------- thread_routine ---------- \n");

    /* set window start location.  */
    glutInitWindowPosition (100, 100);


    /* set width and height of window */
    glutInitWindowSize ((int)((SYWindow *)p)->width, (int)((SYWindow *)p)->height);


    /*  create glut window  */
    win = glutCreateWindow ("md++");
    glutDisplayFunc (display);


    /* set window reshape func.  */
    glutReshapeFunc (&event_WinReshape);

    /* keyboard event processing  */
    glutKeyboardFunc (event_KeyboardProc);

    /* special proc func (arrow key, functions , etc). */    
    glutSpecialFunc (&event_SpecialProc);

    /* set idle func.  */    
    glutIdleFunc (event_IdleProc);

    /* mouse button proc func. */
    glutMouseFunc (&event_Butp);

    /* mouse motion proc func. */    
    glutMotionFunc (&event_Motion);

    win_table[num_win][0] = glutGetWindow();
    win_table[num_win][1] = (int)p;
    num_win++;

    /* enable a single OpenGL light. */    
    glLightfv (GL_LIGHT0,  GL_DIFFUSE, diffuse);
    glLightfv (GL_LIGHT0, GL_POSITION, position);
    glLightfv (GL_LIGHT0, GL_AMBIENT, ambient);
    glEnable (GL_LIGHT0);
    glEnable (GL_LIGHTING);

    /* use depth buffering for hidden surface elimination. */
    glEnable (GL_DEPTH_TEST);
    
    /* ensures constant intensity when zooming  */    
    glEnable (GL_NORMALIZE);

    /* enable colors for polygons */
    glEnable (GL_COLOR_MATERIAL);
    
    viewport_Set (win, (int)((SYWindow *)p)->width, (int)((SYWindow *)p)->height);    

    /* process events  */
    glutMainLoop();

    return NULL;
}


void SYWindow::Run() //initiate a new thread to run Routine()
{
    int r;
    
    pthread_t thread;
    
    r = pthread_create (&thread, NULL, thread_routine, this);
    
#ifdef _GENERAL_H
    if(r!=0) FATAL("Fail to create thread: ["<<strerror(errno)<<"]");
#else
    if(r!=0) fprintf(stderr,"Fail to create thread\n");
#endif
}
#endif


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  Draw3DLine  ****                  *
 *                                                            *
 *------------------------------------------------------------*/

void SYWindow::Draw3DLine(YLine line)
{
    unsigned c, attr;
    double xs1, ys1, zs1, xs2, ys2, zs2;
    float rc, gc, bc;
    
    /**************
     ***  body  ***
     **************/
    
    //fprintf (stderr,"\n---------- Draw3DLine ---------- \n");
    //fflush (stderr);

    c   = line.c;
    attr = line.attr;
    xs1 = line.x0;
    ys1 = line.y0;
    zs1 = line.z0;
    
    PBCSHIFT(xs1, pbcshift[0], 0);
    PBCSHIFT(ys1,-pbcshift[1], 0);
    PBCSHIFT(zs1, pbcshift[2], 0);
    
    xs2 = line.x1;
    ys2 = line.y1;
    zs2 = line.z1;
    
    PBCSHIFT(xs2, pbcshift[0], xs1);
    PBCSHIFT(ys2,-pbcshift[1], ys1);
    PBCSHIFT(zs2, pbcshift[2], zs1);
    
    //fprintf (stderr,">>>>>> pt1 (%g %g %g) \n", xs1, ys1, zs1);
    //fprintf (stderr,">>>>>> pt2 (%g %g %g) \n", xs2, ys2, zs2);
    
    rc = ((c & 0xff0000) >> 16) / 255.0;
    gc = ((c & 0x00ff00) >> 8)  / 255.0;
    bc = ((c & 0x0000ff))       / 255.0;
    
    glBegin (GL_LINES);
    glColor3f (rc, gc, bc);
    glVertex3d (xs1, ys1, zs1);
    glVertex3d (xs2, ys2, zs2);
    glEnd();    
}


/*------------------------------------------------------------*
 *                                                            *
 *                   ****  Draw3DPixel ****                   *
 *                                                            *
 *------------------------------------------------------------*/

void SYWindow::Draw3DPixel (YPoint point)
{
    double x, y, r, xs, ys, zs; 
    unsigned c, attr, attr_real;
    int dlist;
    float rc, gc, bc;

    /**************
     ***  body  ***
     **************/
    
    //fprintf (stderr,"\n---------- Draw3DPixel ---------- \n");
    //fflush (stderr);
    
    xs = point.x;
    ys = point.y;
    zs = point.z;
    r  = point.r;
    c  = point.c;
    attr=point.attr;
    
    PBCSHIFT(xs,pbcshift[0],0);
    PBCSHIFT(ys,-pbcshift[1],0);
    PBCSHIFT(zs,pbcshift[2],0);
    
    if (r > maxpointradius) r=maxpointradius;
    
    if (((x<-r)||(x>(int)width+r))||((y<-r)||(y>(int)height+r))) return;

    if(sphere_render==0)
        attr_real = attr;
    else
        attr_real = sphere_render;
    
    if (attr_real == 1) 
        dlist = SceneDlist.disc_dlist;
    else
        dlist = SceneDlist.surf_dlist;

    //fprintf (stderr,">>>>>> r [%g] \n", r);
    
    glPushMatrix();
    glTranslated (xs, ys, zs);
    glScaled (r, r, r);

    glRotatef (-SceneXform.rot[2], 0.0, 0.0, 1.0);
    glRotatef (-SceneXform.rot[1], 0.0, 1.0, 0.0);
    glRotatef (-SceneXform.rot[0], 1.0, 0.0, 0.0);

    rc = ((c & 0xff0000) >> 16) / 255.0;
    gc = ((c & 0x00ff00) >> 8)  / 255.0;
    bc = ((c & 0x0000ff))       / 255.0;

    glColor3f (rc, gc, bc); 
    
    glCallList (dlist);

    glPopMatrix();
}


/*------------------------------------------------------------*
 *                                                            *
 *                      ****  paint  ****                     *
 *                                                            *
 *------------------------------------------------------------*/

void SYWindow::paint()
{

    int i, j, k, dlist, res, num_stack, ncp;
    float r, x, y, t, dt, cpts[100][2];
    float rc, gc, bc;

    /**************
     ***  body  ***
     **************/
    
    /*
      fprintf (stderr,"\n---------- SYWindow paint ---------- \n");
      fprintf (stderr,">>>>>> nL [%d] \n", nL);
      fprintf (stderr,">>>>>> nP [%d] \n", nP);
      fflush (stderr);
    */

    /* set background color */
    rc = ((bgcolor & 0xff0000) >> 16) / 255.0;
    gc = ((bgcolor & 0x00ff00) >> 8)  / 255.0;
    bc = ((bgcolor & 0x0000ff))       / 255.0;
    glClearColor (rc,gc,bc,1.0);
    
    glLoadIdentity ();
    glPushMatrix ();
    
    /* transform in reverse order  */
    
    glTranslatef (SceneXform.trans[0]+SceneXform.center[0], 
                  SceneXform.trans[1]+SceneXform.center[1],
                  SceneXform.trans[2]+SceneXform.center[2]);
    
    glRotatef (SceneXform.rot[0], 1.0, 0.0, 0.0);
    glRotatef (SceneXform.rot[1], 0.0, 1.0, 0.0);
    glRotatef (SceneXform.rot[2], 0.0, 0.0, 1.0);
    
    glScalef (SceneXform.scale[0], SceneXform.scale[0], SceneXform.scale[0]);
    
    glTranslatef (-SceneXform.center[0], -SceneXform.center[1], -SceneXform.center[2]);
    
    for (i = 0; i < nL; i++) {
        Draw3DLine (Lines[i]);
    }    
    
    /*  create sphere display lists  */    
    /* determines # of polygons that make up a sphere  */
    //res = 24; num_stack = 16; 
    res = 12; num_stack = 8;

    if (!SceneDlist.disc_dlist) {
        r = 1.0;
        ncp = 15;
        dt = (2.0 * M_PI) / (ncp - 1);
        t = 0.0;
        
        for (i = 0; i < ncp; i++) {
            cpts[i][0] = sin(t);
            cpts[i][1] = cos(t);
            t += dt;
        }
        
        dlist = glGenLists (1);
    
        glNewList (dlist, GL_COMPILE);
        
        glBegin (GL_TRIANGLES);
        glNormal3f(0.0, 0.0, 1.0);
        
        for (j = 0; j < ncp; j++)
        {
            if (j == ncp - 1) {
                k = 0;
            }
            else {
                k = j + 1;
            }
            
            x = r*cpts[j][0];
            y = r*cpts[j][1];
            glVertex3d (x, y, 0.0);
            x = r*cpts[k][0];
            y = r*cpts[k][1];
            glVertex3d (x, y, 0.0);
            glVertex3d (0.0, 0.0, 0.0);
        }
        
        glEnd ();
        
        glBegin (GL_LINES);
        glColor3f (0, 0, 0);
        
        for (j = 0; j < ncp; j++) {
            if (j == ncp - 1) {
                k = 0;
            }
            else {
                k = j + 1;
            }
            
            x = r*cpts[j][0];
            y = r*cpts[j][1];
            glVertex3d (x, y, 0.0);
            x = r*cpts[k][0];
            y = r*cpts[k][1];
            glVertex3d (x, y, 0.0);
        }
        
        glEnd ();
        glEndList();
        SceneDlist.disc_dlist = dlist;
    }
    
    if (!SceneDlist.surf_dlist) {
        r = 1.0;
        dlist = glGenLists (1);
        glNewList (dlist, GL_COMPILE);
        opengl_SurfGenSphere (SPHERE_REP_SOLID, r, res, num_stack);
        glEndList();
        SceneDlist.surf_dlist = dlist;
    }
    
    for (i = 0; i < nP;i++) {
        Draw3DPixel (Points[i]);
    }

    if (drawframe)
    {
        float frame[12][6] = { {-1,-1,-1,-1,-1, 1},
                               {-1, 1,-1,-1, 1, 1},
                               { 1,-1,-1, 1,-1, 1},
                               { 1, 1,-1, 1, 1, 1},
                               {-1,-1,-1,-1, 1,-1},
                               {-1,-1, 1,-1, 1, 1},
                               { 1,-1,-1, 1, 1,-1},
                               { 1,-1, 1, 1, 1, 1},
                               {-1,-1,-1, 1,-1,-1},
                               {-1,-1, 1, 1,-1, 1},
                               {-1, 1,-1, 1, 1,-1},
                               {-1, 1, 1, 1, 1, 1} };
        glColor3f (1.0, 1.0, 1.0);
        for(i=0;i<12;i++)
        {
            glBegin(GL_LINES);
            glVertex3d (frame[i][0],frame[i][1],frame[i][2]);
            glVertex3d (frame[i][3],frame[i][4],frame[i][5]);
            glEnd();
        }
    }
    
    glPopMatrix ();
}



int SYWindow::identify(int px, int py)
{
    int i;
    int hit;
    double xs, ys, zs, Z;
    hit=-1;
    for(i=0;i<nP;i++)
    {
        int x, y, r;
        if(Points[i].dscrpt[0]!=0)
        {
            xs = Points[i].x;
            ys = Points[i].y;
            zs = Points[i].z;
            
            PBCSHIFT(xs,pbcshift[0],0);
            PBCSHIFT(ys,-pbcshift[1],0);
            PBCSHIFT(zs,pbcshift[2],0);
    
            Z=A31*xs+A32*ys+A33*zs;
            x=(int)CX(A11*xs+A12*ys+A13*zs,Z);
            y=(int)CY(A21*xs+A22*ys+A23*zs,Z);
            r=(int)CR(Points[i].r,Z);
            if((px<x+r)&&(px>x-r)&&(py<y+r)&&(py>y-r))
            {
                hit=i;
#ifdef _GENERAL_H
                INFO("point("<<px<<","<<py<<"): "<<Points[i].dscrpt);
#else
                fprintf(stdout,"point (%d,%d): %s\n",px,py,Points[i].dscrpt);
#endif
            }
        }
    }
    return hit;
}


void SYWindow::update()
{
    fprintf (stderr,"\n---------- SYWindow update ---------- \n");
    fflush (stderr);
    
    if(!alive) return;
    Lock();
    
    paint();

    Unlock();
}


/*======================================================*
 *=============  o p e n g l   s t u f f  ==============*
 *======================================================*/

/*  stores button action functions.  */

static void (*action_func[BUTTON_MAX])(int, int, int, int, int);


/*------------------------------------------------------------*
 *                                                            *
 *                      ****  display  ****                   *
 *                                                            *
 * redisplay the graphics objects.                            *
 *------------------------------------------------------------*/

void display (void)
{

  int win;

  int wobj;

 /**************
  ***  body  ***
  **************/

  win = glutGetWindow();
  event_WinObjGet (win, &wobj);
  ((SYWindow *)wobj)->update(); 
}


/*------------------------------------------------------------*
 *                                                            *
 *                ****  viewport_Set  ****                    *
 *                                                            *
 * set opengl viewport. this is done when a window reshape    *
 * event occurs.                                              *
 *------------------------------------------------------------*/

void viewport_Set (int win, int w, int h)
{

  float xmin, xmax, ymin, ymax, zmin, zmax;

  float dz;

  int wobj;

 /**************
  ***  body  ***
  **************/

  /*
  fprintf (stderr, "\n---------- viewport_Set ---------- \n");
  fprintf (stderr, ">>>>>> win [%d]  \n", win);
  fprintf (stderr, ">>>>>> w [%d]  h [%d] \n", w, h);
  */

  event_WinObjGet (win, &wobj);

  xmin = ((SYWindow *)wobj)->VpExtent.xmin;
  ymin = ((SYWindow *)wobj)->VpExtent.ymin;
  zmin = ((SYWindow *)wobj)->VpExtent.zmin;
  xmax = ((SYWindow *)wobj)->VpExtent.xmax;
  ymax = ((SYWindow *)wobj)->VpExtent.ymax;
  zmax = ((SYWindow *)wobj)->VpExtent.zmax;

  /*
  fprintf (stderr, ">>>>>> xmin [%f]   xmax [%f] \n", xmin, xmax);
  fprintf (stderr, ">>>>>> ymin [%f]   ymax [%f] \n", ymin, ymax);
  */

  glViewport (0, 0, w, h);

  glMatrixMode (GL_PROJECTION);

  glLoadIdentity();

  dz = zmax - zmin;
  zmin -= 10*dz;
  zmax += 10*dz;

  glOrtho (xmin, xmax, ymin, ymax, zmin, zmax);

  glMatrixMode (GL_MODELVIEW);
}


/*------------------------------------------------------------*
 *                                                            *
 *              ****  opengl_SurfGenSphere  ****              *
 *                                                            *
 * create the geometry for a half a sphere. the resolution of *
 * the sphere is determined by 'slices' and 'stacks'.         *
 *------------------------------------------------------------*/

void opengl_SurfGenSphere (SphereRep dtype, float radius, int slices, int stacks)
{

  float rho, drho, theta, dtheta;

  float x, y, z;

  int i, j, imin, imax;

 /**************
  ***  body  ***
  **************/

  /*
  if (dtype == GR_GEOMETRY_DISPLAY_POINT) {
    slices = 6;
    stacks = 2;
    }
  */

  drho = M_PI / (float) stacks;
  dtheta = 2.0 * M_PI / (float) slices;

  imin = 1;
  imax = stacks / 2;

  if (dtype == SPHERE_REP_POINT) {
    glBegin (GL_POINTS);

    for (j = 0; j <= slices; j++) {
      theta = (j == slices) ? 0.0 : j * dtheta;
      x = (float)(-sin(theta) * sin(drho));
      y = (float)(cos(theta) * sin(drho));
      z = (float)(cos(drho));
      glVertex3f (x * radius, y * radius, z * radius);
      }

    for (i = imin; i < imax; i++) {
      rho = i * drho;

      for (j = 0; j <= slices; j++) {
        theta = (j == slices) ? 0.0f : j * dtheta;
        x = (float)(-sin(theta) * sin(rho));
        y = (float)(cos(theta) * sin(rho));
        z = (float)(cos(rho));
        glVertex3f (x * radius, y * radius, z * radius);
        }
      }

    glEnd();
    return;
  }

  if (dtype == SPHERE_REP_LINE) {
    for (i = imin; i < imax; i++) {
      rho = i * drho;
      glBegin (GL_LINE_LOOP);

      for (j = 0; j <= slices; j++) {
        theta = (j == slices) ? 0.0f : j * dtheta;
        x = (float)(-sin(theta) * sin(rho));
        y = (float)(cos(theta) * sin(rho));
        z = (float)(cos(rho));
        glVertex3f (x*radius, y*radius, z*radius);
        }

      glEnd();

      glBegin (GL_LINES);

      for (j = 0; j <= slices; j++) {
        theta = (j == slices) ? 0.0 : j * dtheta;
        x = (float)(-sin(theta) * sin(rho));
        y = (float)(cos(theta) * sin(rho));
        z = (float)(cos(rho));
        glVertex3f (0.0, 0.0, 0.0);
        glVertex3f (x*radius, y*radius, z*radius);
        }

      glEnd();
    }

    return;
  }


  /* draw end as a triangle fan */

  glBegin (GL_TRIANGLE_FAN);
  glNormal3f (0.0, 0.0, 1.0);
  glVertex3f (0.0, 0.0, radius);

  for (j = 0; j <= slices; j++) {
    theta = (j == slices) ? 0.0 : j * dtheta;
    x = (float)(-sin(theta) * sin(drho));
    y = (float)(cos(theta) * sin(drho));
    z = (float)(cos(drho));

    glNormal3f (x, y, z);
    glVertex3f (x * radius, y * radius, z * radius);
  }

  glEnd();

  
  /* draw intermediate stacks as quad strips  */

  for (i = imin; i < imax; i++) {
    rho = i * drho;

    glBegin (GL_QUAD_STRIP);

    for (j = 0; j <= slices; j++) {
      theta = (j == slices) ? 0.0f : j * dtheta;
      x = (float)(-sin(theta) * sin(rho));
      y = (float)(cos(theta) * sin(rho));
      z = (float)(cos(rho));
      glNormal3f(x, y, z);
      glVertex3f(x * radius, y * radius, z * radius);

      x = (GLfloat)(-sin(theta) * sin(rho + drho));
      y = (GLfloat)(cos(theta) * sin(rho + drho));
      z = (GLfloat)(cos(rho + drho));
      glNormal3f(x, y, z);
      glVertex3f(x * radius, y * radius, z * radius);
      }

    glEnd();
  }
}



/*============================================================*
 *============= e v e n t  p r o c e s s i n g  ==============*
 *============================================================*/

/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_WinObjGet  ****               *
 *                                                            *
 *------------------------------------------------------------*/

void event_WinObjGet (int win, int *p)
{
    
  int i;

 /**************
  ***  body  ***
  **************/

  for (i = 0; i < num_win; i++) {
    if (win == win_table[i][0]) { 
      *p = win_table[i][1];
      }
    }
}


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_WinReshape  ****              *
 *                                                            *
 * process a window reshape event.                            *
 *------------------------------------------------------------*/

void event_WinReshape (int width, int height)
{

  int win;

  float xmin, xmax, ymin, ymax, zmin, zmax;

  float cx, cy;

  float dx, dy;

  float f;

  int wobj;

 /**************
  ***  body  ***
  **************/

  //fprintf (stderr, "\n---------- event_WinReshape ---------- \n");
  win = glutGetWindow();
  event_WinObjGet (win, &wobj);

  xmin = ((SYWindow *)wobj)->VpExtent.xmin;
  ymin = ((SYWindow *)wobj)->VpExtent.ymin;
  zmin = ((SYWindow *)wobj)->VpExtent.zmin;
  xmax = ((SYWindow *)wobj)->VpExtent.xmax;
  ymax = ((SYWindow *)wobj)->VpExtent.ymax;
  zmax = ((SYWindow *)wobj)->VpExtent.zmax;

  cx = (xmax + xmin) / 2.0;
  cy = (ymax + ymin) / 2.0;
  dx = xmax - xmin;
  dy = ymax - ymin;


  /* rescale scene extent in order to
     retain aspect ratio.               */

  if (width <= height) {
    f = (float)height / (float)width;
    dy = dx * f;
    ymin = cy - dy / 2.0;
    ymax = cy + dy / 2.0;
    }
  else {
    f = (float)width / (float)height;
    dx = dy * f;
    xmin = cx - dx / 2.0;
    xmax = cx + dx / 2.0;
    }

  ((SYWindow *)wobj)->VpExtent.xmin = xmin;
  ((SYWindow *)wobj)->VpExtent.xmax = xmax;
  ((SYWindow *)wobj)->VpExtent.ymin = ymin;
  ((SYWindow *)wobj)->VpExtent.ymax = ymax;

  viewport_Set (win, width, height);

  ((SYWindow *)wobj)->update(); 
}


/*------------------------------------------------------------*
 *                                                            *
 *                ****  event_KeyboardProc  ****              *
 *                                                            *
 * process keyboard events.                                   *
 *------------------------------------------------------------*/

void event_KeyboardProc (unsigned char ch, int x, int y)
{
  int win, wobj;
 /**************
  ***  body  ***
  **************/

  if (ch == KEY_ESCAPE) {
    exit (0);
  }
  
  win = glutGetWindow();
  event_WinObjGet (win, &wobj);
  
//  fprintf(stderr,"key = %d\n",ch);

  switch (ch) {
    case 'p':
        ((SYWindow *)wobj)->TogglePause();
#ifdef _GENERAL_H
        INFO("Pause="<<((SYWindow *)wobj)->pause);return;
#else
        fprintf(stdout,"Pause=%d",(int)pause);return;
#endif
        break;

  case 'b':
      ((SYWindow *)wobj)->sphere_render++;
      ((SYWindow *)wobj)->sphere_render%=3;
#ifdef _GENERAL_H
        INFO("Sphere render scheme ="<<((SYWindow *)wobj)->sphere_render);
#else
        fprintf(stdout,"Sphere render scheme ="<<((SYWindow *)wobj)->sphere_render);
#endif
        display ();
        break;

  case 'h':      
      const char helpstr[]=
          "YWindow based on OpenGL:\n"
          "Mouse drag to rotate\n"
          "Hot Keys:\n"
          "F1        : display this message\n"
          "b         : change sphere rendering scheme\n"
          "F10       : output postscript\n";
#ifdef _GENERAL_H
      INFO("YWindow: "<<helpstr);
#else
      fprintf(stdout,"YWindow: %s\n",helpstr);
#endif
      break;

  }
}


/*------------------------------------------------------------*
 *                                                            *
 *               ****  event_SpecialProc  ****                *
 *                                                            *
 * process special function keys.                             *
 *------------------------------------------------------------*/

void event_SpecialProc (int key, int x, int y)
{
  int win, wobj;
 /**************
  ***  body  ***
  **************/

  win = glutGetWindow();
  event_WinObjGet (win, &wobj);
  
//  static char *key_name[] = {"left", "up", "right", "down"};
//  int i;
//  
//  if ((key >= KEY_ARROW_LEFT) && (key <= KEY_ARROW_DOWN)) {
//    i = key - KEY_ARROW_LEFT;
//    printf (" arrow = %s \n", key_name[i]);
//  }
//
//  fprintf(stderr,"special key = %d\n",key);

  switch (key)
  {
    case 10: /* F10 */
        ((SYWindow *)wobj)->writeps();
        fprintf(stderr,"not fully implemented yet\n");
    break;
  }
  
}


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_IdleProc  ****                *
 *                                                            *
 * idle process.                                              *
 *------------------------------------------------------------*/

void event_IdleProc()
  {

  int i;

  void *p;

 /**************
  ***  body  ***
  **************/

  for (i = 0; i < num_win; i++) {
    p = (void*)win_table[i][1];

    if (((SYWindow *)p)->dirty) {
      //fprintf (stderr, "\n>>>>>> idle dirty \n");
      ((SYWindow *)p)->update(); 
      ((SYWindow *)p)->dirty = false;
      }
    }
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_Butp  ****                    *
 *                                                            *
 * process mouse button press event.                          *
 *------------------------------------------------------------*/

void event_Butp (int button, int state, int x, int y)
{

  int mx, my;

  int key_mod, shift, ctrl;

  int win;

  EventButtonId button_id;

  int dx, dy;

  int wobj;

 /**************
  ***  body  ***
  **************/

#ifdef dbg_EventButp
  fprintf (stderr, "\n --------- EventButP --------- \n");
  fprintf (stderr, " >>>>>> button [%d] \n", button);
  fprintf (stderr, " >>>>>> state  [%d] \n", state);
#endif

  win = glutGetWindow();
  event_WinObjGet (win, &wobj);

  if (state == GLUT_UP) {
    event_Butr (x, y, button_id);
    return;
    }


  /*  determine key modifiers  */

  key_mod = glutGetModifiers();
  shift = (key_mod == GLUT_ACTIVE_SHIFT);
  ctrl = (key_mod == GLUT_ACTIVE_CTRL);

  switch (button) {
    case GLUT_LEFT_BUTTON:
      if (shift)
        button_id = BUTTON_SHIFT_LEFT;
      else if (ctrl)
        button_id = BUTTON_CTRL_LEFT;
      else
        button_id = BUTTON_LEFT;
    break;

    case GLUT_MIDDLE_BUTTON:
      if (shift)
        button_id = BUTTON_SHIFT_MIDDLE;
      else if (ctrl)
        button_id = BUTTON_CTRL_MIDDLE;
      else
        button_id = BUTTON_MIDDLE;
    break;

    case GLUT_RIGHT_BUTTON:
      if (shift)
        button_id = BUTTON_SHIFT_RIGHT;
      else if (ctrl)
        button_id = BUTTON_CTRL_RIGHT;
      else
        button_id = BUTTON_RIGHT;
    break;
    }

  mx = x;
  my = y;

  ((SYWindow *)wobj)->Event_info.type = EVENT_BUTTON_PRESS;
  ((SYWindow *)wobj)->Event_info.button_id = button_id;

  if (((SYWindow *)wobj)->Event_info.mouse_pos_set) {
    dx = mx - ((SYWindow *)wobj)->Event_info.mx;
    dy = my - ((SYWindow *)wobj)->Event_info.my;

    if ((button_id >= 0) && (button_id < BUTTON_MAX)) {
      (*action_func[button_id])(win, mx, my, dx, dy);
      }
    }

  else {
    ((SYWindow *)wobj)->Event_info.mouse_pos_set = 1;
    dx = 0;
    dy = 0;

    if ((button_id >= 0) && (button_id < BUTTON_MAX)) {
      (*action_func[button_id])(win, mx, my, dx, dy);
      }
    }

  ((SYWindow *)wobj)->Event_info.mx = mx;
  ((SYWindow *)wobj)->Event_info.my = my;
}


/*------------------------------------------------------------*
 *                                                            *
 *               ****  event_Motion  ****                     *
 *                                                            *
 * process mouse motion event. this will cause a rotation,    *
 * translation or scaling to be performed.                    *
 *------------------------------------------------------------*/

void event_Motion (int mx, int my)
{

  int win;

  int dx, dy;

  int wobj;

 /**************
  ***  body  ***
  **************/

#ifdef EventMotion
  fprintf (stderr, "\n --------- EventMotion --------- \n");
#endif

  win = glutGetWindow();
  event_WinObjGet (win, &wobj);
  ((SYWindow *)wobj)->Event_info.type = EVENT_MOUSE_MOTION;

  if (((SYWindow *)wobj)->Event_info.mouse_pos_set) {
    dx = mx - ((SYWindow*)wobj)->Event_info.mx;
    dy = my - ((SYWindow*)wobj)->Event_info.my;

    if (((SYWindow *)wobj)->Event_info.action != WINDOW_ACTION_NONE) {
      (*((SYWindow *)wobj)->Event_info.action_func)(win, mx, my, dx, dy);
      }
    }
  else {
    ((SYWindow *)wobj)->Event_info.mx = mx;
    ((SYWindow *)wobj)->Event_info.my = my;
    ((SYWindow *)wobj)->Event_info.mouse_pos_set = 1;
    }

  ((SYWindow *)wobj)->Event_info.mx = mx;
  ((SYWindow *)wobj)->Event_info.my = my;
}


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_Butr  ****                    *
 *                                                            *
 * process mouse button release event.                        *
 *------------------------------------------------------------*/

void event_Butr (int mx, int my, EventButtonId button_id)
{

  int win, wobj;

 /**************
  ***  body  ***
  **************/

  //fprintf (stderr, "---------- EventButr ---------- \n");

  win = glutGetWindow();
  event_WinObjGet (win, &wobj);

  ((SYWindow *)wobj)->Event_info.mx = mx;
  ((SYWindow *)wobj)->Event_info.my = my;
  ((SYWindow *)wobj)->Event_info.mouse_pos_set = 0;
  ((SYWindow *)wobj)->Event_info.type = EVENT_NONE;
  ((SYWindow *)wobj)->Event_info.action = WINDOW_ACTION_NONE;
}


/*------------------------------------------------------------*
 *                                                            *
 *              ****  event_WinActionRotxy  ****              *
 *                                                            *
 * process rotation about x-y axes.                           *
 *------------------------------------------------------------*/

void event_WinActionRotxy (int win, int mx, int my, int dx, int dy)
{

  float rx, ry;

  int wobj;

 /**************
  ***  body  ***
  **************/

  event_WinObjGet (win, &wobj);

  //fprintf (stderr, "---------- event_WinActionRotxy -------- \n");
  //fprintf (stderr, " dx dy (%d %d) \n", dx, dy);

  ((SYWindow *)wobj)->Event_info.action = WINDOW_ACTION_ROTATE_XY;
  ((SYWindow *)wobj)->Event_info.action_func = event_WinActionRotxy;

  if (((SYWindow *)wobj)->Event_info.type == EVENT_BUTTON_PRESS) {
    return;
    }

  rx = (float)dy * ((SYWindow *)wobj)->SceneXformMap.rot;
  ry = (float)dx * ((SYWindow *)wobj)->SceneXformMap.rot;

  if ((fabs(rx) < 1.0) && (fabs(ry) < 1.0)) {
    return;
    }

  ((SYWindow *)wobj)->SceneXform.rot[0] += rx;
  ((SYWindow *)wobj)->SceneXform.rot[1] += ry;

  display ();
}


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_WinActionRotz  ****           *
 *                                                            *
 * process rotation about z axes.                             *
 *------------------------------------------------------------*/

void event_WinActionRotz (int win, int mx, int my, int dx, int dy)
{

  int wobj;

  float rz;

 /**************
  ***  body  ***
  **************/

  //fprintf (stderr, "---------- event_WinActionRotz -------- \n");
  //fprintf (stderr, " dx dy (%d %d) \n", dx, dy);
  event_WinObjGet (win, &wobj);

  ((SYWindow *)wobj)->Event_info.action = WINDOW_ACTION_ROTATE_Z;
  ((SYWindow *)wobj)->Event_info.action_func = event_WinActionRotz;

  rz = (float)(dx + dy) * ((SYWindow *)wobj)->SceneXformMap.rot;
  ((SYWindow *)wobj)->SceneXform.rot[2] += rz;

  display ();
}


/*------------------------------------------------------------*
 *                                                            *
 *              ****  event_WinActionPick  ****               *
 *                                                            *
 * process picking.                                           *
 *------------------------------------------------------------*/

void event_WinActionPick (int win, int mx, int my, int dx, int dy)
{

 /**************
  ***  body  ***
  **************/

}


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_WinActionTransxy  ****        *
 *                                                            *
 * process a translate in x-y.                                *
 *------------------------------------------------------------*/

void event_WinActionTransxy (int win, int mx, int my, int dx, int dy)
{

  int wobj;

  float tx, ty;

 /**************
  ***  body  ***
  **************/

  event_WinObjGet (win, &wobj);

  ((SYWindow *)wobj)->Event_info.action = WINDOW_ACTION_TRANSLATE_XY;
  ((SYWindow *)wobj)->Event_info.action_func = event_WinActionTransxy;

  tx = (float)dx  * ((SYWindow *)wobj)->SceneXformMap.translate[0];
  ty = (float)-dy * ((SYWindow *)wobj)->SceneXformMap.translate[1];

  ((SYWindow *)wobj)->SceneXform.trans[0] += tx;
  ((SYWindow *)wobj)->SceneXform.trans[1] += ty;

  display ();
}


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_WinActionTransz  ****         *
 *                                                            *
 * process a translate in z.                                  *
 *------------------------------------------------------------*/

void event_WinActionTransz (int win, int mx, int my, int dx, int dy)
{

  int wobj;

  float tz;

 /**************
  ***  body  ***
  **************/

  event_WinObjGet (win, &wobj);

  ((SYWindow *)wobj)->Event_info.action = WINDOW_ACTION_TRANSLATE_Z;
  ((SYWindow *)wobj)->Event_info.action_func = event_WinActionTransz;

  tz = (float)(dx + dy) * ((SYWindow *)wobj)->SceneXformMap.translate[2];

  ((SYWindow *)wobj)->SceneXform.trans[2] += tz;

  display ();
}

/*------------------------------------------------------------*
 *                                                            *
 *                  ****  event_WinActionScale  ****          *
 *                                                            *
 * process a scaling action.                                  *
 *------------------------------------------------------------*/

void event_WinActionScale (int win, int mx, int my, int dx, int dy)
{

  int wobj;

  float scale;

  float d;

 /**************
  ***  body  ***
  **************/

  //fprintf (stderr, "---------- event_WinActionScale -------- \n");
  //fprintf (stderr, " dx dy (%d %d) \n", dx, dy);

  event_WinObjGet (win, &wobj);

  ((SYWindow *)wobj)->Event_info.action = WINDOW_ACTION_SCALE;
  ((SYWindow *)wobj)->Event_info.action_func = event_WinActionScale;
  d = (float)(dx + dy);

  if (d < 0.0) {
    scale = 1.0 / ((SYWindow *)wobj)->SceneXformMap.scale;
    }
  else if (d == 0.0) {
    scale = 1.0;
    }
  else {
    scale = ((SYWindow *)wobj)->SceneXformMap.scale;
    }

  ((SYWindow *)wobj)->SceneXform.scale[0] *= scale;
  ((SYWindow *)wobj)->SceneXform.scale[1] *= scale;
  ((SYWindow *)wobj)->SceneXform.scale[2] *= scale;

  display ();
}



void SYWindow::newGraph()
{
    if (!alive) {
        alive = true;
    }

    dirty = true;
}

void SYWindow::Evolve()
{
    fprintf (stderr,"\n----------SYWindow evolve ---------- \n");
    fflush (stderr);

    if(!alive) return;

    if(dirty)
    {
        dirty=false;
        update();
    }
}

void SYWindow::DrawPoint(double x, double y, double z, double r,
                        unsigned long c, unsigned long attr)
{
    if(nP<MaxPoints-1)
    {
        Points[nP].x=x;
        Points[nP].y=y;
        Points[nP].z=z;
        Points[nP].r=r;
        Points[nP].c=c;
        Points[nP].attr=attr;
        Points[nP].dscrpt[0]=0;
        nP++;
    }
    else if(nP==MaxPoints-1)
    {
#ifdef _GENERAL_H
        WARNING("SYWindow: Too many points (more than "
                 << (int)MaxPoints
                 << "), ignoring");
#else
        fprintf(stderr,"SYWindow: Too many points (more than %d"
                 "), ignoring\n\n", (int)MaxPoints );
#endif
        Points[nP].x=x;
        Points[nP].y=y;
        Points[nP].z=z;
        Points[nP].r=r;
        Points[nP].c=c;
        Points[nP].attr=attr;
        Points[nP].dscrpt[0]=0;
        nP++;
    }
}

void SYWindow::DrawPoint(double x, double y, double z, double r,
                        unsigned long c, char *ds, unsigned long attr)
{
    if(nP<MaxPoints-1)
    {
        Points[nP].x=x;
        Points[nP].y=y;
        Points[nP].z=z;
        Points[nP].r=r;
        Points[nP].c=c;
        Points[nP].attr=attr;
        strncpy(Points[nP].dscrpt,ds,DSCLEN-1);
        nP++;
    }
    else if (nP==MaxPoints-1)
    {
#ifdef _GENERAL_H
        WARNING("SYWindow: Too many points (more than "
                 << (int)MaxPoints
                 << "), ignoring");
#else
        fprintf(stderr,"SYWindow: Too many points (more than %d"
                 "), ignoring", (int)MaxPoints );
#endif
        Points[nP].x=x;
        Points[nP].y=y;
        Points[nP].z=z;
        Points[nP].r=r;
        Points[nP].c=c;
        Points[nP].attr=attr;
        strncpy(Points[nP].dscrpt,ds,DSCLEN-1);
        nP++;
    }

}

void SYWindow::DrawLine(double x0, double y0, double z0,
                        double x1, double y1, double z1, unsigned long c,
                        double r, unsigned long attr)
{
    if(nL<MaxLines-1)
    {
        Lines[nL].x0=x0;
        Lines[nL].y0=y0;
        Lines[nL].z0=z0;
        Lines[nL].x1=x1;
        Lines[nL].y1=y1;
        Lines[nL].z1=z1;
        Lines[nL].c=c;
        Lines[nL].r=r;
        Lines[nL].attr=attr;
        nL++;
    }
    else if(nL==MaxLines-1)
    {
#ifdef _GENERAL_H
        WARNING("SYWindow: Too many lines (more than "
                 << (int)MaxLines
                 << "), ignoring");
#else
        fprintf(stderr,"SYWindow: Too many lines (more than %d"
                 "), ignoring\n\n", (int)MaxLines );
#endif
        Lines[nL].x0=x0;
        Lines[nL].y0=y0;
        Lines[nL].z0=z0;
        Lines[nL].x1=x1;
        Lines[nL].y1=y1;
        Lines[nL].z1=z1;
        Lines[nL].c=c;
        Lines[nL].r=r;
        Lines[nL].attr=attr;
        nL++;
    }
}

SYWindow::SYWindow(int width_hint, int height_hint,
                   const char *winname, bool s, bool so, bool fm)
{
#ifdef USE_X
    XSizeHints myhint;
    unsigned long mask; 
    Cursor arrow;
#endif

//    int i,j;

  //fprintf (stderr,"\n---------- SYWindow initialize ---------- \n");
  //fprintf (stderr,">>>>>> width_hint [%d] \n", width_hint);
  //fflush (stderr);
    
    //fprintf(stderr,"SYWindow initializer");fflush(stderr);
    initRot();copyRot();

  int argc = 1;

  char *argv = "md++";

  float xmin, xmax, ymin, ymax, zmin, zmax;

  float d, dx, dy, dz;

  static float f = 0.001;


  glutInit (&argc, &argv);

  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  width  = width_hint;
  height = height_hint;
  ksfeedback = 0; // key stroke feedback to main program
    

  /* initialize event info */

  Event_info.mouse_pos_set = 0;
  Event_info.type = EVENT_NONE;
  Event_info.button_id = BUTTON_UNKNOWN;
  Event_info.button_state = EVENT_BUTTON_STATE_NONE;
  Event_info.action = WINDOW_ACTION_NONE;
  Event_info.mx = 0;
  Event_info.my = 0;


  /* set scene extent  */

  d = 1.5;
  VpExtent.xmin = -d;
  VpExtent.ymin = -d;
  VpExtent.zmin = -10*d;
  VpExtent.xmax =  d;
  VpExtent.ymax =  d;
  VpExtent.zmax =  10*d;

  SceneXform.rot[0] = 0.0;
  SceneXform.rot[1] = 0.0;
  SceneXform.rot[2] = 0.0;

  SceneXform.trans[0] = 0.0;
  SceneXform.trans[1] = 0.0;
  SceneXform.trans[2] = 0.0;

  SceneXform.scale[0] = 1.0;
  SceneXform.scale[1] = 1.0;
  SceneXform.scale[2] = 1.0;

  SceneXform.center[0] = 0.0;
  SceneXform.center[1] = 0.0;
  SceneXform.center[2] = 0.0;

  SceneDlist.disc_dlist = 0;
  SceneDlist.point_dlist = 0;
  SceneDlist.line_dlist = 0;
  SceneDlist.surf_dlist = 0;


  /*  map mouse buttons to scene transformations  */

  action_func[BUTTON_LEFT] = event_WinActionRotxy;
  action_func[BUTTON_SHIFT_LEFT] = event_WinActionRotz;
  action_func[BUTTON_CTRL_LEFT] = event_WinActionPick;
  action_func[BUTTON_RIGHT] = event_WinActionTransxy;
  action_func[BUTTON_SHIFT_RIGHT] = event_WinActionTransz;
  action_func[BUTTON_CTRL_RIGHT] = event_WinActionScale;


  /*  set mapping of device units to
      mouse transformations.           */

  SceneXformMap.rot = 0.5;
  SceneXformMap.scale = 1.1;

  xmin = VpExtent.xmin;
  ymin = VpExtent.xmin;
  zmin = VpExtent.zmin;
  xmax = VpExtent.xmax;
  ymax = VpExtent.xmax;
  zmax = VpExtent.zmax;

  dx = xmax - xmin;
  dy = ymax - ymin;
  dz = zmax - zmin;

  SceneXformMap.translate[0] = f * dx;
  SceneXformMap.translate[1] = f * dy;
  SceneXformMap.translate[2] = f * dz;

  sphere_render = 0;
  
    alive=false;
    pause=false;

    autowritegif=false;
    
    framecolor = cwhite;

    newGraph();


    nL=nP=0;
    
    //For lock/unlock
//    semID=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
//    printf("Init semID=%d\n",semID);
//    if(semID==-1) printf("semget failure!\n");
    InitSem();
    Unlock();
    //For lock/unlock
//    semID2=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
//    printf("Init semID2=%d\n",semID2);
//    if(semID2==-1) printf("semget failure!\n");
    InitSem2();
//    UnlockWritegif();
    //set bool square
    square = s;
    sort = so;
    drawframe = fm;
    X0=Y0=X00=Y00=0;Scale=Scale0=1;Aspr=Aspr0=1;
    //projection
    D=D0=10000;

    gifcount=0;pscount=0;

//    lastDrag = 0xFFFFFFFF;
    lastDrag = (unsigned long) (-1);

    rinterval=0;
    msrsp=ROTATION;

    scalepoints = 1;
    enable_pbc = 0;
    pbcshift[0]=pbcshift[1]=pbcshift[2]=0 ;
    maxpointradius = 100;
    maxlinewidth = 100;
}

SYWindow::~SYWindow()
{
//    if(alive)
//    {
//        //alive=false;
//        XFreeGC(theDisplay, gc);
//        XFreePixmap(theDisplay, pixmap);
//    }        
//    XFreeColormap(theDisplay, cmap);
//    XDestroyWindow(theDisplay, theWindow);
//    XCloseDisplay(theDisplay);

    FreeResource();
}

void SYWindow::reversergb()
{
    int j;
    unsigned long L;
    
    for(j=0;j<6;j++)
    {
        L=CCT[0][j];
        CCT[0][j]=CCT[2][j];
        CCT[2][j]=L;
    }
}


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  AllocRGBColor  ****               *
 *                                                            *
 *------------------------------------------------------------*/

unsigned long 
SYWindow::AllocRGBColor (unsigned r, unsigned g, unsigned b)
  {

#ifdef USE_X
    XColor c;
    c.red=r, c.green=g, c.blue=b;
//    unsigned long p;double red,green,blue;
#ifdef _GENERAL_H
    if(XAllocColor(theDisplay, cmap, &c)==0)
        WARNING("Error allocating color ("<<r<<", "<<g<<", "<<b<<")");
#else
    if(XAllocColor(theDisplay, cmap, &c)==0)
        fprintf(stderr,"Error allocating color (%d,%d,%d)\n",r,g,b);
#endif

//    p=c.pixel;
//    r=c.red;g=c.green;b=c.blue;
//    red=((p&0x7C00)>>10)*1.0/(0x001F);
//    green=((p&0x03E0)>>5)*1.0/(0x001F);
//    blue=(p&0x001F)*1.0/(0x001F);
//    INFO_Printf("c.pixel=%x c.red=%x(%f) c.green=%x(%f) c.blue=%x(%f)\n",
//                p,r,red,g,green,b,blue);

    return c.pixel;

#else

  unsigned long val;
  val = b + 256 * (g + 256*r);

  fprintf (stderr, " SYWindow::AllocRGBColor: rgb (%d %d %d) \n", r, g, b);
  return val;
  
#endif
  
}


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  AllocShortRGBColor  ****          *
 *                                                            *
 *------------------------------------------------------------*/

unsigned long SYWindow::AllocShortRGBColor(unsigned r, unsigned g, unsigned b)
{
    unsigned long val;
    val = b + 256 * (g + 256*r);
    return val;
}


/*------------------------------------------------------------*
 *                                                            *
 *              ****  SYWindow::AllocNamedColor  ****         *
 *                                                            *
 *------------------------------------------------------------*/

unsigned long SYWindow::AllocNamedColor(char *name)
{
    fprintf(stderr, "SYWindow::AllocNamedColor not supported in Windows version\n"); 
    return 0;
}

void SYWindow::testcolor()
{
    /* not valid for Windows version */
}

void SYWindow::importgif()
{
    fprintf(stderr, "SYWindow::importgif not supported in Windows version\n"); 
}

void SYWindow::Draw3DLinetoPS(FILE *file,YLine line)
{
    double xs1,ys1,zs1,xs2,ys2,zs2,Z,r;unsigned attr;
    double red,green,blue;
    double x1, y1, x2, y2;
    attr = line.attr;
    xs1 = line.x0;
    ys1 = line.y0;
    zs1 = line.z0;

    PBCSHIFT(xs1,pbcshift[0],0);
    PBCSHIFT(ys1,-pbcshift[1],0);
    PBCSHIFT(zs1,pbcshift[2],0);
    
    Z=A31*xs1+A32*ys1+A33*zs1;
    x1=CX(A11*xs1+A12*ys1+A13*zs1,Z);
    y1=CY(A21*xs1+A22*ys1+A23*zs1,Z);
    xs2 = line.x1;
    ys2 = line.y1;
    zs2 = line.z1;

    PBCSHIFT(xs2,pbcshift[0],xs1);
    PBCSHIFT(ys2,-pbcshift[1],ys1);
    PBCSHIFT(zs2,pbcshift[2],zs1);
        
    Z=A31*xs2+A32*ys2+A33*zs2;
    x2=CX(A11*xs2+A12*ys2+A13*zs2,Z);
    y2=CY(A21*xs2+A22*ys2+A23*zs2,Z);
    r=CR(line.r,0);
    if(r>maxlinewidth) r=maxlinewidth;
    if(r<0.5) r=0.5;
    
    //need some criteria to avoid waste
    //temporary
    if((((x1<0)||(x1>(int)width))||((y1<0)||(y1>(int)height)))
       &&(((x2<0)||(x2>(int)width))||((y2<0)||(y2>(int)height))))
        return;
//    red=GETRED(line.c);
//    green=GETGREEN(line.c);
//    blue=GETBLUE(line.c);

    red=REDany(line.c);
    green=GREENany(line.c);
    blue=BLUEany(line.c);
    
    fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb"
            " %5.2f %5.2f %5.2f %5.2f m l s\n",
            r,red,green,blue,
            x1,height-y1,x2,height-y2);
    if(attr!=0)
    {
        double dy,dx,dr;
        dy=y2-y1;dx=x2-x1;dr=sqrt(dy*dy+dx*dx);dy/=dr;dx/=dr;
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb"
                " %5.2f %5.2f %5.2f %5.2f m l s\n",
                0.5,0.,0.,0.,
                x1+dy*r/2,height-(y1-dx*r/2),
                x2+dy*r/2,height-(y2-dx*r/2));
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb"
                " %5.2f %5.2f %5.2f %5.2f m l s\n",
                0.5,0.,0.,0.,
                x1-dy*r/2,height-(y1+dx*r/2),
                x2-dy*r/2,height-(y2+dx*r/2));
    }
}

void SYWindow::Draw3DPixeltoPS(FILE *file,YPoint point)
{
    double x,y,r; 
//#define _ALLCIRCLE
#ifndef _ALLCIRCLE
    unsigned attr;
#endif
    double xs, ys, zs, Z;
    double red,green,blue;
    xs = point.x;
    ys = point.y;
    zs = point.z;

    PBCSHIFT(xs,pbcshift[0],0);
    PBCSHIFT(ys,-pbcshift[1],0);
    PBCSHIFT(zs,pbcshift[2],0);
    
#ifndef _ALLCIRCLE
    attr = point.attr;
#endif
    Z=A31*xs+A32*ys+A33*zs;
    x=CX(A11*xs+A12*ys+A13*zs,Z);
    y=CY(A21*xs+A22*ys+A23*zs,Z);
    r=CR(point.r,Z);
    if(((x<-r)||(x>(int)width+r))||((y<-r)||(y>(int)height+r))) return;
//    red=GETRED(point.c);
//    green=GETGREEN(point.c);
//    blue=GETBLUE(point.c);
    red=REDany(point.c);
    green=GREENany(point.c);
    blue=BLUEany(point.c);

    if(r>maxpointradius) r=maxpointradius;

//#if DEPTH_REAL == 8
//    red=green=blue=1;
//#endif
//#ifdef _ALLCIRCLE
//    if (0)
//#else
//    if (attr==0)
//#endif
    if(attr==0)
    { /* black border with white disc */
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                1.,1.,1.,
                x,height-y,r);
        fprintf(file,"disc\n");
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                0.,0.,0.,x,height-y,r);
        fprintf(file,"circle\n");
    }
    else if(attr==2)
    { /* black border with colored disc */
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                red,green,blue,
                x,height-y,r);
        fprintf(file,"disc\n");
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                0.,0.,0.,x,height-y,r);
        fprintf(file,"circle\n");
    }
    else
    { /* colored circle with empty disc */
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                red,green,blue,
                x,height-y,r);
        fprintf(file,"circle\n");
    }
}

void SYWindow::writeps()
{
    int i;
    
    const char fname[]="Yshot";
    char extname[100],tmp[100];
    char f1[100];
    FILE *file;

    char head1[500]="%!PS-Adobe-3.0  EPSF-3.0\n"
        "%%Pages: (atend)\n";
    char head2[500]="%%BoundingBox:";
    char head3[5000]=
        "%%EndComments\n"
        "/l {lineto} def\n"
        "/m {moveto} def\n"
        "/t {translate} def\n"
        "/slw {setlinewidth} def\n"
        "/srgb {setrgbcolor} def\n"
        "/np {newpath} def\n"
        "/s {stroke} def\n"
        "/disc { 0 360 arc fill } def\n"
        "/circle { 0 360 arc s } def\n";
        
    char tail[500]="showpage\n";

    sprintf(tmp,"%04d",pscount);
    strcpy(extname,fname);
    strcat(extname,tmp);
    sprintf(f1,"%s.ps",extname); //filename 1
    file=fopen(f1,"w");
    pscount++;

    fprintf(file,"%s",head1);
#ifdef _GENERAL_H
    INFO("writeps -> "<<f1);
    INFO("write to file");
#else
    fprintf(stdout,"writeps -> %s\nwrite to file\n",f1);
#endif
    fprintf(file,"%s %d %d %d %d\n",head2,0,0,width,height);
    fprintf(file,"%s",head3);

    if (sort)
    {
        for(i=0;i<nL+nP;i++)
        {
            int n=Z_ind[i].index;
            if(n<nL) Draw3DLinetoPS(file,Lines[n]);
            else Draw3DPixeltoPS(file,Points[n-nL]);
        }
    }
    else
    {
        for(i=0;i<nL;i++)
        {
            Draw3DLinetoPS(file,Lines[i]);
        }
        for(i=0;i<nP;i++)
        {
            Draw3DPixeltoPS(file,Points[i]);
        }
    }        

    fprintf(file,"%s",tail);
    fclose(file);
}

void YWindow::help()
{
    /* moved to GL handling */
}

void YWindow::printWinSpec()
{
    char tmp[1000];
    sprintf(tmp,"YWindow: Window Specification\n"
            "width=%d\theight=%d\n"
            "X0=%d\tY0=%d\n"
            "Scale=%f\tD=%f\n"
            "rotation matrix=[%10f %10f %10f\n"
            "                 %10f %10f %10f\n"
            "                 %10f %10f %10f]\n\n"
            ,width,height,X0,Y0,Scale,D
            ,A11,A12,A13,A21,A22,A23,A31,A32,A33
            );
#ifdef _GENERAL_H
    INFO(tmp);
#else
    fprintf(stdout,"%s\n",tmp);
#endif
}

void YWindow::setWinSpec(int x0,int y0,double s,double d,double a[3][3])
{
    X0=x0;Y0=y0;Scale=s;D=d;
    A11=a[0][0];A12=a[0][1];A13=a[0][2];
    A21=a[1][0];A22=a[1][1];A23=a[1][2];
    A31=a[2][0];A32=a[2][1];A33=a[2][2];
}


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  update  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

void YWindow::update()
{
    
    /**************
     ***  body  ***
     **************/
    
    /*
      fprintf (stderr,"\n---------- Ywindow update ---------- \n");
      fflush (stderr);
    */
    
    Lock();
    
    /* clear z-buffer and color bits  */
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    paint();
    
    glutSwapBuffers();

    Unlock();
}


void YWindow::drawBoxFrame()
{
    /* moved into GL paint */
}    


//Basic Rotations----------------------------------------
void YWindow::applyRotate()
{
    double d1,d2,d3;
    d1=B11*A11+B12*A21+B13*A31;
    d2=B21*A11+B22*A21+B23*A31;
    d3=B31*A11+B32*A21+B33*A31;
    A11=d1; A21=d2; A31=d3;
    d1=B11*A12+B12*A22+B13*A32;
    d2=B21*A12+B22*A22+B23*A32;
    d3=B31*A12+B32*A22+B33*A32;
    A12=d1; A22=d2; A32=d3;
    d1=B11*A13+B12*A23+B13*A33;
    d2=B21*A13+B22*A23+B23*A33;
    d3=B31*A13+B32*A23+B33*A33;
    A13=d1; A23=d2; A33=d3;
    dirty=true;
}

void YWindow::horizontalRot(double arc)
{
    double d1,d3;
    double c, s;
    c=cos(arc);
    s=sin(arc);
    d1=c*A11+s*A31;
    d3=-s*A11+c*A31;
    A11=d1;A31=d3;
    d1=c*A12+s*A32;
    d3=-s*A12+c*A32;
    A12=d1;A32=d3;
    d1=c*A13+s*A33;
    d3=-s*A13+c*A33;
    A13=d1;A33=d3;
    dirty=true;
}

void YWindow::verticalRot(double arc)
{
    double d2,d3;
    double c,s;
    s=sin(arc);
    c=cos(arc);
    d2=c*A21+s*A31;
    d3=-s*A21+c*A31;
    A21=d2;A31=d3;
    d2=c*A22+s*A32;
    d3=-s*A22+c*A32;
    A22=d2;A32=d3;
    d2=c*A23+s*A33;
    d3=-s*A23+c*A33;
    A23=d2;A33=d3;
    dirty=true;
}

void YWindow::spinRot(double arc)
{
    double d1,d2;
    double c,s;
    s=sin(arc);
    c=cos(arc);
    d1=c*A11+s*A21;
    d2=-s*A11+c*A21;
    A11=d1;A21=d2;
    d1=c*A12+s*A22;
    d2=-s*A12+c*A22;
    A12=d1;A22=d2;
    d1=c*A13+s*A23;
    d2=-s*A13+c*A23;
    A13=d1;A23=d2;
    dirty=true;
}

void YWindow::zoom(double z)
{
//    INFO("Scale="<<Scale<<"  z="<<z);
    Scale*=z;
//    INFO("  Scale="<<Scale);
}

void YWindow::scaleTo(XEvent ev)
{
    Scale*=exp((ev.xmotion.x-Xp+ev.xmotion.y-Yp)*0.001);
    if(Scale<5e-2) Scale=5e-2;
    if(Scale>1e2) Scale=1e2;
//    DUMP("YWindow: Scale="<<Scale);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::translateTo(XEvent ev)
{
    X0+=ev.xmotion.x-Xp;Y0+=ev.xmotion.y-Yp;
//    DUMP("YWindow: X0="<<X0<<" ,Y0="<<Y0);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::pbcshiftTo(XEvent ev, int dir)
{
    pbcshift[dir]+=(0.0+ev.xmotion.x-Xp+ev.xmotion.y-Yp)/(2.0*B*Scale);
//    printf("pbcshift[%d]=%e, B=%d\n",dir,pbcshift[dir],B);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::pbcglideTo(XEvent ev)
{
    pbcglide((ev.xmotion.x-Xp)/(2.0*B),(ev.xmotion.y-Yp)/(2.0*B));
//    printf("pbcshift[%d]=%e, B=%d\n",dir,pbcshift[dir],B);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::pbcglide(double dx, double dy)
{
    double Ainv11,Ainv12,Ainv13,Ainv21,Ainv22,Ainv23,Ainv31,Ainv32,Ainv33;
    /* Ainv = inv (A) */
    Ainv11=A22*A33-A23*A32;
    Ainv22=A33*A11-A31*A13;
    Ainv33=A11*A22-A12*A21;
    Ainv12=A23*A31-A21*A33;
    Ainv23=A31*A12-A32*A11;
    Ainv31=A12*A23-A13*A22;
    Ainv13=A21*A32-A31*A22;
    Ainv21=A32*A13-A12*A33;
    Ainv32=A13*A21-A23*A11;
        
    pbcshift[0]+=Ainv11*dx+Ainv12*dy;
    pbcshift[1]+=Ainv21*dx+Ainv22*dy;
    pbcshift[2]+=Ainv31*dx+Ainv32*dy;
}

void YWindow::projectTo(XEvent ev)
{
    D*=exp(-(ev.xmotion.x-Xp+ev.xmotion.y-Yp)*0.001);
    if(D<0.2) D=0.2;
    if(D>1e2) D=1e2;
//    DUMP("YWindow: D="<<D);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::rotateTo(XEvent ev)
{
    
    double a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3;
    double r,rp,xx,yy,xp,yp;
#define _toler 1.0e-6
    xx=((double)ev.xmotion.x-width/2-X0)/R;
    yy=-((double)ev.xmotion.y-height/2-Y0)/R;
    xp=((double)Xp-width/2-X0)/R;
    yp=-((double)Yp-height/2-Y0)/R;
    
    if((xp-xx)*(xp-xx)+(yp-yy)*(yp-yy)<_toler) return;
//    DUMP("YWindow: Rotate from ("<<xp<<","<<yp<<") to ("<<xx<<","<<yy<<")");
    
    rp=sqrt(xp*xp+yp*yp);
    r=sqrt(xx*xx+yy*yy);
    if(r>=1 || rp>=1)
    {
        d1=xp/rp;
        e1=yp/rp;
        d2=xx/r;
        e2=yy/r;
        
        d3=cos((r-rp)*M_PI);
        e3=sin((r-rp)*M_PI);
        
        B11=d3*d1*d2+e1*e2;
        B12=d3*e1*d2-d1*e2;
        B13=e3*d2;
        B21=d3*d1*e2-e1*d2;
        B22=d3*e1*e2+d1*d2;
        B23=e3*e2;
        B31=-e3*d1;
        B32=-e3*e1;
        B33=d3;
    }
    else
    {
        a1=xp;a2=yp;a3=sqrt(1-a1*a1-a2*a2);
        b1=xx;b2=yy;b3=sqrt(1-b1*b1-b2*b2);
        c1=a2*b3-a3*b2;
        c2=a3*b1-a1*b3;
        c3=a1*b2-a2*b1;
        r=sqrt(c1*c1+c2*c2+c3*c3);
        c1/=r;c2/=r;c3/=r;
        
        d1=a2*c3-a3*c2;
        d2=a3*c1-a1*c3;
        d3=a1*c2-a2*c1;
        e1=b2*c3-b3*c2;
        e2=b3*c1-b1*c3;
        e3=b1*c2-b2*c1;
        
        B11=b1*a1+e1*d1+c1*c1;
        B12=b1*a2+e1*d2+c1*c2;
        B13=b1*a3+e1*d3+c1*c3;
        B21=b2*a1+e2*d1+c2*c1;
        B22=b2*a2+e2*d2+c2*c2;
        B23=b2*a3+e2*d3+c2*c3;
        B31=b3*a1+e3*d1+c3*c1;
        B32=b3*a2+e3*d2+c3*c2;
        B33=b3*a3+e3*d3+c3*c3;

    }

//    xp=xx;
//    yp=yy;

    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
    
    applyRotate();
}

void YWindow::Evolve()
{
    if(!alive) return;
    
    if(enableRot) applyRotate();
    if(dirty)
    {
        dirty=false;
        update();
    }
}

//============ Test suite

#ifdef _TEST

int main(int argc, char *argv[])
{
    int i, j, k;
    int ni=6,nj=6,nk=6;
    double x,y,z,x1,y1,z1,r,br,dx,dy,dz,dr;
    unsigned c, attr;
    char s[100];
    YWindow *win;
    
    win=new YWindow(400,400,"Test Window Display",true);
    
#define yw (*win)
    
    yw.Lock();
    yw.Clear();
    ni=nj=nk=6;

    yw.bgcolor=yw.AllocShortRGBColor(100,100,100);
    for(i=0;i<ni;i++)for(j=0;j<nj;j++)for(k=0;k<nk;k++)
    {
        sprintf(s,"Ball:%d,%d,%d",i,j,k);
        attr=1;
        x=-1+1./ni+2.*i/ni;y=-1+1./nj+2.*j/nj;z=-1+1./nk+2.*k/nk;r=.5/ni;br=r/4;
        c=yw.AllocShortRGBColor(i*0x33, j*0x33, k*0x33);
        yw.DrawPoint(x,y,z,r,c,s,attr);
        if(i>0)
        {
            x1=x-2./ni;y1=y;z1=z;
            dx=x1-x;dy=y1-y;dz=z1-z;dr=sqrt(dx*dx+dy*dy+dz*dz);dx/=dr;dy/=dr;dz/=dr;
            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,1);
//            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,0);
        }
        if(j>0)
        {
            x1=x;y1=y-2./nj;z1=z;
            dx=x1-x;dy=y1-y;dz=z1-z;dr=sqrt(dx*dx+dy*dy+dz*dz);dx/=dr;dy/=dr;dz/=dr;
            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,1);
//            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,0);
        }
        if(k>0)
        {
            x1=x;y1=y;z1=z-2./nk;
            dx=x1-x;dy=y1-y;dz=z1-z;dr=sqrt(dx*dx+dy*dy+dz*dz);dx/=dr;dy/=dr;dz/=dr;
            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,1);
//            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,0);
        }
    }
    
    yw.Unlock();
    yw.Refresh();
#ifdef NO_THREAD
    yw.Routine();
#else
    yw.Run();
    i=99;
#ifndef NO_GENERAL
    _IO << "I can live for <<"<<i<<" seconds.\n" ;
#else
    printf("I can live for %d seconds.\n",i);
#endif
    for(j=0;j<i;j++)
    {
        sleep(1);
#ifndef NO_GENERAL
        _IO << i-j << ' ';
#else
        printf("%d \n",i-j);
#endif
    }
#ifndef NO_GENERAL
    _IO << "Bye.\n";
#else
    printf("Bye.\n");
#endif
#endif
    return 0;
}

#endif //_TEST
