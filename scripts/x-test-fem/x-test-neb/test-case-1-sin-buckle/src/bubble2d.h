/*
  bubble2d.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Thu Apr 10 2014

  FUNCTION  :  bubble 2d simulation
*/

#ifndef _BUBBLE2D_H
#define _BUBBLE2D_H

#include "mdparallel.h"

class BUBBLE2DFrame : public MDPARALLELFrame /* Bubbles with quadratic potential */
{
    /* Bubble 2d potential parameters */
    double _radius0;                 /* original radius of bubble */
    double _sigma0;                  /* surface tension */
    double _kwall, _slope, _channel; /* external potential */
    
public:
    BUBBLE2DFrame(): _radius0(1),_sigma0(1000),_kwall(1000),_slope(0.1),_channel(10){};
    void bubble_2d();
    virtual void potential();
    virtual void initvars();

    virtual void initparser();
    virtual int exec(const char *name);
    virtual void calcprop();
    
};

#endif // _BUBBLE2D_H

