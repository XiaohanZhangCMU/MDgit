/*
  heaziz.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Aug 22 2011 
  Notes: with Maurice de Koning

  FUNCTION  :  He MD simulation using Aziz Potential
               Parameters       Phys. Rev. Lett. 74, 1586 (1995) 
               Functional Form  Mol. Phys. 77, 321 (1992)

*/

#ifndef _HEAZIZ_H
#define _HEAZIZ_H

#include "mdparallel.h"


class HEAZIZFrame: public MDPARALLELFrame
{

    double _ESCALE, _RM, _A, _ALPHA, _BETA, _C6, _C8, _C10, _D;
    double AZIZ_RC;

public:

#if 1 /* 95 version */
    HEAZIZFrame():_ESCALE(9.306723E-04), /* 10.956 K * KB (eV/K) */
                  _RM(2.9683), _A(1.86924404E+05), _ALPHA(10.5717543), _BETA(-2.07758779),
                  _C6(1.35186623), _C8(0.41495143), _C10(0.17151143), _D(1.438),
                  AZIZ_RC(8.0) {};
#endif

#if 0 /* older version http://phycomp.technion.ac.il/~phsorkin/Thesis_html/node10.html */
    HEAZIZFrame():_ESCALE(9.30672288E-04), /* 10.8 K * KB (eV/K) */
                  _RM(2.9673), _A(0.54485046e6), _ALPHA(13.35338), _BETA(0),
                  _C6(1.3732412), _C8(0.4253785), _C10(0.1781), _D(1.241314),
                  AZIZ_RC(8.0) {};
#endif

    virtual void initvars();
    virtual void potential();
    void he_aziz();
};        



#endif // _HEAZIZ_H

