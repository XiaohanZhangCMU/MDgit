NOTES:
This is the wiki page for some preliminary test cases of fiber+fem and NEB+fem. 
(X)FEMFrame implements a new method that solves finite element (discretized) problem within a MD framework. 
It currently requires all elements in the reference configuration have the same shape and size. 
It needs the following data:
1) mesh connectivity, 2) FE node coordinates 3) FE coefficients. dF/dx. 4) boundary information
The matlab code that generates all these informations is thin_film_refined_relaxation.m Change the mesh size, xy/xz plane, box size, and refinement around line 9-13


COMPILE:

For fiber+fem:

make cgmd SYS=mc2 build=R
cgmd class is a child class of XFEMFrame and LJBONDFrame. XFEMFrame is a child of FEMFrame
LJBONDFrame uses ljbond3.cpp and ljbond3.h. 

For NEB+fem:
compile: make xfem SYS=mc2_mpich build=R

----------------------------------------------------------------------

TEST CASE:

fe-mesh-data:
Five different meshes and their corresponding data that will be needed for code, and all the following test cases
M1) membrane2d-5x1element.cn (.bdy, .ele,  CPE4-Q4.dat)
   5 elements, xy plane, [-4,1]x [0,1], boxsize = 10
M2) memebrane2d-8x5element.cn (.bdy, .ele, CPE4-Q4.dat)
   8x5 elements, xy plane. [-4,4] x [-2.5,2.5], boxsize = 20
M3) memebrane2d-8x5element-xz.cn (.bdy, .ele, CPE4-Q4.dat)
   same as 2), but on xz plane
M4) memebrane2d-8x5element-xz.cn (.bdy, .ele, CPE4-Q4.dat)
   8x5 elements, xz plane. [-3500,3500] x [-3500,3500], boxsize = 8000
M5) membrane2d-16x3element.cn (.bdy, .ele, CPE4-Q4.dat)
   16x3 elements, xy plane, [-4,4] x [0,1], boxsize = 10

nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn

NEB-FEM TEST CASES: [ status , problem id, mesh-id ]

1) xfem neb-fem.tcl 0 1  1 --- 5
This tests M1 --- M5 under extension by moving boundary nodes.

2) xfem neb-fem.tcl 0 2  1 --- 5
This tests M1 --- M5 under extension by increase box size.

3) xfem neb-fem.tcl 0 3  1 --- 5
This tests M1 --- M5 under compression by increase box size.

4) xfem neb-fem.tcl 0 4  7
This compresses M5 and creates a buckled shape
change perturbation force and creates the other states of buckled shape

fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FIBER-FEM TEST CASES:
1) cgmd fiber-fem.tcl 0 1  4
A ) The result gives the final relaxed energy of fiber+membrane+interaction

B) comment out: 

initialize_substrate and 
MD++ dockbeads. 

gives final relaxed energy of fiber only

C) comment out only:

MD++ dockbeads

D) comment out 
initialize_config

The verfication is energy of A = B+C. D =0

2) cgmd fiber-fem.tcl 0 2 4 

extension is applied by box increasing

A) the result gives the final relaxed energy under simple extension. 

B) the result gives the final relaxed energy with only fiber

C) the result gives the final relaxed energy with only membrane

check : C = NEB-TEST:: 2)-4. 

check: A >= B+C

3) A full simulation during cycle loading with interactions between fibers and membrane.
