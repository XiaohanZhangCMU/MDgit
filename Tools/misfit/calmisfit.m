%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: calmisfit.m
%%
%% Last Modified : Mon Oct 17 15:03:56 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate misfit energy of a perfect crystal on one plane
%
% to use this matlab script, MD++ has to be compiled with static library
% e.g.:
%       make alglue SPEC=-static build=R

%size of the supercells
Nx = 5;
Ny = 3;
Nz = 3;

%crystalstructure = 'face-centered-cubic';
%latticeconst = 4.05; % in Angstrom (for Al)
%executable = 'alglue_gpp';
%scriptfile = 'al_misfit.script';
%rundir = 'al_misfit';
%latticesize = [ 1 -1  0 6   %(n=0)
%               -1 -1  1 3
%                1  1  2 3 ];
%latticesize = [ 1 -1  0 6  %(n=infinity)
%               -1 -1 -1 3
%               -1 -1  2 3 ];

crystalstructure = 'face-centered-cubic';
latticeconst = 5.260 % in Angstrom (for Ar)
executable = 'lj_gpp';
scriptfile = 'ar_misfit.script';
rundir = 'ar_misfit';
latticesize = [  1  -1  0  6
                 1   1  1  3
                 1   1 -2  3 ];

%crystalstructure = 'body-centered-cubic';
%latticeconst = 3.1472 % in Angstrom (for Mo)
%executable = 'fs_gpp';
%scriptfile = 'mo_misfit.script';
%rundir = 'mo_misfit';
%latticesize = [ .5 .5 .5   4 
%                -1  0  1   4 
%                 1 -2  1   4 ]

normaldir = 2; %surface normal (along y direction)

dxyz = [ 0 0.005 0.5
         0 0.05 0
         0 0.05 0 ];

enable_plot = 1;

dx = dxyz(1,2);
dy = dxyz(2,2);
dz = dxyz(3,2);

vx = latticesize(1,1:3);
vy = latticesize(2,1:3);
vz = latticesize(3,1:3);

Lx = latticeconst*norm(vx);
Ly = latticeconst*norm(vz);
Lz = latticeconst*norm(vy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call MD++ to compute Edata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
callmdpp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save mo_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotEmisfit
