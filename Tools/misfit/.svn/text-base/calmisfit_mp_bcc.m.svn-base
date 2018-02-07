%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: calmisfit_mp_bcc.m
%%
%% Last Modified : Sat Oct  8 01:03:56 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate misfit energy of a perfect crystal on a set of planes
% the parameter set is designed for BCC crystals
%
% to use this matlab script, MD++ has to be compiled with static library
% e.g.:
%       make alglue SPEC=-static build=R

%size of the supercells

crystalstructure = 'body-centered-cubic';
latticeconst = 3.1472 % in Angstrom (for Mo)
rlist = 4.6; % neighborlist radius
executable = 'fs_gpp';
scriptfile = 'mo_misfit.script';
rundir = 'mo_misfit';

N=[0:10,25];
Sdata=zeros(length(N),2);  %ideal shear strength
dist =zeros(size(N));      %interplanar distance

for iter=1:length(N),

  m=1;
  n=N(iter);

  vx = [1 1 1];
  vz = [1 -2 1]*m + [-1 -1 2]*n;
  vy = cross(vx,vz)/3;
  
  dist(iter) = latticeconst/norm(vy);
  
  nx = ceil(2*rlist/(norm(vx)*latticeconst));
  ny = ceil(2*rlist/(norm(vy)*latticeconst));
  nz = ceil(2*rlist/(norm(vz)*latticeconst));
  
  disp(sprintf('m = %d  n = %d',m,n));

  latticesize = [ vx   nx+2
                  vy   ny 
                  vz   nz ]
  
  normaldir = 2; %surface normal (along y direction)
  
  dxyz = [ 0 0.005 0.5
           0 0.05 0
           0 0.05 0 ];
  
  enable_plot = 1;
  
  dx = dxyz(1,2);
  dy = dxyz(2,2);
  dz = dxyz(3,2);
  
  Lx = latticeconst*norm(vx);
  Ly = latticeconst*norm(vz);
  Lz = latticeconst*norm(vy);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % call MD++ to compute Edata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  callmdpp
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % visualization
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %plotEmisfit
  %drawnow
  
  sigmax=0;
  sigmin=0;
  if (sum([Mx,My,Mz]>1)==1)
    disp('E is a 1-dimensional array');  
    E = Edata(:,4);
    %find maximum slope dE/dx
    sigmax = max(E(2:end)-E(1:end-1))/(dx*Lx) * 160.2;
    sigmin = min(E(2:end)-E(1:end-1))/(dx*Lx) * 160.2;
    disp(sprintf('sigmax = %e (GPa)',sigmax));
  end  
  Sdata(iter,:) = [sigmax,sigmin];
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save mo_mp_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotEmisfit_mp

