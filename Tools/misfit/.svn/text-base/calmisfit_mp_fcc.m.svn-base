%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: calmisfit_mp_fcc.m
%%
%% Last Modified : Tue Oct 18 10:42:35 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate misfit energy of a perfect crystal on a set of planes
% the parameter set is designed for FCC crystals
%Sda
% to use this matlab script, MD++ has to be compiled with static library
% e.g.:
%       make alglue SPEC=-static build=R

%size of the supercells

%crystalstructure = 'face-centered-cubic';
%latticeconst = 4.05; % in Angstrom (for Al)
%rlist = 6.2; %neighborlist radius
%executable = 'alglue_gpp';
%scriptfile = 'al_misfit.script';
%rundir = 'al_misfit';

crystalstructure = 'face-centered-cubic';
latticeconst = 5.260 % in Angstrom (for Ar)
rlist = 9;
executable = 'lj_gpp';
scriptfile = 'ar_misfit.script';
rundir = 'ar_misfit';

%N=[-1, 0, 1:25];
%M=1;

N = [-1 0 1 2 3 4 5 6 7 8 9 10 25 3 4 5 5 5 -1 -1 -1 -1 -2 -2 -2 -3 -3 2 3];
M = [ 1 1 1 1 1 1 1 1 1 1 1  1  1 2 3 2 3 4  1  2  3  4  3  5  7  4  5 3 4];


Sdata=zeros(length(N),2);  %ideal shear strength
dist =zeros(size(N));      %interplanar distance
Esave = zeros(length(N),101);

for iter=1:length(N),
  if(length(M)==1)
    m=M;
  else
    m=M(iter);
  end
  n=N(iter);

  vx = [1 -1 0];
  vy = [-1 -1 1]*m + [-1 -1 -1]*n;
  vz = ([1 1 2]*m + [-1 -1 2]*n)/2;

  if(mod(m+n,2)==0)
    vy = vy/2;
    vz = vz/2;
  end
  
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
  Esave(iter,:) = E';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save al_mp_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotEmisfit_mp

