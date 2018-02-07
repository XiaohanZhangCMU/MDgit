%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: calmisfit_mp_fcc.m
%%
%% Last Modified : Wed Oct 12 18:21:31 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate misfit energy of a perfect crystal on a set of planes
% the parameter set is designed for FCC crystals
%
% to use this matlab script, MD++ has to be compiled with static library
% e.g.:
%       make alglue SPEC=-static build=R

%size of the supercells

crystalstructure = 'face-centered-cubic';
latticeconst = 4.05; % in Angstrom (for Al)
rlist = 6.2; %neighborlist radius
executable = 'alglue_gpp';
scriptfile = 'al_misfit.script';
rundir = 'al_misfit';

%crystalstructure = 'face-centered-cubic';
%latticeconst = 5.260 % in Angstrom (for Ar)
%rlist = 9;
%executable = 'lj_gpp';
%scriptfile = 'ar_misfit.script';
%rundir = 'ar_misfit';

%N=[-1, 0, 1:10, 25];
%N=[0];
%M = 1;
%N = ones(size(N));

N = [3 4 5 5 5 -1 -1 -1];
M = [2 3 2 3 4  2  3  4];

Sdata=zeros(length(N),4);  %ideal shear strength
dist =zeros(size(N));      %interplanar distance

for iter=1:length(N),

  m=M(iter);
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
  
  dxyz0 = [ 0 0.005 0.5
        -0.05 0.01 0.15
           0 0.05 0 ];
  
  enable_plot = 1;
  
  Lx = latticeconst*norm(vx);
  Ly = latticeconst*norm(vz);
  Lz = latticeconst*norm(vy);

  dxyz=dxyz0;
  dxyz(2,:)=dxyz0(2,:)*7/Ly;    
  
  dx = dxyz(1,2);
  dy = dxyz(2,2);
  dz = dxyz(3,2);    
  
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
  sigmax0=0;
  sigmin0=0;
  if (sum([Mx,My,Mz]>1)==1)
    disp('E is a 1-dimensional array');  
    E = Edata(:,4);
    %find maximum slope dE/dx
    sigmax = max(E(2:end)-E(1:end-1))/(dx*Lx) * 160.2;
    sigmin = min(E(2:end)-E(1:end-1))/(dx*Lx) * 160.2;
    disp(sprintf('sigmax = %e (GPa)',sigmax));
  elseif (sum([Mx,My,Mz]>1)==2)
    disp('E is a 2-dimensional array');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert data to relaxed energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = reshape(Edata(:,4),Mx,My);
    Erlx = zeros(Mx,1);
    Zrlx = zeros(Mx,1);
    
    normalaxis = dxyz(normaldir,1):dxyz(normaldir,2):(dxyz(normaldir,3)-0*dxyz(normaldir,2));
    normalgrid = dxyz(normaldir,1):dxyz(normaldir,2)/1000:(dxyz(normaldir,3)-0*dxyz(normaldir,2));
    
    a = zeros(1,My);
    for i = 1:Mx,
      a(:) = E(i,:);
      a_interp = interp1(normalaxis,a,normalgrid,'spline');
      [Erlx(i), ind] = min(a_interp);
      Zrlx(i) = normalgrid(ind);
      if ((Zrlx(i)<=normalaxis(1))|(Zrlx(i)>=normalaxis(end)))
        disp(sprintf('Zrlx(%d) = %e exceeds bound (%e, %e) !\n',i,Zrlx(i),normalaxis(1),normalaxis(end)));
      end
    end
    %find maximum slope dErlx/dx
    sigmax = max(Erlx(2:end)-Erlx(1:end-1))/(dx*Lx) * 160.2;
    sigmin = min(Erlx(2:end)-Erlx(1:end-1))/(dx*Lx) * 160.2;

    %retriev unrelaxed energy
    Eunrlx=E(:,find(abs(normalaxis)<1e-8));
    %find maximum slope dEunrlx/dx
    sigmax0 = max(Eunrlx(2:end)-Eunrlx(1:end-1))/(dx*Lx) * 160.2;
    sigmin0 = min(Eunrlx(2:end)-Eunrlx(1:end-1))/(dx*Lx) * 160.2;

    disp(sprintf('sigmax = %e   sigmax0 = %e (GPa)',sigmax,sigmax0));

  end

  Sdata(iter,:) = [sigmax,sigmin,sigmax0,sigmin0];
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save al_mp_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotEmisfit_mp

save al_rlx_mp_data_100_nm2
