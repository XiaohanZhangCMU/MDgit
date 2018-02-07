%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: calmisfit_mp_bcc.m
%%
%% Last Modified : Tue Oct 11 15:00:34 2005
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
%N=[0];

Sdata=zeros(length(N),4);  %ideal shear strength
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
  
  dxyz0 = [ 0 0.005 0.5
        -0.05 0.01 0.1
           0 0.05 0 ];
  
  enable_plot = 1;
  
  Lx = latticeconst*norm(vx);
  Ly = latticeconst*norm(vy);
  Lz = latticeconst*norm(vz);

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
%save mo_mp_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotEmisfit_mp

save mo_rlx_mp_data_100
