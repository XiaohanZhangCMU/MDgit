%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: plot_eam.m
%%
%% Chris Weinberger, cweinber@stanford.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fid = fopen('../../potentials/EAMDATA/auu4.eam');
fid = fopen('../../../LAMMPS.svn/potentials/Au_u3.eam');

comment = fgets(fid);

ielem = fscanf(fid,'%d',1);
amass = fscanf(fid,'%f',1);
blatt = fscanf(fid,'%f',1);
lat   = fscanf(fid,'%s',1);

nrho = fscanf(fid,'%d',1);
drho = fscanf(fid,'%f',1);
nr   = fscanf(fid,'%d',1);
dr   = fscanf(fid,'%f',1);
rcut = fscanf(fid,'%f',1);

ielem,amass,blatt,lat,nrho,drho,nr,dr,rcut

Frho = zeros(nrho,1);
Zr   = zeros(nr,1);
rhor = zeros(nr,1);

Frho = fscanf(fid,'%f',nrho);
Zr   = fscanf(fid,'%f',nr);
rhor = fscanf(fid,'%f',nr);

fclose(fid);

figure(2); 
subplot(1,3,1); plot([0:nrho-1]*drho,Frho); xlabel('\rho'); ylabel('F(\rho)');
subplot(1,3,2); plot([0:nr-1]*dr,Zr);       xlabel('r');    ylabel('Z(r)');
subplot(1,3,3); plot([0:nr-1]*dr,rhor);     xlabel('r');    ylabel('\rho(r)');
