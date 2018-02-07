% plot fs potentials in eam format

%Mo         d        A        beta  c    c0          c1         c2        B      alpha  b0
fs_data = [ 4.114825 1.887117 0     3.25 43.4475218 -31.9332978 6.0804249 1223.0 3.90   2.7255];

%Ta         d        A        beta  c    c0          c1         c2        B      alpha  b0
%fs_data = [ 4.076980 2.591061 0     4.20 1.2157373   0.027147  -0.121735  91.0   1.05   2.8629];

%Nb         d        A        beta  c    c0          c1         c2        B      alpha  b0
%fs_data = [ 3.915354 3.013789 0     4.20 -1.5640104  2.0055779 -0.4663764 48.0   0.8   2.8585];

d = fs_data(1);
A = fs_data(2);
beta = fs_data(3);
c = fs_data(4);
c0= fs_data(5);
c1= fs_data(6);
c2= fs_data(7);
B = fs_data(8);
alpha= fs_data(9);
b0= fs_data(10);

fac = 1.0;
%fac = 1.8871; %conversion factor between fs and LAMMPS (Mo.AT1.fs)

%LAMMPS potential file to compare with
if (~exist('elemdata'))
  [griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs('./Mo.stanford.fs');
  %[griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs('./Nb.stanford.fs');
  %[griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs('../../../LAMMPS.svn/potentials/Mo.mod.fs');
  %[griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs('../../../LAMMPS.svn/potentials/Mo.AT1.fs');
  %[griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs('../../../LAMMPS.svn/potentials/Fe.eam.fs');
end

rho = [0:griddata.Nrho-1]*griddata.drho;
r   = [0:griddata.Nr-1  ]*griddata.dr;

%analytic expression of FS potential (Ackland)
Frho = -A*sqrt(rho);
rhor = (r<d).*((r-d).^2 + beta*(r-d).^3/d);
phir = (r<c).*(r-c).^2.*(c0+c1.*r+c2.*r.^2) + (r<b0).*(B*(b0-r).^3.*exp(-alpha*r));

figure(1);
subplot(1,3,1);
plot(rho,Frho/fac,rho,Frhodata{1},'--'); xlabel('\rho'); ylabel('F(\rho)');
subplot(1,3,2);
plot(r,rhor*fac^2,r,rhordata{1}{1},'--'); xlabel('r'); ylabel('\rho(r)');
subplot(1,3,3);
plot(r,phir.*r,r,phirdata{1}{1},'--'); xlabel('r'); ylabel('\phi(r) \cdot r');

