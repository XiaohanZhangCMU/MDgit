function writeMDPPeam(potfilename,griddata,elemdata,rhor,phirtimesr,Frho)
% write EAM potential for single element into MD++ format
% the derivative of the data needs to be computed first
% test:
% [griddata, elemdata, Frhodata, rhordata, phirdata] = ...
%    readLAMMPSeamfs('../../../../LAMMPS.svn/potentials/NiAlH_jea.eam.fs');
% writeMDPPeam('eamdata.Ni.LAMMPS.jea.eam.fs',griddata, ...
%    elemdata{1},rhordata{1}{1},phirdata{1}{1},Frhodata{1});
% this program will convert phi(r)*r into phi(r)

Nr   = griddata.Nr;
Nrho = griddata.Nrho;
if (Nr~=Nrho)
    disp(sprintf('Error: Nr = %d  Nrho = %d  not equal to each other!',Nr,Nrho));
    return;
end

fid = fopen(potfilename,'w');

fprintf(fid,'EAM potential file for MD++ written by writeMDPPeam.m\n');
fprintf(fid,' %4d %20.6E %20.6E\n',elemdata.Z,elemdata.atommass,elemdata.latticeconst);
fprintf(fid,'%26.16E %26.16E %26.16E %26.16E\n',griddata.dr,griddata.drho,griddata.rcut,griddata.rmin);

% compute derivatives
r = [0:Nr-1]*griddata.dr + griddata.rmin;  rho = [0:Nrho-1]*griddata.drho;
cs = spline(r,rhor);  
csp=cs; csp.order=3; csp.coefs=[cs.coefs(:,1)*3,cs.coefs(:,2)*2,cs.coefs(:,3)];
rhopr = ppval(csp,r);

if (r(1) == 0)
    phir = phirtimesr(2:end)./r(2:end); phir = phir([1,1:end]);
    phir(1) = 2*phir(2)-phir(3);
else
    phir = phirtimesr./r;
end
cs = spline(r,phir);  
csp=cs; csp.order=3; csp.coefs=[cs.coefs(:,1)*3,cs.coefs(:,2)*2,cs.coefs(:,3)];
phipr = ppval(csp,r);

cs = spline(rho,Frho);  
csp=cs; csp.order=3; csp.coefs=[cs.coefs(:,1)*3,cs.coefs(:,2)*2,cs.coefs(:,3)];
Fprho = ppval(csp,rho);

% write rho(r) and rhop(r)
for i=1:Nr,
  fprintf(fid,'%26.16E %26.16E %26.16E %26.16E\n',rhor(i),rhopr(i),rhor(i),rhopr(i));
end

% write phi(r) and phip(r)
for i=1:Nr,
  fprintf(fid,'%26.16E %26.16E %26.16E %26.16E\n',phir(i),phipr(i),phir(i),phipr(i));
end

% [do not write phix(r) and phixp(r)]

% write F(rho) and Fp(rho)
for i=1:Nr,
  fprintf(fid,'%26.16E %26.16E %26.16E %26.16E\n',Frho(i),Fprho(i),Frho(i),Fprho(i));
end

fclose(fid);

% plot data
figure(100);
p1=plot(r,rhor,r,rhopr); set(p1(1),'LineWidth',3); 
xlabel('r');  ylabel('\rho(r)');
figure(101);
p2=plot(r,phir,r,phipr); set(p2(1),'LineWidth',3);
xlabel('r');  ylabel('\phi(r)');
figure(102);
p3=plot(rho,Frho,rho,Fprho); set(p3(1),'LineWidth',3);
xlabel('\rho');  ylabel('F(\rho)');
