%plot FS2.Fe

load FS2Fe.out

x=FS2Fe(:,1);
rho=FS2Fe(:,2);
V=FS2Fe(:,3);
drho=FS2Fe(:,4);
dV=FS2Fe(:,5);


eV = 1.60219e-19;
angstrom = 1.0e-10;

rho=rho/(eV*eV);
V=V/(eV);
drho=drho*angstrom/(eV*eV);
dV=dV*angstrom/eV;
 
subplot(2,1,1);
plot(x,V,'b',x,dV/10,'r');
%plot(V);
subplot(2,1,2);
plot(x,rho,'b',x,drho,'r');
