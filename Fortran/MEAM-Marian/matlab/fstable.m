%compute the table for FS potential

%Fe-Fe (Ackland et al. PMA 75, 713, 1997)
a0 = 2.8665;

a = [
    -36.559853
     62.416005
    -13.155649
     -2.721376
      8.761986
    100.0
]/a0^3;

A = [ 72.868366
    -100.944815 ]/a0^3;

r = [ 1.18
      1.15
      1.08
      0.99
      0.93
      0.866025 ]*a0;
  
R = [ 1.3
      1.2 ]*a0;
  
x1 = 0.90;
x2 = 1.90;
B0 = 7.14705133;
B1 = 0.69010282;
B2 =-4.16604662;
B3 = 1.06871772;
Z = 26;

ab = 0.529177; %Bohr radius (A)
as=0.88534*ab/(Z^(2/3)+Z^(2/3))^(1/2);


rmin=0.100000001490116 ;
%rmax=3.38247000000000 ;
rmax = 3.73042000000000 ;
dr = 6.564939997019770e-4 ;
%dr = dr/10;

x=[rmin:dr:rmax]';

c1=7.615551295656224e+000
c2=9.849458675137477e-001;

V = zeros(size(x));
dV = V;
for i=1:length(a),
    V = V+a(i)*(r(i)-x).^3.*(x<r(i)).*(x>=x2);
    dV=dV+(-3)*a(i)*(r(i)-x).^2.*(x<r(i)).*(x>=x2);
end
V = V+c1*exp(B0+B1*x+B2*x.*2+B3*x.*3).*(x>x1).*(x<x2);
dV=dV+c1*exp(B0+B1*x+B2*x.*2+B3*x.*3).*(B1+2*B2.*x+3*B3.*x.^2).*(x>x1).*(x<x2);

%universal potential
y=x/as;
eps0=8.854187816e-12;
k=1/4/pi/eps0*1.602e-19*1e10;

phi=0.1818*exp(-3.2*y)+0.5099*exp(-0.9423*y)+0.2802*exp(-0.4029*y)...
    +0.02817*exp(-0.2016*y);
dphi=(0.1818*(-3.2)*exp(-3.2*y)...
    +0.5099*(-0.9423)*exp(-0.9423*y)...
    +0.2802*(-0.4029)*exp(-0.4029*y)...
    +0.02817*(-0.2016)*exp(-0.2016*y))/as;
V = V+c2* k*Z*Z./x.*phi.*(x<=x1);
dV=dV+ (c2* k*Z*Z*(-1)./(x.^2).*phi + c2* k*Z*Z./x.*dphi).*(x<=x1);

Phi=A(1)*(x<R(1)).*(R(1)-x).^3 + A(2)*(x<R(2)).*(R(2)-x).^3;
Phi=Phi.*(Phi>0);
dPhi=A(1)*(x<R(1)).*(-3).*(R(1)-x).^2 + A(2)*(x<R(2)).*(-3).*(R(2)-x).^2;
dPhi=dPhi.*(Phi>0);

subplot(2,1,1);
plot(x,V,'b',x,dV/10,'r');
%plot(V);
subplot(2,1,2);
plot(x,Phi,'b',x,dPhi,'r');



y=x1/as;
phi=0.1818*exp(-3.2*y)+0.5099*exp(-0.9423*y)+0.2802*exp(-0.4029*y)...
    +0.02817*exp(-0.2016*y);
V11 = c2*k*Z*Z./x1.*phi;
V12 = c1*exp(B0+B1*x1+B2*x1.*2+B3*x1.*3);

V21 = c1*exp(B0+B1*x2+B2*x2.*2+B3*x2.*3);
V22 = a(1)*(r(1)-x2).^3+a(2)*(r(2)-x2).^3+a(3)*(r(3)-x2).^3 ...
     +a(4)*(r(4)-x2).^3+a(5)*(r(5)-x2).^3+a(6)*(r(6)-x2).^3 ; 
[V11,V12;V21,V22]

eV = 1.60219e-19;
angstrom = 1.0e-10;

data=[x,Phi*(eV*eV),V*(eV),dPhi*(eV*eV/angstrom),dV*(eV/angstrom)];

fp=fopen('fs2.out','w');
for i=1:length(x),
    fprintf(fp,'%20.16e %20.16e %20.16e %20.16e %20.16e\n',...
        data(i,1),data(i,2),data(i,3),data(i,4),data(i,5));
end
fclose(fp);