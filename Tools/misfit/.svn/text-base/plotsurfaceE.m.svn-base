%function plotmisfit(s)
% Calculate and plot the misfit energy.
% Usage: plotmisfit('datafile')
%    eg. plotmisfit('~/Codes/MD++/runs/mo_misfit/Emisfit.out')
% Variables
%  E0 : potential energy of perfect crystal per unit area (eV/A^2)
%  Emisfit : misfit energy (eV/A^2)
% Keonwook Kang, kwkang@stanford.edu
% Apr 28 2005

%clear all

%s='~/Codes/MD++/runs/mo_misfit/Emisfit.out';
filename='~/Codes/MD++.old/runs/si_misfit/SurfE_SW_110C.out';
%filename='~/Codes/MD++/runs/si_misfit/Emisfit.out';
%filename='~/Codes/MD++/runs/ge_misfit/Emisfit.out';
%filename='~/Codes/MD++/runs/ge_misfit/SurfE_SW_100.out';
%filename='~/Codes/MD++/runs/mo_misfit/Emisfit.out';

% calculate the lattice vector magnitude along the normal (or Z) direction
latticeconstant = 5.4309; % in Angstrom for Si
%latticeconstant = 5.6536; % in Angstrom for Ge
%latticeconstant = 3.1472; % in Angstrom for Mo 
onelatticevector = sqrt(2)*latticeconstant; % sqrt(abs([110]))*a 
%onelatticevector = sqrt(3)*latticeconstant; % sqrt(abs([111]))*a
%onelatticevector = latticeconstant; % sqrt(abs([100]))*a 
% for Mo
%onelatticevector = sqrt(3)*latticeconstant; % sqrt(abs([111]))*a

% read data
line1 = []; line2 = []; cndata = [];
fid = fopen(filename,'r');
lines = fgetl(fid);
lines = fgetl(fid);
[line1] = fscanf(fid, '%d %d %d %lf %lf %lf\n', [6 1]);
stepx = line1(1,1); stepy = line1(2,1); stepz = line1(3,1);
z0 = line1(4,1); z1 = line1(5,1); E0 = line1(6,1);
garbage = fscanf(fid, '%s',[1 1]); ind1 = fscanf(fid,'%d\t',[1 1]);
garbage = fscanf(fid, '%s',[1 1]); ind2 = fscanf(fid,'%d\t',[1 1]);
garbage = fscanf(fid, '%s',[1 1]); ind3 = fscanf(fid,'%d\t',[1 1]);
garbage = fscanf(fid, '%s\n',[1 1]);
line3 = fgetl(fid);
N = (stepx+1)*(stepy+1)*(stepz+1);
[cndata] = fscanf(fid,'%d %d %d %lf',[4 N]);
                     % x  y  z  PotE
if(0)                   
if i == 1
    if j == 2
         x = cndata(1,:); y = cndata(2,:); z = cndata(3,:);
    elseif j == 3
         x = cndata(1,:); z = cndata(2,:); y = cndata(3,:);
    end
elseif i == 2
    if j==3
         y = cndata(1,:); z = cndata(2,:); x = cndata(3,:);
    elseif j==1
         y = cndata(1,:); x = cndata(2,:); z = cndata(3,:);
    end
elseif i == 3
    if j==1
         z = cndata(1,:); x = cndata(2,:); y = cndata(3,:);
    elseif j==2
         z = cndata(1,:); y = cndata(2,:); x = cndata(3,:);
    end
end
end

Emisfit = cndata(4,:)-E0;

fclose(fid);

%figure(1);
%size(Emisfit)

 
X = [0:stepx]/stepx;
Y = [0:stepy]/stepy;
Z = [0:stepz]/stepz;
%
%I = 0;
%EmisfitI = reshape(Emisfit((stepx+1)*(stepy+1)*I+1:(stepx+1)*(stepy+1)*(I+1)), ...
%              [stepx+1 stepy+1]);
Emisfit_reshaped = reshape(Emisfit, [ stepx+1 stepy+1 stepz+1]);

Emisfitrelax = zeros(stepx+1, stepy+1);
Eminind = Emisfitrelax;

zaxis = (1-Z)*z0 + Z*z1;

%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amax = onelatticevector*z1; %
amin = onelatticevector*z0;
%amax = planeseparationdistance_eq*z1; %
%amin = planeseparationdistance_eq*z0;
a = linspace(amin,amax,length(Emisfit_reshaped(1,1,:)));
SurfaceE =zeros(1,length(Emisfit_reshaped(1,1,:)));
%a = linspace(amin,amax,length(Emisfit_reshaped(1,:,1)));
%SurfaceE =zeros(1,length(Emisfit_reshaped(1,:,1)));
for i=1:stepz+1
    SurfaceE(i) = Emisfit_reshaped(1,1,i);
    %SurfaceE(i) = Emisfit_reshaped(1,i,1);
    %zmin = fminbnd( 'Emsftz', z0, z1, [], zaxis, a);
end

step = a(2)-a(1);
ai = a(1):amax/100:amax;
stress = zeros(length(SurfaceE),1);
stress(1) = (SurfaceE(2)-SurfaceE(1))/step;
for i=2:length(SurfaceE)-1
     stress(i) = (SurfaceE(i+1)-SurfaceE(i-1))/(2*step);
end
stress(end) = (SurfaceE(end)-SurfaceE(end-1))/step;
stress_interp = interp1(a, stress, ai,'spline');
[tensile, I] = max(stress_interp);

disp(['Surface Energy = ', num2str(SurfaceE(end)), '(eV/Angs^2)'])
disp(['Tensile Strength = ', num2str(tensile), '(eV/Angs^3)'])
disp(['The distance where the max. stress is obatained = ', num2str(ai(I)),'(Angs)'])

if(0)
figure(1)
plot(a,SurfaceE,'.-')
set(gca,'fontsize',12)
xlabel('distance, a (Angstrom)','fontsize',14)
ylabel('Surface Energy (eV/A^2)','fontsize',14) % primitive lattice vector a1, a2, and a3
title('(110) Surface Energy of Ge bulk','fontsize',16)

figure(2)
plot(a,stress,'.',ai,stress_interp)
set(gca,'fontsize',12)
xlabel('distance, a (Angstrom)','fontsize',14)
ylabel('Stress (eV/A^3)','fontsize',14) % primitive lattice vector a1, a2, and a3
title('Tensile stress along [110] direction','fontsize',16)
end

% n=20;
% a=a/n;
% SurfaceE=SurfaceE/n;
 C=polyfit(a,SurfaceE,2);

figure(2);
plot(a,SurfaceE,'o', a,polyval(C,a),'-');

figure(1);
clf
hl1 = line(a,SurfaceE,'Color','b','LineStyle','-','marker','o');
ax1 = gca;
set(ax1,'YColor','b')
ax2 = axes('Position',get(ax1,'Position'),...
           'YAxisLocation','right',...
           'Color','none',...
           'YColor','r');
hl2 = line(ai,stress_interp,'Color','r','Parent',ax2);
hl3 = line(a,stress,'Color','r','linestyle','none','marker','o','Parent',ax2);

xlabel('distance, a (Angstrom)','fontsize',14)
set(get(ax1,'Ylabel'),'String','(110) Surface Energy (eV/A^2)','fontsize',14)
set(get(ax2,'Ylabel'),'String','[110] Tensile Strength (eV/A^3)','fontsize',14)
title('Ge Bulk Using Tersoff Potential','fontsize',16)
