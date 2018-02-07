%function plotmisfit(s)
% Calculate and plot the misfit energy.
% Usage: plotmisfit('datafile')
%    eg. plotmisfit('~/Codes/MD++/runs/mo_misfit/Emisfit.out')
% Variables
%  E0 : potential energy of perfect crystal per unit area (eV/A^2)
%  Emisfit : misfit energy (eV/A^2)
% Keonwook Kang, kwkang@stanford.edu
% Apr 28 2005

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the output data file

%s='~/Codes/MD++/runs/mo_misfit/Emisfit.out';
%filename='~/Codes/MD++/runs/mo_misfit/Emisfit.out';
filename='~/Codes/MD++/runs/al_misfit/Gamma110_v2.out';
%filename='~/Codes/MD++/runs/si_misfit/Emisfit.out';
%filename='~/Codes/MD++/runs/ge_misfit/Emisfit.out';
fid = fopen(filename,'r');

line1 = []; line2 = []; cndata = [];
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

stepx0 = stepx; stepy0 = stepy; stepz0 = stepz;
if ind3 == 1 % if normal to misfit plane is x-direction
    tmp = stepz; stepz = stepx; stepx = stepy; stepy = tmp;
end
if ind3 == 2 % if normal to misfit plane is y-direction
    tmp = stepz; stepz = stepy; stepy = stepx; stepx = tmp;
end
                     
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
zi = z0:(z1-z0)/1000:z1; % added for my way of finding min.
a = zeros(length(Emisfit_reshaped(1,1,:)),1);
for i = 1:stepx+1
    for j=1:stepy+1
         a(:) = Emisfit_reshaped(i,j,:);
         % This part coded by Prof. Cai. but can't caputre min sometimes
         %zmin = fminbnd( 'Emsftz', z0, z1, [], zaxis, a);
         %Emisfitrelax(i,j)=Emsftz(zmin,zaxis,a);
         a_interp = interp1(zaxis,a,zi,'spline');
         [Emisfitrelax(i,j), I] = min(a_interp);
         Eminind(i,j)=zi(I);
         %Eminind(i,j)     =zmin;
%        if Emisfitmin(i,j) >= Emisfit
%        [Emisfitrelax(i,j), Eminind(i,j)] = min(Emisfit_reshaped(i,j,:));
    end
end

if(0)
figure(1)
mesh(Y,X,Emisfit_reshaped(:,:,1))
xlabel(['a',num2str(ind2)])
ylabel(['a',num2str(ind1)]) % primitive lattice vector a1, a2, and a3
%xlabel('dx'), ylabel('dz')
zlabel('Misfit Energy (eV/A^2)')
title(['Misfit Energy at d(a', num2str(ind3),')=0'])
end
if(1)
figure(2)
mesh(Y,X,Emisfitrelax)
%xlabel(['a',num2str(ind2)])
%ylabel(['a',num2str(ind1)]) % primitive lattice vector a1, a2, and a3
%xlabel('[1\overbar{1}0]'), 
text('Interpreter', 'latex', 'string', '$$[1 \overline{1} 0]$$', ...
    'position',[0.13 -0.7 0],'fontsize',13)
ylabel('[001]')
zlabel('\gamma_{SF} (eV/A^2)')
title('Generalized Stacking Fault Energy')
%save Mo_Gamma110 Emisfitrelax
end

if(0)
figure(3)
mesh(Y,X,Eminind)
xlabel(['a',num2str(ind2)])
ylabel(['a',num2str(ind1)]) % primitive lattice vector a1, a2, and a3
% xlabel('dx'), ylabel('dz')
zlabel(['Altitude in a',num2str(ind3),' direction'])
%zlabel('zmin (A)')
title('Altitude of atoms at which the min. energy is obtained')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the lattice vector magnitude along the shear (or Z) direction
%latticeconstant = 5.4309; % in Angstrom for Si
latticeconstant = 5.6536; % in Angstrom for Ge
%latticeconstant = 3.1472; % in Angstrom for Mo 
onelatticevector = sqrt(2)*latticeconstant; % sqrt(abs([1 0 -1]))*a
% For Mo
%onelatticevector = sqrt(3)*latticeconstant; % sqrt(abs([1 1 1]))*a

amax = onelatticevector; %
a = linspace(0,amax,length(Emisfitrelax));
misfitE_1D =zeros(1,length(Emisfitrelax));
stress = misfitE_1D;
for i=1:stepx0+1
    %misfitE_1D(i) = Emisfitrelax(i,1);
    misfitE_1D(i) = Emisfitrelax(1,i);
end

step = a(2)-a(1);
ai = a(1):amax/100:amax;
%stress = zeros(length(SurfaceE),1);
stress(1) = (misfitE_1D(2)-misfitE_1D(end-1))/step;
for i=2:length(misfitE_1D)-1
     stress(i) = (misfitE_1D(i+1)-misfitE_1D(i-1))/(2*step);
end
stress(end) = (misfitE_1D(2)-misfitE_1D(end-1))/step;
stress_interp = interp1(a, stress, ai,'spline');
[shear, I] = max(stress_interp);

%disp(['Surface Energy = ', num2str(SurfaceE(end)), '(eV/Angs^2)'])
disp(['Shear Strength = ', num2str(shear), '(eV/Angs^3)'])
disp(['The distance where the max. stress is obatained = ', num2str(ai(I)),'(Angs)'])

figure(4);
hl1 = line(a,misfitE_1D,'Color','b','LineStyle','-','marker','o');
ax1 = gca;
set(ax1,'YColor','b')
ax2 = axes('Position',get(ax1,'Position'),...
           'YAxisLocation','right',...
           'Color','none',...
           'YColor','r');
%hl2 = line(ai,stress_interp,'Color','r','Parent',ax2);
%hl3 = line(a,stress,'Color','r','linestyle','none','marker','o','Parent',ax2);
hl2 = line(ai,stress_interp*160,'Color','r','Parent',ax2);
hl3 = line(a,stress*160,'Color','r','linestyle','none','marker','o','Parent',ax2);

xlabel('distance, a (Angstrom)','fontsize',14)
%set(get(ax1,'Ylabel'),'String','(111) Misfit Energy along [10-1](eV/A^2)','fontsize',14)
set(get(ax1,'Ylabel'),'String','(-1 0 1) Misfit Energy along [1 1 1](eV/A^2)','fontsize',14)
set(get(ax2,'Ylabel'),'String','Shear Stress (GPa)','fontsize',14)
%title('Si Bulk Using Tersoff Potential','fontsize',16)
%title('Ge Bulk Using Tersoff Potential','fontsize',16)
title('Mo Bulk Using FS Potential','fontsize',16)



