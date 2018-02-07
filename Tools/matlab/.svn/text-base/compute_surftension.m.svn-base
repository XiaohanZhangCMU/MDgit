%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: surfacetension.m
%%
%% compute surface tension
%%
%% Last Modified : Fri Nov 26 2010
%% Seunghwa Ryu, shryu@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% inputs
file1='PtPn_info.out';
file2='Pnorm.out';
file3='Ptran.out';
file4='slot_coordinate.out';
basicinfo=load(file1);
Pnorm_=load(file2);
Ptran_=load(file3);
coordinate=load(file4);

height = basicinfo(2);
slotvolume = basicinfo(3);
timestep = basicinfo(4);
coordinate=(coordinate+0.5)*height;
distance=0.0:height/1000:height;

% divide by slotvolume, it is raw data in eV
Pnorm=zeros(1,length(coordinate)+1);
Ptran=zeros(1,length(coordinate)+1);
Pnorm(1:end-1)=mean(Pnorm_(:,2:end)/slotvolume);
Pnorm(end)=Pnorm(1);
Ptran(1:end-1)=mean(Ptran_(:,2:end)/slotvolume);
Ptran(end)=Ptran(1);
tau=Pnorm-Ptran; % tau is surface tension of whole thin film = 2*surface energy
Pnorm_fit=spline([coordinate height],Pnorm,distance);
Ptran_fit=spline([coordinate height],Ptran,distance);
tau_fit=Pnorm_fit-Ptran_fit;

% in MPa from eV/A^3
Pnorm_inMPa=Pnorm*160.2e3; 
Ptran_inMPa=Ptran*160.2e3;
tau_inMPa=Pnorm_inMPa-Ptran_inMPa;
Pnorm_inMPa_fit=spline([coordinate height],Pnorm_inMPa,distance);
Ptran_inMPa_fit=spline([coordinate height],Ptran_inMPa,distance);
tau_inMPa_fit=Pnorm_inMPa_fit-Ptran_inMPa_fit;

figure(1);
plot([coordinate height],Pnorm_inMPa,'o',distance,Pnorm_inMPa_fit);
xlabel('distance (A)');
ylabel('Pn (MPa)');

figure(2);
plot([coordinate height],Ptran_inMPa,'o',distance,Ptran_inMPa_fit);
xlabel('distance (A)');
ylabel('Pt (MPa)');

figure(3);
plot([coordinate height],tau_inMPa,'o',distance,tau_inMPa_fit);
xlabel('distance (A)');
ylabel('Tau (MPa)');

figure(4);
plot(Pnorm_(:,1),Pnorm_(:,round(length(coordinate)/2)),'o-')

TAU_whole=trapz(distance,tau_inMPa_fit)*10^(-4);
Surf_energy=TAU_whole/2;
disp(sprintf('surface energy = %f (J/m^2)',Surf_energy));

%% correlation time part will be implemented later (Fir Nov 26 2010) %%
