%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: free_eng_ha.m
%%
%% compute free energy based on harmonic approximation
%%  from hessian.out
%%
%% Last Modified : Wed Jan 24 16:54:36 2007
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input file
%datafile = '../../runs/si-perf-ha/hessian.out';
%datafile = '../../runs/si-core-ha/hessian.out';
%datafile = '../../runs/si-core-vac-ha/hessian.out';
datafile = '../../runs/si-Hessian-2/hessian.T1700.out';

% parameters
T = 1700;     % in K
m = 28.0855; % in (g/mol), Si

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('\nload file %s',datafile));
data = load(datafile);
N = sqrt(length(data(:,1))/3);
disp(sprintf('number of atoms: N = %d',N));
data = reshape(-data',N*3,N*3);
data = (data+data')/2;  % symmetrize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('diagonalizing...');
Lambda = eig(data);
Lambda = sort(Lambda);
eps = 1e-10;
Lambda = (abs(Lambda)>eps).*Lambda; %in eV/A^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AVO  = 6.0221e+23;
EV   = 1.6022e-19;
BOLZ = 1.381e-23 ;        % (J/K)  
KB   = 8.6173e-5 ;        % (eV/K) 
hbar = 1.0546e-34;        % (m^2 kg / s)

mass = m*1e-3/AVO;        % convert to kg
Lam  = Lambda*EV*1e20;    % convert to J/m^2
omega= sqrt(Lam/mass);    % (1/s)
hw   = hbar*omega;        % (J)
KT   = BOLZ*T;            % (J)

hwoverKT = hw/KT;

%remove zero entries
hwoverKT = (abs(hwoverKT)>eps).*hwoverKT + (abs(hwoverKT)<eps).*1;

F = KB*T*sum(log(hwoverKT));

disp(sprintf('F = %f (eV)  T = %.1f (K)',F,T));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
