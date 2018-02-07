%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: plot_switch.m
%%
%% compute free energy difference from MC adiabatic switching
%%  simulation data (stored in prop.out)
%%
%% Last Modified : Wed Aug 30 14:44:29 2006
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input file
%datafile = '../../runs/si-core-vac-mc/prop.out';
datafile = '../../runs/si-core-vac-sw2ha/prop.out';
%datafile = '../../runs/si-perf-sw2ha/prop.out';
%datafile = '../../runs/si-perf-mc/prop.out';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switching parameters (need to be consistent with MD++ script)
totalsteps = 500000;
savepropfreq = 100;
lambda0 = 1;
lambda1 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('\nload file %s',datafile));
data = load(datafile);

N = totalsteps/savepropfreq+1;    %number of data per switch
Nsw = floor(length(data(:,1))/N); %number of switches
disp(sprintf('Number of switches: Nsw = %d',Nsw));
steps = data(1:N,1);
dHdl  = reshape(data(1:N*Nsw,12),N,Nsw);
if lambda0>lambda1
    dHdl = -dHdl;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = [steps(1):steps(2):totalsteps]'/totalsteps;
t2 = t .*t;
t3 = t2.*t;
t4 = t3.*t;
t5 = t4.*t;
lam  = (t5.*(70.*t4-315.*t3+540.*t2-420.*t+126.));
dldt = (t4.*(630.*t4-2520.*t3+3780.*t2-2520.*t+630.));

ind = steps/savepropfreq + 1;
dHdt= zeros(N,Nsw);
dW  = zeros(1,Nsw);
for j=1:Nsw,
    dHdt(:,j) = dHdl(:,j).*dldt(ind);
    dW(j) = (sum(dHdt(:,j))-0.5*(dHdt(1,j)+dHdt(end,j)))/N;
end
disp(sprintf('dW = %.3f +- %.3f (eV)',mean(dW),std(dW)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(t,lam,'b-',t,dldt,'r-');
set(gca,'FontSize',19);
xlabel('s');

figure(2);
colors = ['b.';'g.';'r.';'c.';'m.';'y.';'k.'];
plot(t(ind),dHdt(:,Nsw),colors(mod(Nsw-1,length(colors))+1,:));
hold on
for j=Nsw-1:-1:1,
    plot(t(ind),dHdt(:,j),colors(mod(j-1,length(colors))+1,:));
end
hold off
set(gca,'FontSize',24);
xlabel('s');
ylabel('dH/ds');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
