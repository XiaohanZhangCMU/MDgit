%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: plotEmisfit.m
%%
%% Last Modified : Tue Oct 10 18:07:19 2006
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relevant information from MD++ script
dirname = '../../../../runs/al-misfit';

latticesize = [ 1 1  2  4
               -1 1  0  4
	        1 1 -1  4 ]

input = [ 3             ...  % surface normal is z direction
          0 0.01  0.50  ...  % x0,dx,x1
          0 0.01  0     ...  % y0,dy,y1
      -0.01 0.01  0.20  ...  % z0,dz,z1
        ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in MD++ output data
Edata = load(sprintf('%s/Emisfit.out',dirname));
Mx=max(Edata(:,1));
My=max(Edata(:,2));
Mz=max(Edata(:,3));

x0 = input(2); dx = input(3); x1 = input(4);
y0 = input(5); dy = input(6); y1 = input(7);
z0 = input(8); dz = input(9); z1 = input(10);

E0=min(Edata(:,4)); % reference energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find relaxed energy with respect to z
%
% surface energy as a function of x and z
E = reshape(Edata(:,4),Mx,Mz);

% after relax wrt z, obtain surface energy as a function of x
Erlx = zeros(Mx,1);
Zrlx = zeros(Mx,1);

% construct a very fine grid to find minimum
finegrid = [z0:dz/1000:z1];

for i = 1:Mx,
  E_interp = interp1([z0:dz:z1],E(i,:),finegrid,'spline');
  [Erlx(i), ind] = min(E_interp);
  Zrlx(i) = finegrid(ind);
end
Erlx = Erlx - E0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot relaxed energy
fs = 17;
figure(1);
plot([x0:dx:x1], Erlx);
set(gca,'FontSize',fs);
xlabel('x');
ylabel('E_{relaxed}  (eV/A^2)');

figure(2);
plot([x0:dx:x1], Zrlx);
set(gca,'FontSize',fs);
xlabel('x');
ylabel('z_{min}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
