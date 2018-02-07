%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: diffusion.m
%%
%% Last Modified : Wed Aug 30 14:44:29 2006
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


istart = 10;
iend   = 70;
%dirname = 'si_au_md_T0600';
%dirname = 'si_au_md_T0700';
%dirname = 'si_au_md_T0800';
%dirname = 'si_au_md_T1000';
%dirname = 'si_au_md_quench_3';
%dirname = 'si-example';

%dirname = 'si_au_nve_T0830_s0044';
%dirname = 'si_au_nve_T0790_s0044';
%dirname = 'si_au_nve_T0700_s0044';
%dirname = 'si_au_nve_T0610_s0044';
%dirname = 'si_au_nve_T0540_s0044';
dirname = 'si_au_nve_T0460_s0044';


savecnfreq = 2500;
dt = 1e-4; %ps

fs = 17;

%first frame
propfilename = sprintf('../../runs/%s/prop.out',dirname);
prop = load(propfilename);
Ekin=prop(:,2); Epot=prop(:,3); Etot=Ekin+Epot; P=prop(:,4); T=prop(:,10); 
HP=prop(:,9); N=length(prop(:,1));
figure(1);
subplot(3,1,1);
plot([0:N-1],T,'-');
ylabel('T');
ylim([0 1800]);
grid
subplot(3,1,2);
%plot([0:N-1],prop(:,2),[0:N-1],prop(:,3),[0:N-1],prop(:,2)+prop(:,3));
%plot([0:N-1],Etot-Etot(1));
plot([0:N-1],HP-HP(1));
ylabel('E');
subplot(3,1,3);
plot([0:N-1],P*160);
ylabel('P (GPa)');

%return;

iframe = istart;
cnfilename = sprintf('../../runs/%s/inter%04d.cn',dirname,iframe);
[np,sext,h]=loadcn2(cnfilename);
species = sext(:,10);
ind = find (species==0);  % Au atoms
sAu0 = sext(ind,1:3);
rAu0 = (h*sAu0')';
ind = find (species==1);  % Si atoms
sSi0 = sext(ind,1:3);
rSi0 = (h*sSi0')';
%[np,s0,h]=loadcn(cnfilename);
%r0 = (h*s0')';

msdSi = zeros(iend-istart+1,1);
msdAu = zeros(iend-istart+1,1);

for iframe = istart+1:iend
  cnfilename = sprintf('../../runs/%s/inter%04d.cn',dirname,iframe);
  [np,sext,h]=loadcn2(cnfilename);
  species = sext(:,10);
  ind = find (species==0);  % Au atoms
  sAu1 = sext(ind,1:3);
  rAu1 = (h*sAu1')';
  ind = find (species==1);  % Si atoms
  sSi1 = sext(ind,1:3);
  rSi1 = (h*sSi1')';
  %[np,s1,h]=loadcn(cnfilename);
  %r1 = (h*s1')';

  drSi = rSi1-rSi0;
  drSi2 = sum(drSi.*drSi,2);
  msdSi(iframe-istart+1) = mean(drSi2);
  drAu = rAu1-rAu0;
  drAu2 = sum(drAu.*drAu,2);
  msdAu(iframe-istart+1) = mean(drAu2);

end

%figure(1);
%plot3(r0(:,1),r0(:,2),r0(:,3),'bo');
%axis equal
%
%figure(2);
%plot3(r1(:,1),r1(:,2),r1(:,3),'ro');
%axis equal

%figure(3);
%plot(sort(dr2));
%hist(dr2, [0:10:400]);
%ylim([0 100]);

figure(4)
td = [0:length(msdSi)-1]*savecnfreq*dt;
skip=10;
%
pSi=polyfit(td(skip:end)',msdSi(skip:end),1);  DSi=pSi(1);
pAu=polyfit(td(skip:end)',msdAu(skip:end),1);  DAu=pAu(1);
%
%DSi = (msdSi(end)-msdSi(skip))/(td(end)-td(skip));
%DAu = (msdAu(end)-msdAu(skip))/((length(msdAu)-skip)*savecnfreq*dt);
plot(td,msdSi, 'b.-', td,msdAu, 'r.-', ...
     td,pSi(2)+td*DSi,'b--', ...
     td,pAu(2)+td*DAu,'r--' );
set(gca,'FontSize',17);
legend('Si','Au',2);
xlabel('t  (ps)');
ylabel('msd (A^2)');

disp(sprintf(' T_avg = %f K  P_avg = %f GPa',mean(T(skip:end)),mean(P(skip:end))*160));
disp(sprintf('D_Si = %.5f       D_Au = %.5f     (A^2/ps)',DSi,DAu));
disp(sprintf('D_Si = %.4e   D_Au = %.4e (cm^2/s)',DSi*1e-4,DAu*1e-4));
disp(sprintf('D_Si = %.4e   D_Au = %.4e (m^2/s)\n',DSi*1e-8,DAu*1e-8));

disp(sprintf('    %f %f     %.4e  %.4e',mean(T(skip:end)),mean(P(skip:end))*160,DSi*1e-8,DAu*1e-8));

