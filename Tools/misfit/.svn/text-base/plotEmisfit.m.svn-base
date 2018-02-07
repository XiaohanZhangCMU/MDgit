%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: plotEmisfit.m
%%
%% Last Modified : Mon Oct 17 15:02:53 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load mo_data

if (sum([Mx,My,Mz]>1)==1)
  disp('E is a 1-dimensional array');  
  E = Edata(:,4);
  figure(1);
  plot([0:dx:dxyz(1,3)],E,'.-');
  xlabel('dx');
  ylabel('E (eV/A^2)');
  %find maximum slope dE/dx
  sigmax = max(E(2:end)-E(1:end-1))/(dx*Lx) * 160.2;
  disp(sprintf('sigmax = %e (GPa)',sigmax));
  
elseif (sum([Mx,My,Mz]>1)==2)
  disp('E is a 2-dimensional array');  
  MM=-sort(-[Mx,My,Mz]);
  E = reshape(Edata(:,4),MM(1),MM(2));
  figure(1);
  mesh(E);
  %figure(2);
  %contour(E,30);

  if(exist('Erlx'))
    figure(2);
    plot(Erlx-Erlx(1),'r.-');
    hold on
    plot(Eunrlx-Eunrlx(1),'b-');    
    hold off
    legend('relaxed','unrelaxed');
    figure(3);
    plot([1:Mx],Zrlx,'r.-',[1 Mx],[1 1]*dxyz(normaldir,1),'b-',[1 ...
                    Mx],[1 1]*dxyz(normaldir,3),'b-');

    sigmax = max(Erlx(2:end)-Erlx(1:end-1))/(dx*Lx) * 160.2;
    sigmin = min(Erlx(2:end)-Erlx(1:end-1))/(dx*Lx) * 160.2;

    sigmax0 = max(Eunrlx(2:end)-Eunrlx(1:end-1))/(dx*Lx) * 160.2;
    sigmin0 = min(Eunrlx(2:end)-Eunrlx(1:end-1))/(dx*Lx) * 160.2;

    disp(sprintf('sigmax = %e   sigmax0 = %e (GPa)',sigmax,sigmax0));
    
  end    
else
  E = reshape(Edata(:,4),Mx,My,Mz);
  disp('E is a 3-dimensional array');
end


