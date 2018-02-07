%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: plotEmisfit_manyplane.m
%%
%% Last Modified : Tue Oct 11 14:59:48 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load mo_mp_data.mat
%Sdata=max(abs(Sdata),[],1);
[ms,ns]=size(Sdata);

%fit data
if(strcmp(executable, 'fs_gpp'))
  plot(dist,Sdata(:,1),'bd');  
  x=[0:0.01:1];
  f=Sdata(1,1)+sin(x*pi)*(max(Sdata(:,1))-Sdata(1,1));
  hold on
  plot(x*max(dist),f,'g-');
  if(ns>2)
    f=Sdata(1,3)+sin(x*pi)*(max(Sdata(:,3))-Sdata(1,3));
    plot(x*max(dist),f,'b-');
    plot(dist,Sdata(:,3),'ro');
  end
  %from Richard Baran
  olddata=[24.7355 28.5633 28.0482 27.4421 26.9891 26.6561 26.4052 26.2106 26.0559 25.9300 25.8259 25.2060];
  plot(dist,olddata/2,'x');
  hold off
  legend('MD++ (relaxed)','MD++ (unrelaxed)','old data by Richard');
  title('BCC Finnis-Sinclair Mo Ideal Shear Strength');
  xlabel('d (A)');
  ylabel('\sigma (GPa)');
elseif(strcmp(executable, 'alglue_gpp'))
  M = ones(size(N));
  vy_all = M'*[-1 -1 1] + N'*[-1 -1 -1];
  dist2  = latticeconst./sqrt(sum(vy_all.*vy_all,2));
  if(ns>2)
    plot(dist2,Sdata(:,1),'b+',dist2([1:5,end]),Sdata([1:5,end],1),'gd',dist2(3:end),Sdata(3:end,1),'b-');
    hold on
    ind = 3;
  else
    ind=1;
  end
  plot(dist2,Sdata(:,ind),'b.',dist2([1:5,end]),Sdata([1:5,end],ind),'ro',dist2(3:end),Sdata(3:end,ind),'b-');
  hold off
  t1 =text(dist2(1),Sdata(1,ind)-0.2,'n=-1');
  t1a=text(dist2(1),Sdata(1,ind)-0.4,'[0 0 1]');

  t2 =text(dist2(2),Sdata(2,ind)-0.2,'n=0');
  t2a=text(dist2(2),Sdata(2,ind)-0.4,'[-1 -1 1]');

  t3 =text(dist2(3),Sdata(3,ind)-0.2,'n=1');
  t3a=text(dist2(3),Sdata(3,ind)-0.4,'[-1 -1 0]');

  t4 =text(dist2(4),Sdata(4,ind)-0.2,'n=2');
  t4a=text(dist2(4),Sdata(4,ind)-0.4,'[-3 -3 -1]');

  t5 =text(dist2(5),Sdata(5,ind)-0.2,'n=3');
  t5a=text(dist2(5),Sdata(5,ind)-0.4,'[-2 -2 -1]');

  te =text(dist2(end),Sdata(end)-0.2,sprintf('n=%d',N(end)));

  
  xlim([-0.05 3]);
  ylim([4  11.5]);

  title('FCC Ercolessi-Adams Al Ideal Shear Strength');
  xlabel('d (A)');
  ylabel('\sigma (GPa)');

  for i=1:length(N),
    disp(sprintf(' %4d  %14.9f %14.9f',N(i),Sdata(i,1),Sdata(i,3)));
  end
else
  plot(N,Sdata);
end
