%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: calmsd.m
%%
%% Last Modified : Tue Oct 10 16:05:06 2006
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nstart = 1;
nend = 646;

data = zeros(size([nstart:nend]),3);

id = 10;

filename = sprintf('../../runs/ar-md/perf.cn');
[np,s0,h] = loadcn(filename);
r0 = (h*s0')';

for i=nstart:nend,
  filename = sprintf('../../runs/ar-md/inter%04d.cn',i);
  %filename
  [np,s,h] = loadcn(filename);
  r = (h*s')';
  data(i,:) = r(id,:)-r0(id,:);
end
  
ravg = [mean(data(:,1)), mean(data(:,2)), mean(data(:,3))];
rstd = [ std(data(:,1)),  std(data(:,2)),  std(data(:,3))];

ravg
rstd
