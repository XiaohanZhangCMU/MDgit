% find dislocation velocity from cfg files
% example:
%   filepath = '../../runs/ta-screw-md-300K-600MPa';
%   spec = [ 20 1020 100 98.77 100*1e-3 2.0 10 -1 1 0.5];
%   find_disl_vel

if (~exist('filepath'))
  filepath = '../../runs/ta-edge-10';
  %filepath = '../../runs/ta-edge-100';
  %filepath = '../../runs/ta-71-md-100';
end

if (~exist('spec'))
  spec = [ 10 410 10 74.34 1000*1e-3 4.5 10 0.5 0.516 0.7];
end

nstart = spec(1);
nend   = spec(2);
nstep  = spec(3);
x0     = spec(4);      %initial dislocation position
%x0 = 98;              %initial dislocation position
x0 = 74.34;            %initial dislocation position
frameno = nstart:nstep:nend;  nframe = length(frameno);
t = frameno * spec(5);  % time in ps
csdmin = spec(6);
csdmax = spec(7);
symin  = spec(8);
symax  = spec(9); 
avgstart=spec(10);

xavg = zeros(size(frameno));

for i=1:nframe,
    % load atom configuration
    ifile = frameno(i);
    filename = sprintf('%s/auto%04d.cfg',filepath,ifile);
    [np,data,h]=loadcfg(filename);
    
    % obtain scaled coordinates and cental symmetry paramter
    s = data(:,1:3); csd = data(:,4); 

    % view histogram of central symmetry parameter
    %figure(3); hist(csd);

    % select dislocation core atoms
    ind = find(csd>csdmin & csd<csdmax & s(:,2) > symin & s(:,2) < symax);
    %ind = find(csd>4.5 & csd<10);
    s = data(ind,1:3); r = s*h'; 

    % make sure all atoms are nearest images
    sx = s(:,1); sx = sx - round(sx-sx(1)); s(:,1) = sx; r = s*h';
    
    [sz ind] = sort(s(:,3)); s = s(ind,:); r = r(ind,:);
    
    % compute average dislocation position
    xavg(i) = mean(r(:,1)) - x0;
    if(i>1) xavg(i) = xavg(i) - round((xavg(i)-xavg(i-1))/h(1,1))*h(1,1); end

    % remove debris (left by screw dislocations)
    ind = find(abs(r(:,1)-mean(r(:,1))) < 20);
    if length(ind) < length(r(:,1)),
       r = r(ind,:);
       xavg(i) = mean(r(:,1)) - x0;
       if(i>1) xavg(i) = xavg(i) - round((xavg(i)-xavg(i-1))/h(1,1))*h(1,1); end
    end
        
    % plot dislocation position
    figure(1); plot(r(:,3),r(:,1),'*-'); axis equal
    xlim([0 h(3,3)]); ylim([0 h(1,1)]);
    title(sprintf('load %s',filename));
    drawnow; pause(0.2);
end

% estimate dislocation velocity 
%  (from second half of data)
nmid = round(nframe*avgstart);
vavg = (xavg(nframe) - xavg(nmid)) / ...
       (t(nframe) - t(nmid));  % in A/ps
vavg_in_m_s = vavg * 100;  % in m/s

fs = 12;
figure(2);
plot(t, xavg, '.-', ...
     t(nmid:end), (t(nmid:end)-t(nmid))*vavg+xavg(nmid),'r--');
set(gca,'FontSize',fs);
xlabel('t (ps)');
ylabel('x (angstrom)');
title(sprintf('vavg = %.1f (m/s)',vavg_in_m_s));
