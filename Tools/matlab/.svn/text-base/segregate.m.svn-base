% determine whether or not Y atom segregates at GB

istart = 30;
iend   = 50;
dirname = 'yszruns/c';

fs = 17;

%first frame
iframe = 1;
cnfilename = sprintf('../../runs/%s/inter%04d.cn',dirname,iframe);
[np,sext,h]=loadcn2(cnfilename);
species = sext(:,10);
ind = find (species==3);  % Y atoms
sY0 = sext(ind,1:3);
rY0 = (h*sY0')';

%average over a set of later frames
if(1)
 sY = [];
 species = [];
 for iframe = istart:iend,
    cnfilename = sprintf('../../runs/%s/inter%04d.cn',dirname,iframe);
    [np,sext,h]=loadcn2(cnfilename);
    species = sext(:,10);
    ind = find (species==3);  % Y atoms
    
    if length(sY)==0
        sY = sext(ind,1:3);
    else
        sY = [sY; sext(ind,1:3)];
    end
 end
end    
rY = (h*sY')';

figure(1);
subplot(2,1,1);
%plot3(rY(:,1),rY(:,2),rY(:,3),'b.'); %view([0 45]);
plot(rY(:,1),1.0*rY(:,2)+0.0*rY(:,3),'b.');
 hold on; plot([-15 -15],[-15 15],'r-'); plot([15 15],[-15 15],'r-'); hold off
axis equal 
set(gca,'FontSize',fs);
%xlabel('x');
ylabel('y');
%zlabel('z');
xlim([-40 40]);
title(sprintf('%s  from %d to %d',dirname,istart,iend));

subplot(2,1,2);
%hist(rY(:,1),[-40:2.5:40]);
% hold on; plot([-15 -15],[0 10],'r-'); plot([15 15],[0 10],'r-'); hold off
%xlim([-40 40]);
%figure(2);
x = [-40:0.1:40];  f = zeros(size(x)); f0 = f;
a = 0.5;  %half-width of the Gaussian distribution

%first frame
for i=1:length(rY0(:,1))
    f0 = f0+exp(-(x-rY0(i,1)       ).^2/2/a^2);
    f0 = f0+exp(-(x-rY0(i,1)-h(1,1)).^2/2/a^2);
    f0 = f0+exp(-(x-rY0(i,1)+h(1,1)).^2/2/a^2);
end
f0 = f0/sqrt(2*pi)/a * ((max(rY(:,1))-min(rY(:,1)))/length(rY0(:,1)));
for i=1:length(rY(:,1))
    f = f+exp(-(x-rY(i,1)       ).^2/2/a^2);
    f = f+exp(-(x-rY(i,1)-h(1,1)).^2/2/a^2);
    f = f+exp(-(x-rY(i,1)+h(1,1)).^2/2/a^2);
end
f = f/sqrt(2*pi)/a * ((max(rY(:,1))-min(rY(:,1)))/length(rY(:,1)));
p=plot(x,f0,'r-',x,f,'b-',[-15 -15],[0 2],'r--', ...
    [15 15],[0 2],'r--',x,ones(size(x)),'b--');
set(p(2),'LineWidth',2);
set(gca,'FontSize',fs);
xlim([-40 40]);
xlabel('x');
ylabel('atom density');
leg=legend('initial','final'); set(leg,'FontSize',12');

%print -djpeg -r150 ysz-b