%analyze dislocation nucleation in Si NW

%find the slice to analyze from unstretched NW
[np0,s0,h0]=loadcn('../../runs/sinw_5D_long_seq_MEAM/Rc6_Siz/900K/relaxed2.cn');
r0= (h0*s0')';
%find atoms in the slice
ex = [0 0 1]; ey = [1 -1 0]; ez = cross(ex,ey); 
ex = ex/norm(ex); ey = ey/norm(ey); ez = ez/norm(ez);
ez1 = [1 1 1]; ex1 = [-1 -1 2]; ey1 = cross(ez1,ex1);
%ez1 = [1 1 1]; ex1 = [-1 0 1]; ey1 = cross(ez1,ex1);
ex1 = ex1/norm(ex1); ey1 = ey1/norm(ey1); ez1 = ez1/norm(ez1);
T = [dot(ex,ex1), dot(ex,ey1), dot(ex,ez1)
     dot(ey,ex1), dot(ey,ey1), dot(ey,ez1)
     dot(ez,ex1), dot(ez,ey1), dot(ez,ez1)];
r0p = r0*T;

ez1 = [1 1 -1]; ex1 = [1 1 2]; ey1 = cross(ez1,ex1);
ex1 = ex1/norm(ex1); ey1 = ey1/norm(ey1); ez1 = ez1/norm(ez1);
T2 = [dot(ex,ex1), dot(ex,ey1), dot(ex,ez1)
     dot(ey,ex1), dot(ey,ey1), dot(ey,ez1)
     dot(ez,ex1), dot(ez,ey1), dot(ez,ez1)];

%inds = find(abs(r0p(:,3)-130)<2.3);
%inds = find(abs(r0p(:,3)+176)<8);
inds = find(abs(r0p(:,3)+177)<15);
%inds = find((abs(r0p(:,3)+177)<15)&(r0(:,2)>0));


%if(~exist('s'))
  [np0,s0,h0]=loadcn('../../runs/sinw_5D_long_seq_MEAM/Rc6_Siz/900K/intercn900K_14p_ts20ps0001.cn');
  [np, s ,h ]=loadcn('../../runs/sinw_5D_long_seq_MEAM/Rc6_Siz/900K/intercn900K_14p_ts20ps0040.cn');
%  [np,s,h]=loadcn('../../runs/sinw_5D_long_seq_MEAM/Rc6_Siz/900K/intercn900K_14p_ts20ps0030.cn');
%end

%compute displacement between two frames
%ind0 = 11720;
ds = s-s0;  
%ds=ds-ones(length(s(:,1)),1)*(s(ind0,:)-s0(ind0,:));
dr = (h*ds')';

r = (h*s')';
h0(3,3)=h0(3,3)*1.0;
r0= (h0*s0')';

R0 = sqrt(r0(:,1).*r0(:,1)+r0(:,2).*r0(:,2));
drnorm = sqrt(sum(dr.*dr,2));

%dravg=sum(r(inds,:)-r0(inds,:),1)/length(inds);
%ind0=14149;
%dravg=r(ind0,:)-r0(ind0,:);
%dravg = [0 0 0.2];
dravg = [0 0 0];

if(1)
V = [    1.4210         0         0   -0.7105
         0    0.9994    0.0361   -0.5178
         0    0.0349   -1.0342   10.6119
         0         0         0    1.0000
];
else
V = [    0.0000   -1.8349   -0.0000    0.9174
   -0.0511   -0.0000    0.9998   -0.4744
    2.9283    0.0000    0.0175   16.5161
         0         0         0    1.0000
];
end

figure(3);
plot3(r0(inds,3),r0(inds,1),r0(inds,2),'b.');
axis equal
view([0 90]);
xlabel('z'); ylabel('x'); zlabel('y');
view(V);

figure(2);
plot3(r(inds,3),r(inds,1),r(inds,2),'r.');
axis equal
view([0 90]);
xlabel('z'); ylabel('x'); zlabel('y');
view(V);

figure(1);
plot3(r0(inds,3),r0(inds,1),r0(inds,2),'b.');
%plot3(r0(inds,3),r0(inds,1),r0(inds,2),'b.',r(inds,3),r(inds,1),r(inds,2),'r.');
r = r-ones(np,1)*dravg;
%clf
hold on
for i=1:length(inds),
    plot3([r0(inds(i),3),r(inds(i),3)],[r0(inds(i),1),r(inds(i),1)],[r0(inds(i),2),r(inds(i),2)],'r-');
end
hold off
axis equal
view([0 90]);
xlabel('z'); ylabel('x'); zlabel('y');
view(V);


return;


























if(0)
figure(1);
%hist(drnorm,[0:0.02:4]);
hist(drnorm,50);

figure(2);
%ind = find((drnorm>1)&(R0<21)&(r0(:,3)>100)&(r0(:,3)<200));
ind = find((drnorm>0.8)&(R0<21)&(r0(:,3)>100)&(r0(:,3)<200));
length(ind)
%plot3(r0(ind,1),r0(ind,2),r0(ind,3),'b.',r(ind,1),r(ind,2),r(ind,3),'r.');
%plot3(r0(ind,3),r0(ind,1),r0(ind,2),'b.',r(ind,3),r(ind,1),r(ind,2),'r.');
plot3(r0(ind,3),r0(ind,1),r0(ind,2),'b.');
hold on
for i=1:length(ind),
%    plot3([r0(ind(i),1),r(ind(i),1)],[r0(ind(i),2),r(ind(i),2)],[r0(ind(i),3),r(ind(i),3)],'k-');
    plot3([r0(ind(i),3),r(ind(i),3)],[r0(ind(i),1),r(ind(i),1)],[r0(ind(i),2),r(ind(i),2)],'k-');
end
hold off
axis equal

%view([0 0]);
view([90 0]);
%view([90 90]);
end