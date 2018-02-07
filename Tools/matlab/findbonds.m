%analyze dislocation nucleation in Si NW

%find the slice to analyze from unstretched NW
[np0,s0,h0]=loadcn('../../runs/sinw_5D_long_seq_MEAM/Rc6_Siz/900K/relaxed2.cn');
r0= (h0*s0')';
%find all bonds

ns = length(inds);

%bondtype = 1;
%bondtype = 2;

if(~exist('nn'))
nn = zeros(ns,1);
nbrs = zeros(ns,6);
nn2=nn; nbrs2=nbrs;
rc = 2.6;

rc2=rc*rc;
for i=1:ns,
    for j=i+1:ns,
        dr = r0(inds(i),:)-r0(inds(j),:);
            drp = dr*T;
            drp2 = dr*T2;        
        if ((sum(dr.*dr)<rc2) & (abs(drp(3))>1.5) & (abs(drp(1))<0.1) & (abs(drp(2))<0.1))
            nn(i) = nn(i)+1;
            nbrs(i,nn(i)) = inds(j);
        end
        if ((sum(dr.*dr)<rc2) & (abs(drp2(3))>1.5) & (abs(drp2(1))<0.1) & (abs(drp2(2))<0.1))
            nn2(i) = nn2(i)+1;
            nbrs2(i,nn2(i)) = inds(j);
        end
    end
end
end

figure(4);
plot3(r0(inds,3),r0(inds,1),r0(inds,2),'b.');
hold on
for i=1:ns,
    for j=1:nn(i),
      plot3([r0(inds(i),3),r0(nbrs(i,j),3)],[r0(inds(i),1),r0(nbrs(i,j),1)],[r0(inds(i),2),r0(nbrs(i,j),2)],'-');
    end
end
hold off
axis equal
view([0 90]);
xlabel('z'); ylabel('x'); zlabel('y');

V=[    1.4246         0         0   -0.7123
         0    1.0318    0.0000   -0.5159
         0    0.0000   -1.0000   10.6171
         0         0         0    1.0000
];
view(V);

[np0,s0,h0]=loadcn('../../runs/sinw_5D_long_seq_MEAM/Rc6_Siz/900K/intercn900K_14p_ts20ps0001.cn');
r0=(h0*s0')';
figure(5);
%plot3(r(inds,3),r(inds,1),r(inds,2),'r.');
clf
hold on
for i=1:ns,
    for j=1:nn(i),
        dr = r(inds(i),:)-r(nbrs(i,j),:);
            drp = dr*T;
        if ((norm(drp-[0 0 2.4])>1.9)&((norm(drp+[0 0 2.4])>1.9)))
      plot3([r(inds(i),3),r(nbrs(i,j),3)],[r(inds(i),1),r(nbrs(i,j),1)],[r(inds(i),2),r(nbrs(i,j),2)],'r-');
      plot3([r0(inds(i),3),r0(nbrs(i,j),3)],[r0(inds(i),1),r0(nbrs(i,j),1)],[r0(inds(i),2),r0(nbrs(i,j),2)],'m-');
        else
%      plot3([r(inds(i),3),r(nbrs(i,j),3)],[r(inds(i),1),r(nbrs(i,j),1)],[r(inds(i),2),r(nbrs(i,j),2)],'b-');
        end
    end
end
for i=1:ns,
    for j=1:nn2(i),
        dr = r(inds(i),:)-r(nbrs2(i,j),:);
            drp2 = dr*T2;        
        if ((norm(drp2-[0 0 2.4])>1.9)&((norm(drp2+[0 0 2.4])>1.9)))
      plot3([r(inds(i),3),r(nbrs2(i,j),3)],[r(inds(i),1),r(nbrs2(i,j),1)],[r(inds(i),2),r(nbrs2(i,j),2)],'g-');
      plot3([r0(inds(i),3),r0(nbrs2(i,j),3)],[r0(inds(i),1),r0(nbrs2(i,j),1)],[r0(inds(i),2),r0(nbrs2(i,j),2)],'b-');
        else
%      plot3([r(inds(i),3),r(nbrs2(i,j),3)],[r(inds(i),1),r(nbrs2(i,j),1)],[r(inds(i),2),r(nbrs2(i,j),2)],'b-');
        end
    end
end
hold off
axis equal
view([0 90]);
xlabel('z'); ylabel('x'); zlabel('y');

V=[
    1.6084    0.0000    0.0000   -0.8042
   -0.0000    0.9986    0.0531   -0.5259
    0.0000    0.0523   -1.0129   11.2225
         0         0         0    1.0000
];    
view(V);