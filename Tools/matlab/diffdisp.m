bccfile  = 'bcc3x5x1.cn';
dislfile = 'quad3x5x1-unrelaxed.cn';
rxy2_min  = 1;
rxy2_max  = 1.2*3.1472^2.*( norm([1 1 0]./2)^2 + norm([1 1 2]./3)^2);
ArrowScale = 30;
% Min and Max of (R_xy)^2 in determining whether to plot the arrow

[np, s, h] =loadcn(bccfile);
[np1,s1,h1]=loadcn(dislfile);

r=s*h; r1=s1*h1;
zd=r1(:,3)-r(:,3);

B=[1 1 1]./2.*3.1472; b=norm(B);

% Plot atoms in 3-D but view in X-Y plane.
clf;  plot(r(:,1),r(:,2),'ro','MarkerSize',2,'LineWidth',5)
axis equal; grid off;
xlabel('X   [112]','FontSize',20); ylabel('Y   [110]','FontSize',20)

t=[];  pstart=[]; pend=[];
% Calculate displacement vectors between atoms and draw arrow lines.
for i=1:(np-1),
  for j=(i+1):np,
    xd=r1(i,1)-r1(j,1); yd=r1(i,2)-r1(j,2); zdd=zd(i)-zd(j); ...
    zdd=zdd-round(zdd./b).*b;
    rxy2 = xd^2+yd^2;
    if (rxy2<rxy2_max) & (rxy2>rxy2_min),
      if zdd > 0,
        pstart = [ pstart; r1(i,1:2)];
        pend   = [ pend;   r1(j,1:2)];
      else,
        pstart = [ pstart; r1(j,1:2)];
        pend   = [ pend;   r1(i,1:2)];
      end
      t = [ t; abs(zdd) ];
    end
  end
end

du_max=max(t),
t=t./du_max;  ta=(1+t)./2; tb=1-ta;
for i = 1:length(t),
  Pstart = pstart(i,:).*ta(i) + pend(i,:).*tb(i);
  Pend = pstart(i,:).*tb(i) + pend(i,:).*ta(i);
  arrow(Pstart,Pend,ArrowScale*t(i));
end;
