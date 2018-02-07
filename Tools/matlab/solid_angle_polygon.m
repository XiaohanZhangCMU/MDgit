function Omega = solid_angle_polygon(rn,xp,yp,zp)
% compute the solid angle of a polygon rn(x1,y1,z1; ...; xn,yn,zn)
% seen by a field point rp = (xp,yp,zp)
%
% to test:
%  x=[-2:0.1:2]; rn=rand(4,3); rn(:,3)=0;
%  plot(x,solid_angle_polygon(rn,x,0.3,0.1))

[Nnodes,M] = size(rn);
if (Nnodes<3)
    disp('solid_angle_polygon: rn should have at least 3 rows!');
    return;
end
if (M~=3)
    disp('solid_angle_polygon: rn should have 3 columns!');
    return;
end

%Omega=solid_angle_triangle(rn(1,:),rn(2,:),rn(3,:),xp,yp,zp);
Omega=SolidAngleTriangle(rn(1,1),rn(1,2),rn(1,3),...
                           rn(2,1),rn(2,2),rn(2,3),...
                           rn(3,1),rn(3,2),rn(3,3),...
                           xp,yp,zp);

for i=4:Nnodes,
%   Omega=Omega+solid_angle_triangle(rn(1,:),rn(i-1,:),rn(i,:),xp,yp,zp);
   Omega=Omega+SolidAngleTriangle(rn(1,1),rn(1,2),rn(1,3),...
                                    rn(i-1,1),rn(i-1,2),rn(i-1,3),...
                                    rn(i,1),rn(i,2),rn(i,3),...
                                    xp,yp,zp);
end
