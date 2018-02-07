function Omega = solid_angle_triangle(r1,r2,r3,xp,yp,zp)
% compute the solid angle of a triangle r1:(x1,y1,z1)-r2:(x2,y2,z2)-r3:(x3,y3,z3)
% seen by a field point rp = (xp,yp,zp)
%
% to test:
%  x=[-2:0.1:2]; plot(x,solid_angle_triangle([-1 -1 0],[1 -1 0],[0 1 0],x,0.3,0.1))

nx=length(xp); ny=length(yp); nz=length(zp);
N=max([nx,ny,nz]);

if ((nx~=ny)&(nx~=1)&(ny~=1)) | ...
   ((ny~=nz)&(ny~=1)&(nz~=1)) | ...
   ((nz~=nx)&(nz~=1)&(nx~=1))
    disp('solid_angle_triangle: x,y,z should have either the same or unit length!');
    return;
end

Omega=zeros(N,1);
for i=1:N,  
  xs = xp( min(i,nx) );
  ys = yp( min(i,ny) );
  zs = zp( min(i,nz) );
  
  R1 = [r1(1)-xs,r1(2)-ys,r1(3)-zs];
  R2 = [r2(1)-xs,r2(2)-ys,r2(3)-zs];
  R3 = [r3(1)-xs,r3(2)-ys,r3(3)-zs];
  
  RR1 = norm(R1); RR2 = norm(R2); RR3 = norm(R3);
  numer = det([R1;R2;R3]);
  denom = RR1*RR2*RR3 + dot(R1,R2)*RR3 + dot(R2,R3)*RR1 + dot(R3,R1)*RR2;
  
  Omega(i) = -2 * atan2(numer, denom);
end