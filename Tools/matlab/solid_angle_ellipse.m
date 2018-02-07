function Omega = solid_angle_ellipse(a,b,x,y,z)
% compute the solid angle of an ellipse (x/a)^2+(y/b)^2=1
% seen by a field point r(3) = [x,y,z]
% need to perform an integral numerically
%
% to test:
%  x=[-2:0.1:2]; plot(x,solid_angle_ellipse(1,2,x,0.3,0.1))

tol = 1e-3;

nx=length(x); ny=length(y); nz=length(z);
N=max([nx,ny,nz]);

if ((nx~=ny)&(nx~=1)&(ny~=1)) | ...
   ((ny~=nz)&(ny~=1)&(nz~=1)) | ...
   ((nz~=nx)&(nz~=1)&(nx~=1))
    disp('solid_angle_ellipse: x,y,z should have either the same or unit length!');
    return;
end

N=max([nx,ny,nz]);
Omega=zeros(N,1);
for i=1:N,  
  xs = x( min(i,nx) );
  ys = y( min(i,ny) );
  zs = z( min(i,nz) );
  
  if (abs(zs)<1e-4)
      inside = ((xs^2/a^2)+(ys^2/b^2) < 1);
      if inside
          Omega = sign(zs)*pi*2;
      else
          Omega = 0;
      end
  else
    
  % initial number of integration points
  Nint = 4;
  maxiter=10;
  integral=zeros(maxiter,1);
  for iter=1:maxiter,
    theta=([0:Nint]/Nint - 0.5)*pi;  dtheta=theta(2)-theta(1);
    weight=ones(size(theta)); weight(1)=0.5; weight(end)=0.5;
    dx=a*sin(theta)-xs; 
    dy1=b*cos(theta)-ys; dy2=-b*cos(theta)-ys;
    r1=sqrt(dx.^2+dy1.^2+zs.^2); r2=sqrt(dx.^2+dy2.^2+zs.^2); 
    integrand = cos(theta)./(dx.^2+zs.^2).*( dy1./r1 - dy2./r2 ); 
    integral(iter) = sum( integrand.*weight )*dtheta;
    Omega(i) = a*zs*integral(iter);
    Nint = Nint*4;
    
    % print intermediate result
    %disp(sprintf('iter=%d Nint=%d Omega=%e',iter,Nint,Omega(i)));
    
    % test convergence
    converged = 0;
    if iter>1
        if abs(integral(iter)-integral(iter-1)) < tol*integral(iter)
            converged = 1;
        end
    end
    if converged
        break;
    end
    
    if iter==maxiter
        disp('error: solid_angle_ellipse: maxiter exceeded');
    end
  end
  
  % plot integrand
  %plot(theta, integrand, '.-');
  end
end