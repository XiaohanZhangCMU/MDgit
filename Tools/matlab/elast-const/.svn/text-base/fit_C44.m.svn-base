function P = fit_C44(P)
% fit elastic constant C11 from data: (a, E)

P.eps = P.data(:,1)-1;
P.Epot = P.data(:,3);
P.V0 = dot(cross(P.C1,P.C2),P.C3)*P.a0^3;

[P.P1,P.S1] = polyfit(P.eps,P.Epot,2);

P.C11p = 2*P.P1(1)/P.V0 *160.2; % conversion from eV/A^3 to GPa

%rotation from cubic axis to [110],[-1 1 0],[0 0 1]
A=0.5; B=0.25; C=0.25;

%How to obtain coefficients A, B, C  change to if(1) to verify
if(0)
% cubic coordinate system
  e1 = [1 0 0]; e2 = [0 1 0]; e3 = [0 0 1];
% new coordinate system
  e1p= [1 1 0]; e2p=[-1 1 0]; e3p= [0 0 1];
% normalization
  e1p=e1p/norm(e1p); e2p=e2p/norm(e2p); e3p=e3p/norm(e3p); 
% rotation matrix
  T = [ dot(e1,e1p) dot(e1,e2p) dot(e1,e3p)
        dot(e2,e1p) dot(e2,e2p) dot(e2,e3p)
        dot(e3,e1p) dot(e3,e2p) dot(e3,e3p) ];
  Cp = rotate_tensor_4th (cubic_elast_stiff(1,0,0), T); A = Cp(1,1,1,1);
  Cp = rotate_tensor_4th (cubic_elast_stiff(0,1,0), T); B = Cp(1,1,1,1)/2;
  Cp = rotate_tensor_4th (cubic_elast_stiff(0,0,1), T); C = Cp(1,1,1,1)/4;
end

%C11p = A*C11 + 2*B*C12 + 4*C*C14;
P.C44 = (P.C11p - A*P.C11 - 2*B*P.C12)/(4*C);

P.eaxis = [0:0.01:1]*(max(P.eps)-min(P.eps))+min(P.eps);
P.Efit = polyval(P.P1,P.eaxis);
