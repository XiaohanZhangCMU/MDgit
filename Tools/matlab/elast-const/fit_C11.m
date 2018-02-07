function P = fit_C11(P)
% fit elastic constant C11 from data: (a, E)

P.eps = P.data(:,1)-1;
P.Epot = P.data(:,3);
P.V0 = dot(cross(P.C1,P.C2),P.C3)*P.a0^3;

[P.P1,P.S1] = polyfit(P.eps,P.Epot,2);

P.C11 = 2*P.P1(1)/P.V0 *160.2; % conversion from eV/A^3 to GPa
P.C12 = (3*P.B - P.C11)/2;

P.eaxis = [0:0.01:1]*(max(P.eps)-min(P.eps))+min(P.eps);
P.Efit = polyval(P.P1,P.eaxis);
