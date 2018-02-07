function P = fit_aEB(P)
% fit lattice constant a0, cohesive energy Ecoh, bulk modulus B
% from data: (a, E)

P.a = P.data(:,1);
P.Epot = P.data(:,3);
P.V = dot(cross(P.C1,P.C2),P.C3)*P.a.^3;

[P.P1,P.S1] = polyfit(P.a,P.Epot,2);
P.a0 = -P.P1(2)/(2*P.P1(1)); 
P.V0 = dot(cross(P.C1,P.C2),P.C3)*P.a0^3;
P.Ecoh = polyval(P.P1,P.a0)/P.N;
[P.P2,P.S2] = polyfit(P.V,P.Epot,2);
P.B = 2*P.V0*P.P2(1) * 160.2; % conversion from eV/A^3 to GPa

P.aaxis = [0:0.01:1]*(max(P.a)-min(P.a))+min(P.a);
P.Efit = polyval(P.P1,P.aaxis);
P.vaxis = dot(cross(P.C1,P.C2),P.C3)*P.aaxis.^3;