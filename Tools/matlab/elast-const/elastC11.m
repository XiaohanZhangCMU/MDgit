% Find elastic constant C11 of Si, SiGe, Ge from VASP calculations

generic = struct('N', 0, 'C1', 0, 'C2', 0, 'C3', 0, 'data', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA Si ECUT = 500 data from ~/Codes/VASP/Si_bulk/LDA/ELatt.C11.dat
LDA_Si = generic;
LDA_Si.N = 8;
LDA_Si.C1 = [1 0 0];
LDA_Si.C2 = [0 1 0];
LDA_Si.C3 = [0 0 1];
LDA_Si.data = [  % latt_const.(A)   Epot(eV) 
0.997 1  -.47799249E+02  -.47800224E+02 -.477992E+02
0.998 1  -.47799632E+02  -.47800604E+02 -.477996E+02
0.999 1  -.47799859E+02  -.47800828E+02 -.477999E+02
1.000 1  -.47799948E+02  -.47800917E+02 -.477999E+02
1.001 1  -.47799866E+02  -.47800836E+02 -.477999E+02
1.002 1  -.47799634E+02  -.47800606E+02 -.477996E+02
1.003 1  -.47799248E+02  -.47800223E+02 -.477992E+02
];
LDA_Si.a0 = 5.389282;   % (A) from bulkmodulus.m
LDA_Si.Ecoh = 5.974940; % (eV)
LDA_Si.B = 95.695892;   % (GPa)
LDA_Si = fit_C11(LDA_Si);
% Result : C11 = 157.953597 (GPa), C12 = 64.567039 (GPa)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SiGe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA SiGe ECUT = 500 data from ~/Codes/VASP/SiGe_bulk/LDA/ELatt.C11.dat
LDA_SiGe = generic;
LDA_SiGe.N = 8;
LDA_SiGe.C1 = [1 0 0];
LDA_SiGe.C2 = [0 1 0];
LDA_SiGe.C3 = [0 0 1];
LDA_SiGe.data = [  % latt_const.(A)   Epot(eV) 
0.997 1  -.44385012E+02  -.44385985E+02 -.443850E+02
0.998 1  -.44385372E+02  -.44386341E+02 -.443854E+02
0.999 1  -.44385586E+02  -.44386553E+02 -.443856E+02
1.000 1  -.44385657E+02  -.44386624E+02 -.443857E+02
1.001 1  -.44385579E+02  -.44386547E+02 -.443856E+02
1.002 1  -.44385356E+02  -.44386326E+02 -.443854E+02
1.003 1  -.44384989E+02  -.44385963E+02 -.443850E+02
];
LDA_SiGe.a0 = 5.510926;   % (A) from bulkmodulus.m
LDA_SiGe.Ecoh = 5.548209; % (eV)
LDA_SiGe.B = 84.684935;   % (GPa)
LDA_SiGe = fit_C11(LDA_SiGe);
% Result : C11 = 139.427665 (GPa), C12 = 57.313570 (GPa)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA Ge ECUT = 500 data from ~/Codes/VASP/Ge_bulk/LDA/ELatt.C11.dat
LDA_Ge = generic;
LDA_Ge.N = 8;
LDA_Ge.C1 = [1 0 0];
LDA_Ge.C2 = [0 1 0];
LDA_Ge.C3 = [0 0 1];
LDA_Ge.data = [  % latt_const.(A)   Epot(eV) 
0.997 1  -.41360230E+02  -.41361558E+02 -.413602E+02
0.998 1  -.41360551E+02  -.41361895E+02 -.413606E+02
0.999 1  -.41360738E+02  -.41362100E+02 -.413607E+02
1.000 1  -.41360789E+02  -.41362174E+02 -.413608E+02
1.001 1  -.41360712E+02  -.41362124E+02 -.413607E+02
1.002 1  -.41360498E+02  -.41361939E+02 -.413605E+02
1.003 1  -.41360150E+02  -.41361625E+02 -.413601E+02
];
LDA_Ge.a0 = 5.644045;   % (A) from bulkmodulus.m
LDA_Ge.Ecoh = 5.170127; % (eV)
LDA_Ge.B = 72.767175;   % (GPa)
%
LDA_Ge = fit_C11(LDA_Ge);
% Result : C11 = 118.930956 (GPa), C12 = 49.685284 (GPa)

fs = 16;
figure(1);
plot(LDA_Si.eps,  (LDA_Si.Epot-min(LDA_Si.Epot))/LDA_Si.V0,'b.', ...
     LDA_Si.eaxis,(LDA_Si.Efit-min(LDA_Si.Epot))/LDA_Si.V0,'g-', ...
     LDA_SiGe.eps,  (LDA_SiGe.Epot-min(LDA_SiGe.Epot))/LDA_SiGe.V0,'r.', ...
     LDA_SiGe.eaxis,(LDA_SiGe.Efit-min(LDA_SiGe.Epot))/LDA_SiGe.V0,'b-', ...
     LDA_Ge.eps,  (LDA_Ge.Epot-min(LDA_Ge.Epot))/LDA_Ge.V0,'b.', ...
     LDA_Ge.eaxis,(LDA_Ge.Efit-min(LDA_Ge.Epot))/LDA_Ge.V0,'r-' );
set(gca,'FontSize',fs); 
xlabel('\epsilon');
ylabel('E - E_{coh} (eV)');

disp(sprintf('LDA_Si:   B = %.6f GPa  C11 = %.6f GPa C12 = %.6f GPa', ...
              LDA_Si.B, LDA_Si.C11, LDA_Si.C12));
disp(sprintf('LDA_SiGe: B = %.6f GPa  C11 = %.6f GPa C12 = %.6f GPa', ...
              LDA_SiGe.B, LDA_SiGe.C11, LDA_SiGe.C12));
disp(sprintf('LDA_Ge:   B = %.6f GPa  C11 = %.6f GPa C12 = %.6f GPa', ...
              LDA_Ge.B, LDA_Ge.C11, LDA_Ge.C12));

