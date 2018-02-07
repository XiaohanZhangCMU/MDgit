% Find elastic constant C44 of Si, SiGe, Ge from VASP calculations

generic = struct('N', 0, 'C1', 0, 'C2', 0, 'C3', 0, 'data', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA Si ECUT = 500 data from ~/Codes/VASP/Si_bulk/LDA/ELatt.C11.dat
LDA_Si = generic;
LDA_Si.N = 16;
LDA_Si.C1 = [ 1 1 0];
LDA_Si.C2 = [-1 1 0];
LDA_Si.C3 = [ 0 0 1];
LDA_Si.data = [  % latt_const.(A)   Epot(eV) 
0.997 3  -.95640920E+02  -.95641413E+02 -.258146E-03
0.998 3  -.95641800E+02  -.95642292E+02 -.113930E-03
0.999 3  -.95642312E+02  -.95642804E+02 -.282626E-04
1.000 1  -.95642450E+02  -.95642942E+02 -.956425E+02
1.001 3  -.95642218E+02  -.95642709E+02 -.278599E-04
1.002 3  -.95641621E+02  -.95642112E+02 -.110732E-03
1.003 3  -.95640660E+02  -.95641150E+02 -.247371E-03
];
LDA_Si.a0 = 5.389282;   % (A) from bulkmodulus.m
LDA_Si.Ecoh = 5.974940; % (eV)
LDA_Si.B  = 95.695892;   % (GPa)
LDA_Si.C11=157.953597;  % (GPa)
LDA_Si.C12= 64.567039;  % (GPa)

LDA_Si = fit_C44(LDA_Si);
% Result : C44 = 77.470106 (GPa)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SiGe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA SiGe ECUT = 500 data from ~/Codes/VASP/SiGe_bulk/LDA/ELatt.C11.dat
LDA_SiGe = generic;
LDA_SiGe.N = 16;
LDA_SiGe.C1 = [ 1 1 0];
LDA_SiGe.C2 = [-1 1 0];
LDA_SiGe.C3 = [ 0 0 1];
LDA_SiGe.data = [  % latt_const.(A)   Epot(eV) 
0.997 3  -.88831724E+02  -.88832210E+02 -.248004E-03
0.998 3  -.88832639E+02  -.88833125E+02 -.109448E-03
0.999 3  -.88833192E+02  -.88833678E+02 -.271593E-04
1.000 1  -.88833383E+02  -.88833869E+02 -.888334E+02
1.001 3  -.88833220E+02  -.88833706E+02 -.267804E-04
1.002 3  -.88832704E+02  -.88833189E+02 -.106417E-03
1.003 3  -.88831833E+02  -.88832319E+02 -.237756E-03
];
LDA_SiGe.a0 = 5.510926;   % (A) from bulkmodulus.m
LDA_SiGe.Ecoh = 5.548209; % (eV)
LDA_SiGe.B  = 84.684935;   % (GPa)
LDA_SiGe.C11=139.427665;   % (GPa)
LDA_SiGe.C12= 57.313570;   % (GPa)
%
LDA_SiGe = fit_C44(LDA_SiGe);
% Result : C44 = 72.358764 (GPa)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA Ge ECUT = 500 data from ~/Codes/VASP/Ge_bulk/LDA/ELatt.C11.dat
LDA_Ge = generic;
LDA_Ge.N = 16;
LDA_Ge.C1 = [ 1 1 0];
LDA_Ge.C2 = [-1 1 0];
LDA_Ge.C3 = [ 0 0 1];
LDA_Ge.data = [  % latt_const.(A)   Epot(eV) 
0.997 3  -.82831720E+02  -.82832214E+02 -.239467E-03
0.998 3  -.82832609E+02  -.82833104E+02 -.105699E-03
0.999 3  -.82833131E+02  -.82833627E+02 -.262475E-04
1.000 1  -.82833319E+02  -.82833816E+02 -.828333E+02
1.001 3  -.82833179E+02  -.82833677E+02 -.258992E-04
1.002 3  -.82832708E+02  -.82833208E+02 -.102894E-03
1.003 3  -.82831907E+02  -.82832409E+02 -.229956E-03
];
LDA_Ge.a0 = 5.644045;   % (A) from bulkmodulus.m
LDA_Ge.Ecoh = 5.170127; % (eV)
LDA_Ge.B  = 72.767175;   % (GPa)
LDA_Ge.C11=118.930956;  % (GPa)
LDA_Ge.C12= 49.685284;  % (GPa)
%
LDA_Ge = fit_C44(LDA_Ge);
% Result : C44 = 64.949594 (GPa)

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

disp(sprintf('LDA_Si:   B = %.6f GPa  C11 = %.6f GPa C12 = %.6f GPa C44 = %.6f GPa', ...
              LDA_Si.B, LDA_Si.C11, LDA_Si.C12, LDA_Si.C44));
disp(sprintf('LDA_SiGe: B = %.6f GPa  C11 = %.6f GPa C12 = %.6f GPa C44 = %.6f GPa', ...
              LDA_SiGe.B, LDA_SiGe.C11, LDA_SiGe.C12, LDA_SiGe.C44));
disp(sprintf('LDA_Ge:   B = %.6f GPa  C11 = %.6f GPa C12 = %.6f GPa C44 = %.6f GPa', ...
              LDA_Ge.B, LDA_Ge.C11, LDA_Ge.C12, LDA_Ge.C44));

