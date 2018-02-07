% Find a0, B, Ecoh of Si, SiGe, Ge from VASP calculations
% a0 : eql. latt. const.
% B  : Bulk modulus (GPa)
% Ecoh: cohesive energy eV

generic = struct('N', 0, 'C1', 0, 'C2', 0, 'C3', 0, 'data', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA ECUT = 500 Si data from ~/Codes/VASP/Si/LDA/ELatt.dat
LDA_Si = generic;
LDA_Si.N = 2;
LDA_Si.C1 = [.5 .5 0];
LDA_Si.C2 = [0 .5 .5];
LDA_Si.C3 = [.5 0 .5];
LDA_Si.data = [  % latt_const.(A)   Epot(eV) 
5.382 1   -.11949688E+02   -.11949930E+02  -.119497E+02
5.383 1   -.11949738E+02   -.11949980E+02  -.119497E+02
5.384 1   -.11949779E+02   -.11950021E+02  -.119498E+02
5.385 1   -.11949813E+02   -.11950055E+02  -.119498E+02
5.386 1   -.11949842E+02   -.11950084E+02  -.119498E+02
5.387 1   -.11949861E+02   -.11950103E+02  -.119499E+02
5.388 1   -.11949873E+02   -.11950115E+02  -.119499E+02
5.389 1   -.11949878E+02   -.11950120E+02  -.119499E+02
5.390 1   -.11949880E+02   -.11950122E+02  -.119499E+02
5.391 1   -.11949870E+02   -.11950113E+02  -.119499E+02
5.392 1   -.11949854E+02   -.11950096E+02  -.119499E+02
5.393 1   -.11949831E+02   -.11950073E+02  -.119498E+02
5.394 1   -.11949800E+02   -.11950042E+02  -.119498E+02
5.395 1   -.11949761E+02   -.11950004E+02  -.119498E+02
5.396 1   -.11949716E+02   -.11949958E+02  -.119497E+02
%5.397 1   -.11949665E+02   -.11949907E+02  -.119497E+02
];
%
LDA_Si = fit_aEB(LDA_Si);
% Result : a0 = 5.389282 (A), Ecoh = 5.974940 (eV), B = 95.695892 (GPa)

%VASP GGA Si ECUT = 500 data from ~/Codes/VASP/Si/GGA/ELatt.dat
GGA_Si = generic;
GGA_Si.N = 2;
GGA_Si.C1 = [.5 .5 0];
GGA_Si.C2 = [0 .5 .5];
GGA_Si.C3 = [.5 0 .5];
GGA_Si.data = [  % latt_const.(A)   Epot(eV) 
5.462 1   -.10846262E+02   -.10846504E+02  -.108463E+02
5.463 1   -.10846309E+02   -.10846552E+02  -.108463E+02
5.464 1   -.10846341E+02   -.10846584E+02  -.108463E+02
5.465 1   -.10846367E+02   -.10846609E+02  -.108464E+02
5.466 1   -.10846390E+02   -.10846632E+02  -.108464E+02
5.467 1   -.10846401E+02   -.10846644E+02  -.108464E+02
5.468 1   -.10846406E+02   -.10846648E+02  -.108464E+02
5.469 1   -.10846410E+02   -.10846652E+02  -.108464E+02
5.470 1   -.10846401E+02   -.10846643E+02  -.108464E+02
5.471 1   -.10846385E+02   -.10846627E+02  -.108464E+02
5.472 1   -.10846362E+02   -.10846605E+02  -.108464E+02
5.473 1   -.10846344E+02   -.10846586E+02  -.108463E+02
5.474 1   -.10846308E+02   -.10846550E+02  -.108463E+02
5.475 1   -.10846265E+02   -.10846507E+02  -.108463E+02
5.476 1   -.10846221E+02   -.10846463E+02  -.108462E+02
];
%
GGA_Si = fit_aEB(GGA_Si);
% Result : a0 = 5.468507 (A), Ecoh = 5.423204 (eV), B = 87.536423 (GPa)

%MEAM Si data from ~/Codes/MD++/runs/si-meam/ELatt.dat
MEAM_Si = generic;
MEAM_Si.N = 512;
MEAM_Si.C1 = [4 0 0];
MEAM_Si.C2 = [0 4 0];
MEAM_Si.C3 = [0 0 4];
MEAM_Si.data = [  % latt_const.(A)   Epot(eV) 
5.426 1 -2.37053609971399e+03
5.427 1 -2.37054471188054e+03
5.428 1 -2.37055140422060e+03
5.429 1 -2.37055618016874e+03
5.430 1 -2.37055904315490e+03
5.431 1 -2.37055999660445e+03
5.432 1 -2.37055904393808e+03
5.433 1 -2.37055618857204e+03
5.434 1 -2.37055143391781e+03
5.435 1 -2.37054478338231e+03
5.436 1 -2.37053624036802e+03
];
%
MEAM_Si = fit_aEB(MEAM_Si);
% Result : a0 = 5.431 (A), Ecoh = 4.63 (eV), B = 97.5 (GPa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SiGe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA ECUT = 500 SiGe data from ~/Codes/VASP/SiGe/LDA/ELatt.dat
LDA_SiGe = generic;
LDA_SiGe.N = 2;
LDA_SiGe.C1 = [.5 .5 0];
LDA_SiGe.C2 = [0 .5 .5];
LDA_SiGe.C3 = [.5 0 .5];
LDA_SiGe.data = [  % latt_const.(A)   Epot(eV) 
5.505 1   -.11096302E+02   -.11096544E+02  -.110963E+02
5.506 1   -.11096337E+02   -.11096579E+02  -.110963E+02
5.507 1   -.11096368E+02   -.11096609E+02  -.110964E+02
5.508 1   -.11096389E+02   -.11096631E+02  -.110964E+02
5.509 1   -.11096404E+02   -.11096646E+02  -.110964E+02
5.510 1   -.11096416E+02   -.11096658E+02  -.110964E+02
5.511 1   -.11096418E+02   -.11096660E+02  -.110964E+02
5.512 1   -.11096413E+02   -.11096655E+02  -.110964E+02
5.513 1   -.11096402E+02   -.11096644E+02  -.110964E+02
5.514 1   -.11096387E+02   -.11096629E+02  -.110964E+02
5.515 1   -.11096363E+02   -.11096605E+02  -.110964E+02
5.516 1   -.11096332E+02   -.11096574E+02  -.110963E+02
5.517 1   -.11096297E+02   -.11096539E+02  -.110963E+02
];
%
LDA_SiGe = fit_aEB(LDA_SiGe);
% Result : a0 = 5.510926 (A), Ecoh = 5.548209 (eV), B = 84.684935 (GPa)

%VASP GGA SiGe ECUT = 500 data from ~/Codes/VASP/SiGe/GGA/ELatt.dat
GGA_SiGe = generic;
GGA_SiGe.N = 2;
GGA_SiGe.C1 = [.5 .5 0];
GGA_SiGe.C2 = [0 .5 .5];
GGA_SiGe.C3 = [.5 0 .5];
GGA_SiGe.data = [  % latt_const.(A)   Epot(eV) 
5.607 1   -.98407838E+01   -.98410255E+01  -.984078E+01
5.608 1   -.98408285E+01   -.98410703E+01  -.984083E+01
5.609 1   -.98408718E+01   -.98411136E+01  -.984087E+01
5.610 1   -.98409048E+01   -.98411465E+01  -.984090E+01
5.612 1   -.98409574E+01   -.98411991E+01  -.984096E+01
5.613 1   -.98409726E+01   -.98412144E+01  -.984097E+01
5.614 1   -.98409820E+01   -.98412237E+01  -.984098E+01
5.615 1   -.98409901E+01   -.98412318E+01  -.984099E+01
5.616 1   -.98409878E+01   -.98412295E+01  -.984099E+01
5.617 1   -.98409796E+01   -.98412213E+01  -.984098E+01
5.618 1   -.98409656E+01   -.98412073E+01  -.984097E+01
5.619 1   -.98409478E+01   -.98411895E+01  -.984095E+01
5.620 1   -.98409222E+01   -.98411639E+01  -.984092E+01
5.622 1   -.98408587E+01   -.98411005E+01  -.984086E+01
5.623 1   -.98408158E+01   -.98410575E+01  -.984082E+01
5.624 1   -.98407670E+01   -.98410088E+01  -.984077E+01
];
%
GGA_SiGe = fit_aEB(GGA_SiGe);
% Result : a0 = 5.615328 (A), Ecoh = 4.920494 (eV), B = 74.557008 (GPa)

%MEAM GGA SiGe data from ~/Codes/MD++/sige-meam/ELatt.dat
MEAM_SiGe = generic;
MEAM_SiGe.N = 512;
MEAM_SiGe.C1 = [4 0 0];
MEAM_SiGe.C2 = [0 4 0];
MEAM_SiGe.C3 = [0 0 4];
MEAM_SiGe.data = [  % latt_const.(A)   Epot(eV) 
5.533 1 -2.16061883941950e+03
5.534 1 -2.16062654102315e+03
5.535 1 -2.16063250420578e+03
5.536 1 -2.16063673206754e+03
5.537 1 -2.16063922770444e+03
5.538 1 -2.16063999420842e+03
5.539 1 -2.16063903466717e+03
5.540 1 -2.16063635216429e+03
5.541 1 -2.16063194977931e+03
5.542 1 -2.16062583058753e+03
5.543 1 -2.16061799766019e+03
];
%
MEAM_SiGe = fit_aEB(MEAM_SiGe);
% Result : a0 = 5.615328 (A), Ecoh = 4.920494 (eV), B = 74.557008 (GPa)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VASP LDA Ge ECUT = 500 data from ~/Codes/VASP/Ge_bulk/LDA/ELatt.dat
LDA_Ge = generic;
LDA_Ge.N = 2;
LDA_Ge.C1 = [.5 .5 0];
LDA_Ge.C2 = [0 .5 .5];
LDA_Ge.C3 = [.5 0 .5];
LDA_Ge.data = [  % latt_const.(A)   Epot(eV) 
5.635 1   -.10340018E+02   -.10340338E+02  -.103400E+02
5.636 1   -.10340068E+02   -.10340390E+02  -.103401E+02
5.637 1   -.10340111E+02   -.10340437E+02  -.103401E+02
5.638 1   -.10340149E+02   -.10340477E+02  -.103401E+02
5.639 1   -.10340181E+02   -.10340512E+02  -.103402E+02
5.640 1   -.10340207E+02   -.10340541E+02  -.103402E+02
5.641 1   -.10340227E+02   -.10340564E+02  -.103402E+02
5.642 1   -.10340243E+02   -.10340583E+02  -.103402E+02
5.643 1   -.10340251E+02   -.10340594E+02  -.103403E+02
5.644 1   -.10340254E+02   -.10340600E+02  -.103403E+02
5.645 1   -.10340254E+02   -.10340603E+02  -.103403E+02
5.646 1   -.10340244E+02   -.10340597E+02  -.103402E+02
5.647 1   -.10340230E+02   -.10340586E+02  -.103402E+02
5.648 1   -.10340209E+02   -.10340568E+02  -.103402E+02
5.649 1   -.10340183E+02   -.10340546E+02  -.103402E+02
5.650 1   -.10340151E+02   -.10340518E+02  -.103402E+02
5.651 1   -.10340113E+02   -.10340483E+02  -.103401E+02
5.652 1   -.10340073E+02   -.10340447E+02  -.103401E+02
5.653 1   -.10340024E+02   -.10340402E+02  -.103400E+02
%5.654 1   -.10339968E+02   -.10340350E+02  -.103400E+02
%5.655 1   -.10339909E+02   -.10340295E+02  -.103399E+02
%5.656 1   -.10339843E+02   -.10340233E+02  -.103398E+02
];
%
LDA_Ge = fit_aEB(LDA_Ge);
% Result : a0 = 5.644045 (A), Ecoh = 5.170127 (eV), B = 72.767175 (GPa)

%VASP GGA Ge ECUT = 500 data from ~/Codes/VASP/Ge_bulk/GGA/ELatt.dat
GGA_Ge = generic;
GGA_Ge.N = 2;
GGA_Ge.C1 = [.5 .5 0];
GGA_Ge.C2 = [0 .5 .5];
GGA_Ge.C3 = [.5 0 .5];
GGA_Ge.data = [  % latt_const.(A)   Epot(eV) 
5.770 1   -.89752731E+01   -.89759038E+01  -.897527E+01
5.771 1   -.89753195E+01   -.89759514E+01  -.897532E+01
5.772 1   -.89753579E+01   -.89759909E+01  -.897536E+01
5.773 1   -.89753915E+01   -.89760252E+01  -.897539E+01
5.774 1   -.89754241E+01   -.89760584E+01  -.897542E+01
5.775 1   -.89754477E+01   -.89760825E+01  -.897545E+01
5.776 1   -.89754665E+01   -.89761014E+01  -.897547E+01
5.777 1   -.89754818E+01   -.89761168E+01  -.897548E+01
5.778 1   -.89754908E+01   -.89761255E+01  -.897549E+01
5.779 1   -.89754949E+01   -.89761293E+01  -.897549E+01
5.780 1   -.89754957E+01   -.89761295E+01  -.897550E+01
5.781 1   -.89754902E+01   -.89761231E+01  -.897549E+01
5.782 1   -.89754798E+01   -.89761117E+01  -.897548E+01
5.783 1   -.89754645E+01   -.89760953E+01  -.897546E+01
5.784 1   -.89754481E+01   -.89760774E+01  -.897545E+01
5.785 1   -.89754232E+01   -.89760510E+01  -.897542E+01
5.786 1   -.89753936E+01   -.89760197E+01  -.897539E+01
5.787 1   -.89753602E+01   -.89759843E+01  -.897536E+01
5.788 1   -.89753211E+01   -.89759431E+01  -.897532E+01
5.789 1   -.89752773E+01   -.89758969E+01  -.897528E+01
%5.790 1   -.89752309E+01   -.89758481E+01  -.897523E+01
%5.791 1   -.89751776E+01   -.89757921E+01  -.897518E+01
];
%
GGA_Ge = fit_aEB(GGA_Ge);
% Result : a0 = 5.779521 (A), Ecoh = 4.487748 (eV), B = 60.081859 (GPa)

%MEAM Ge data from ~/Codes/MD++/runs/ge-meam/ELatt.dat
MEAM_Ge = generic;
MEAM_Ge.N = 512;
MEAM_Ge.C1 = [4 0 0];
MEAM_Ge.C2 = [0 4 0];
MEAM_Ge.C3 = [0 0 4];
MEAM_Ge.data = [  % latt_const.(A)   Epot(eV) 
5.653 1 -1.97118449276057e+03
5.654 1 -1.97119062369550e+03
5.655 1 -1.97119521784151e+03
5.656 1 -1.97119827789813e+03
5.657 1 -1.97119980656143e+03
5.658 1 -1.97119980652383e+03
5.659 1 -1.97119828047419e+03
5.660 1 -1.97119523109798e+03
5.661 1 -1.97119066107690e+03
5.662 1 -1.97118457308924e+03
];
%
MEAM_Ge = fit_aEB(MEAM_Ge);
% Result : a0 = 5.6575 (A), Ecoh = 3.85 (eV), B = 74.99 (GPa)


fs = 16;
figure(1);
subplot(2,1,1);
plot(LDA_SiGe.a,LDA_SiGe.Epot/LDA_SiGe.N - LDA_SiGe.Ecoh,'b.', ...
     LDA_SiGe.aaxis,LDA_SiGe.Efit/LDA_SiGe.N - LDA_SiGe.Ecoh,'r-', ...
     LDA_Si.a,LDA_Si.Epot/LDA_Si.N - LDA_Si.Ecoh,'b.', ...
     LDA_Si.aaxis,LDA_Si.Efit/LDA_Si.N - LDA_Si.Ecoh,'r-', ...
     LDA_Ge.a,LDA_Ge.Epot/LDA_Ge.N - LDA_Ge.Ecoh,'b.', ...
     LDA_Ge.aaxis,LDA_Ge.Efit/LDA_Ge.N - LDA_Ge.Ecoh,'r-', ...
     MEAM_SiGe.a,MEAM_SiGe.Epot/MEAM_SiGe.N - MEAM_SiGe.Ecoh,'m.', ...
     MEAM_SiGe.aaxis,MEAM_SiGe.Efit/MEAM_SiGe.N - MEAM_SiGe.Ecoh,'r-', ...
     MEAM_Si.a,MEAM_Si.Epot/MEAM_Si.N - MEAM_Si.Ecoh,'m.', ...
     MEAM_Si.aaxis,MEAM_Si.Efit/MEAM_Si.N - MEAM_Si.Ecoh,'r-', ...
     MEAM_Ge.a,MEAM_Ge.Epot/MEAM_Ge.N - MEAM_Ge.Ecoh,'m.', ...
     MEAM_Ge.aaxis,MEAM_Ge.Efit/MEAM_Ge.N - MEAM_Ge.Ecoh,'r-'     );
set(gca,'FontSize',fs); 
%xlabel('a  (angstrom)');
text('Interpreter','latex','String','$a \ \ (\AA)$',...
     'Position',[5.52,-9e-5],'FontSize',fs);
ylabel('E - E_{coh} (eV)');
%title('LDA + MEAM');

t1=text(LDA_Si.a0-0.005,10e-5,'Si (LDA)');    set(t1,'FontSize',fs-2);
t2=text(LDA_SiGe.a0-0.010,8e-5,'SiGe (LDA)'); set(t2,'FontSize',fs-2);
t3=text(LDA_Ge.a0-0.005,13e-5,'Ge (LDA)');    set(t3,'FontSize',fs-2);
t1a=text(MEAM_Si.a0-0.005,6e-5,'Si (MEAM)');  set(t1a,'FontSize',fs-2,'Color','blue');
t2a=text(MEAM_SiGe.a0-0.007,5.4e-5,'SiGe (MEAM)');set(t2a,'FontSize',fs-2,'Color','blue');
t3a=text(MEAM_Ge.a0-0.005,4e-5,'Ge (MEAM)');  set(t3a,'FontSize',fs-2,'Color','blue');
xlim([5.3 5.8]);
ylim([-5 15]*1e-5);

%figure(2);
subplot(2,1,2);
plot(GGA_SiGe.a,GGA_SiGe.Epot/GGA_SiGe.N - GGA_SiGe.Ecoh,'b.', ...
     GGA_SiGe.aaxis,GGA_SiGe.Efit/GGA_SiGe.N - GGA_SiGe.Ecoh,'r-', ...
     GGA_Si.a,GGA_Si.Epot/GGA_Si.N - GGA_Si.Ecoh,'b.', ...
     GGA_Si.aaxis,GGA_Si.Efit/GGA_Si.N - GGA_Si.Ecoh,'r-', ...
     GGA_Ge.a,GGA_Ge.Epot/GGA_Ge.N - GGA_Ge.Ecoh,'b.', ...
     GGA_Ge.aaxis,GGA_Ge.Efit/GGA_Ge.N - GGA_Ge.Ecoh,'r-' );
set(gca,'FontSize',fs);
%xlabel('a  (angstrom)');
text('Interpreter','latex','String','$a \ \ (\AA)$',...
     'Position',[5.52,-9e-5],'FontSize',fs);
ylabel('E - E_{coh} (eV)');
%title('GGA');
t4=text(GGA_Si.a0-0.008,11e-5,'Si (GGA)');     set(t4,'FontSize',fs-2);
t5=text(GGA_SiGe.a0-0.020,12e-5,'SiGe (GGA)'); set(t5,'FontSize',fs-2);
t6=text(GGA_Ge.a0-0.050,13e-5,'Ge (GGA)');     set(t6,'FontSize',fs-2);
xlim([5.3 5.8]);
ylim([-5 15]*1e-5);

disp(sprintf('LDA_Si:   a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              LDA_Si.a0, LDA_Si.Ecoh, LDA_Si.B));

disp(sprintf('GGA_Si:   a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              GGA_Si.a0, GGA_Si.Ecoh, GGA_Si.B));

disp(sprintf('MEAM_Si:  a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              MEAM_Si.a0, MEAM_Si.Ecoh, MEAM_Si.B));

disp(' ');

disp(sprintf('LDA_SiGe: a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              LDA_SiGe.a0, LDA_SiGe.Ecoh, LDA_SiGe.B));

disp(sprintf('GGA_SiGe: a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              GGA_SiGe.a0, GGA_SiGe.Ecoh, GGA_SiGe.B));

disp(sprintf('MEAM_SiGe:  a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              MEAM_SiGe.a0, MEAM_SiGe.Ecoh, MEAM_SiGe.B));
              
disp(' ');

disp(sprintf('LDA_Ge:   a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              LDA_Ge.a0, LDA_Ge.Ecoh, LDA_Ge.B));
          
disp(sprintf('GGA_Ge:   a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              GGA_Ge.a0, GGA_Ge.Ecoh, GGA_Ge.B));
              
disp(sprintf('MEAM_Ge:  a0 = %.6f A, Ecoh = %.6f eV, B = %.6f GPa', ...
              MEAM_Ge.a0, MEAM_Ge.Ecoh, MEAM_Ge.B));
              
              
figure(2);
title('experimental lattice constant of Si(1-x)Ge(x)');
x = [0:0.05:1];
a0x = inline('5.431 + 0.20*x + 0.027*x.^2');
plot(x,a0x(x));
disp(sprintf('experiment: Si: a0 = %.3f  SiGe: a0 = %.3f  Ge: a0 = %.3f', ...
              a0x(0),a0x(0.5),a0x(1)));
              

