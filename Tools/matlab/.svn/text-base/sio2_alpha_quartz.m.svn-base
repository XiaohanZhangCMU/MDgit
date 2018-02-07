% create lattice.h files for SiO2 structure in MD++
% alpha-quartz (low quartz)   **** this file ****
% beta-quartz  (high quartz) 
% data from http://cst-www.nrl.navy.mil/lattice/struk.xmol/sio2a.pos

% alpha-quartz
a = 2.4567*2;
c = 5.4052;

A1 = [1/2 -sqrt(3)/2 0]*a;
A2 = [1/2  sqrt(3)/2 0]*a;
A3 = [0    0         1]*c;

% scaled coordinates                  Cartesian coordinates
coords = [
.46990000  0.00000000   .66666667     1.15440333 -1.99948522  3.60346667  %Si
0.00000000   .46990000   .33333333     1.15440333  1.99948522  1.80173333 %Si
-.46990000  -.46990000   .00000000    -2.30880666  -.00000000   .00000000 %Si
.41410000   .26810000   .78540000     1.67596074  -.62124887  4.24524408  %O
-.26810000   .14600000   .45206667     -.29996307  1.76204901  2.44351075 %O
-.14600000  -.41410000   .11873333    -1.37599767 -1.14080014   .64177741 %O
.26810000   .41410000  -.78540000     1.67596074   .62124887 -4.24524408  %O
-.41410000  -.14600000  -.11873333    -1.37599767  1.14080014  -.64177741 %O
.14600000  -.26810000  -.45206667     -.29996307 -1.76204901 -2.44351075  %O
];

H = [A1',A2',A3'];
s_Si = coords(1:3,1:3); %r_Si = coords(1:3,4:6);
r_Si = s_Si*H';

s_O  = coords(4:9,1:3); %r_O  = coords(4:9,4:6);
r_O  = s_O*H';

% Find supercell twice as large as primitive cell
H2 = [A2'+A1',A2'-A1',A3'];

r_Si2 = [r_Si; r_Si + ones(3,1)*A2];
r_O2 = [r_O ; r_O + ones(6,1)*A2];

s_Si2 = r_Si2*inv(H2)';  s_Si2 = s_Si2 - round(s_Si2);
s_O2  = r_O2 *inv(H2)';  s_O2  = s_O2  - round(s_O2);

r_Si2 = s_Si2*H';
r_O2  = s_O2 *H';


% plot primitive cell
figure(1);
p1=plot3(r_Si(:,1),r_Si(:,2),r_Si(:,3),'bo', r_O(:,1),r_O(:,2),r_O(:,3),'ro');
set(p1(1),'MarkerSize',9);
axis equal
grid
hold on
plot3([0 A1(1)],[0 A1(2)],[0 A1(3)],'k-');
plot3([0 A2(1)],[0 A2(2)],[0 A2(3)],'k-');
plot3([0 A3(1)],[0 A3(2)],[0 A3(3)],'k-');
hold off

% plot unit cell
figure(2);
p2=plot3(r_Si2(:,1),r_Si2(:,2),r_Si2(:,3),'bo', r_O2(:,1),r_O2(:,2),r_O2(:,3),'ro');
set(p2(1),'MarkerSize',9);
axis equal
grid


% print out for MD++ source code
disp(' ');
disp('/* This goes to lattice.h */');
disp('    const double alpha_quartz_basis[54] = {');
for i=1:length(s_Si2(:,1)),
    disp(sprintf('\t%14.10f, %14.10f, %14.10f, /* Si */',s_Si2(i,1),s_Si2(i,2),s_Si2(i,3)));
end
for i=1:length(s_O2(:,1)),
    disp(sprintf('\t%14.10f, %14.10f, %14.10f, /* O  */',s_O2(i,1),s_O2(i,2),s_O2(i,3)));
end
disp('    };');
disp(' ');
disp('    const int alpha_quartz_species[18]={0,0,0,0,0,0,');
disp('                                        1,1,1,1,1,1, 1,1,1,1,1,1};');
