%plot polarization degeneracy removal data

data = [
    -8  -0.1013 
    -4  -0.0180
    -2  -0.0062
    -1  -0.0020 
    1    0.0019
    2    0.0070
    4    0.0239
    8    0.1182
] ;

stress = data(:,1);
dE = data(:,2);

A=1.8e-3;
x=[-10:0.1:10];
y=sign(x).*(A*x.^2);

figure(1);
plot(stress,dE, 'o-', x,y);

set(gca,'FontSize',12);
xlabel('\sigma (GPa)');
ylabel('dE (eV)');


%figure(2);
%loglog(abs(stress),abs(dE), 'o-');
%xlim([1e-1, 1e1]);
