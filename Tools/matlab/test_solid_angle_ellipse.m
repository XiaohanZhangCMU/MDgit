% test solid_angle_ellipse.m (SolidAngleEllipse.C)

x = [-5:0.2:5];
y = [-5:0.2:5];
z = 1e-1;

a = 2;
b = 4;

Omega = zeros(length(y),length(x));

for j=1:length(y),
    for i=1:length(x),
        zt=x(i)*0.1+y(j)*0.05;
        Omega(j,i) = SolidAngleEllipse(a,b,x(i),y(j),zt);
    end
end

figure(1);
mesh(x,y,Omega);
xlabel('x');
ylabel('y');

figure(2);
contour(x,y,Omega,50);
xlabel('x');
ylabel('y');
