% test solid_angle_ellipse.m (SolidAngleEllipse.C)

x = [-5:0.2:5];
y = [-5:0.2:5];
z = 1e-1;

r1=[-2 -2 0];
r2=[ 2 -2 0];
r3=[ 0  2 0];

Omega = zeros(length(y),length(x));

for j=1:length(y),
    for i=1:length(x),
        zt=x(i)*0.1+y(j)*0.05;
        %Omega(j,i) = solid_angle_triangle(r1,r2,r3,x(i),y(j),z);
        Omega(j,i) = SolidAngleTriangle(r1(1),r1(2),r1(3),r2(1),r2(2),r2(3),r3(1),r3(2),r3(3),x(i),y(j),z);
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
