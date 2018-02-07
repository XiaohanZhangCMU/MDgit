% test solid_angle_ellipse.m (SolidAngleEllipse.C)

x = [-5:0.2:5];
y = [-5:0.2:5];
z = 1e-1;

r1=[-2 -2 0];
r2=[ 2 -2 0];
r3=[ 1  2 0];
r4=[-1  2 0];
r5=[-4  0.5 0];

Omega = zeros(length(y),length(x));

for j=1:length(y),
    for i=1:length(x),
        zt=x(i)*0.1+y(j)*0.05;
        %Omega(j,i) = solid_angle_polygon([r1;r2;r3],x(i),y(j),z);
        %Omega(j,i) = solid_angle_polygon([r1;r2;r3;r4],x(i),y(j),z);
        Omega(j,i) = solid_angle_polygon([r1;r2;r3;r4;r5],x(i),y(j),z);
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
