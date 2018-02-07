%test cubic spline

load FS2Fe.out

%
%f0=1.0;p0=1.0;f1=3.0;p1=-4.0;dr=3.0;
f0=0;p0=0;f1=100;p1=-100;dr=3.0;

if(0)
    %read from table
eV = 1.60219e-19;
angstrom = 1.0e-10;

data=FS2Fe;
%ind=4000;
ind=4999;
%ind=5000;
f0=data(ind,2)/(eV*eV); 
p0=data(ind,4)*angstrom/(eV*eV);
if(ind<5000)
f1=data(ind+1,2)/(eV*eV);
p1=data(ind+1,4)*angstrom/(eV*eV);
else
    f1=0;
    p1=0;
end
dr = 6.564939997019770e-4;
end

qq=[0:0.01:1]*dr;

a=f0;
b=p0;
A1=(f1-a-b*dr);
A2=(p1-b)*dr;

d=(A2-2*A1)/dr^3;
c=(3*A1-A2)/dr^2;

f=a+b*qq+c*qq.^2+d*qq.^3;

plot(qq,f,'b.-');