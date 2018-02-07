%---------------------------------------------------%
% Convert modified FS potential from MD++ to LAMMPS %
% and then convert to MD++ EAM format               %
%---------------------------------------------------%

%elename{1} = 'Mo';
%elename{1} = 'Ta';
elename{1} = 'Nb';

if (strcmp(elename{1},'Mo'))
 %Mo         d        A        beta  c    c0          c1         c2        B      alpha  b0
 fs_data = [ 4.114825 1.887117 0     3.25 43.4475218 -31.9332978 6.0804249 1223.0 3.90   2.7255];
 %            Z      atommass   lattconst
 mat_data = [ 42     95.94       3.1472 ];  latttype = 'BCC';
 %
 Nele = 1; Nrho = 5000; Nr = 5000;
 potfilename = 'mo.eam.fs';
 outpotfile = 'eamdata.Mo.fs.mdpp';
 %
elseif (strcmp(elename{1},'Ta'))
 %Ta         d        A        beta  c    c0          c1         c2        B      alpha  b0
 fs_data = [ 4.076980 2.591061 0     4.20 1.2157373   0.027147  -0.121735  91.0   1.05   2.8629];
 %            Z      atommass   lattconst
 mat_data = [ 73     180.9479   3.3058 ];  latttype = 'BCC';
 %
 Nele = 1; Nrho = 5000; Nr = 5000;
 potfilename = 'ta.eam.fs';
 outpotfile = 'eamdata.Ta.fs.mdpp';
 %
elseif (strcmp(elename{1},'Nb'))
 %Nb         d        A        beta  c    c0          c1         c2        B      alpha  b0
 fs_data = [ 3.915354 3.013789 0     4.20 -1.5640104  2.0055779 -0.4663764 48.0   0.8   2.8585];
 %            Z      atommass   lattconst
 mat_data = [ 41     92.9064   3.3008 ];  latttype = 'BCC';
 %
 Nele = 1; Nrho = 5000; Nr = 5000;
 potfilename = 'nb.eam.fs';
 outpotfile = 'eamdata.Nb.fs.mdpp';
 % 
else
 disp('unknown elename! do nothing.');
 return;
end


d = fs_data(1);                     % cut-off distance for rho(r)
A = fs_data(2);
beta = fs_data(3);
c = fs_data(4);                     % cut-off distance for phi(r)
c0= fs_data(5);
c1= fs_data(6);
c2= fs_data(7);
B = fs_data(8);
alpha= fs_data(9);
b0= fs_data(10);                    % nearest neighbor distance
Z = mat_data(1);                    % atomic number
atommass = mat_data(2);             % atomic mass
latticeconst = mat_data(3);         % lattice constant


rcut = max(c,d);
rhob0 = (b0-d).^2 + beta*(b0-d).^3/d;

rho_max = round(rhob0*8*1.3* 5);
rmax    = round(rcut * 1.5);

drho = rho_max/Nrho;  %drho = 2.0e-3;
dr   = rmax/Nr;

rho = [0:Nrho-1]*drho;
r   = [0:Nr-1  ]*dr;

%analytic expression of FS potential (Ackland)
Frho = -A*sqrt(rho);
rhor = (r<d).*((r-d).^2 + beta*(r-d).^3/d);
phir = (r<c).*((r-c).^2.*(c0+c1.*r+c2.*r.^2)) + (r<b0).*(B*(b0-r).^3.*exp(-alpha*r));

griddata = struct('Nele',Nele,'Nrho',Nrho,'drho',drho,'Nr',Nr,'dr',dr,'rcut',rcut);

elemdata{1} = struct('Z',Z,'atommass',atommass,'latticeconst',latticeconst,'latttype',latttype);
Frhodata{1} = Frho;
rhordata_iele{1} = rhor; rhordata{1} = rhordata_iele;
phirdata_iele{1} = phir.*r; phirdata{1} = phirdata_iele;

% write FS potential in LAMMPS format
writeLAMMPSeamfs(potfilename,griddata,elename,elemdata,rhordata,phirdata,Frhodata)

% convert LAMMPS format to MD++ format
[griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs(potfilename);
writeMDPPeam(outpotfile,griddata,elemdata{1},rhordata{1}{1},phirdata{1}{1},Frhodata{1});

