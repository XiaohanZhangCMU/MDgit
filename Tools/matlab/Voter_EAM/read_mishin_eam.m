function [r, rho, rhor, phir, Frho] = read_mishin_eam(filename, nr, nrho)

fid = fopen(filename);

comment = fgets(fid);

ielem = fscanf(fid,'%d',1);
amass = fscanf(fid,'%f',1);
blatt = fscanf(fid,'%f',1);

dr   = fscanf(fid,'%f',1);
drho = fscanf(fid,'%f',1);
rcut = fscanf(fid,'%f',1);
rmin = fscanf(fid,'%f',1);

%ielem,amass,blatt,nrho,drho,nr,dr,rcut

data = fscanf(fid,'%f %f %f %f',[4, nr]); data = data';
rhor = data(:,1);
data = fscanf(fid,'%f %f %f %f',[4, nr]); data = data';
phir = data(:,1);
data = fscanf(fid,'%f %f %f %f',[4, nr]); data = data';
Frho = data(:,1);

fclose(fid);

r = [0:nr-1]*dr;
rho = [0:nrho-1]*drho;

return