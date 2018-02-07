function cn2poscar(cnfilename,title,fmt)
%
% convert MD++ cn file to VASP POSCAR file
%
% for example: 
%  cn2poscar('../../runs/si_au/L12-au3si-3x3x3-rand-rlx.cn')
%

%clear all

if(~exist('cnfilename'))
    cnfilename = '../../runs/si_au/L12-au3si-3x3x3-rand-rlx.cn';
end

if(~exist('ttl'))
    ttl = 'Au_3Si liquid';
end

if(~exist('fmt'))
    fmt = 0;  %0: direct, 1: cartesian
end
    
scalefactor = 1;

[token, remain] = strtok(fliplr(cnfilename),'.');
poscarfilename = strcat(fliplr(remain),'POSCAR');   % output filename

fid = fopen(cnfilename,'r');

disp(sprintf('read file: %s',cnfilename));
NP=fscanf(fid,'%d\n',1);
%%%%%%%%%%%%%%%  x y z vx vy vz Etot/ind fixed Topol. species 
cndata = fscanf(fid,'%e %e %e %e %e %e %e %d %e %d',[10 NP]); cndata = cndata';

species = cndata(:,10);
CN = fscanf(fid, '%e %e %e', [3 3]); CN = CN';

disp('convert cndata from scaled to real coordinates');
sr = cndata(:,1:3);

if (fmt~=0)
  r = (CN'*sr')';
  cndata(:,1:3) = r;
end


tol = 1e-10;
disp(sprintf('set H entries less than %e to zero',tol));
for i=1:3,
    for j=1:3,
        if(abs(CN(i,j))<tol)
            CN(i,j)=0;
        end
    end
end

NS = fscanf(fid, '%d', 1);  % number of species
for i=1:NS
    elements{i} = fscanf(fid, '%s', 1);
end

fclose(fid);

sortedcndata = sortrows(cndata,10); % separate cndata according to species
NPS = hist(cndata(:,10),0:NS-1);    % number of each species
NPSindex = [0 cumsum(NPS)];         % index 

disp(sprintf('write file: %s\n',poscarfilename));
fid = fopen(poscarfilename,'w');
fprintf(fid,'%s\n',ttl);
fprintf(fid,'%f\n',scalefactor);
fprintf(fid,'%25.16e %25.16e %25.16e\n %25.16e %25.16e %25.16e\n %25.16e %25.16e %25.16e\n', CN);
fprintf(fid,'%d %d\n', NPS);
%fprintf(fid,'%s\n','Selective Dynamics');

if(fmt==0)
  fprintf(fid,'%s\n','Cartesian');
else
  fprintf(fid,'%s\n','Direct');
end

for i=1:NS
  fprintf(fid, '%25.16e %25.16e %25.16e\n', sortedcndata(NPSindex(i)+1:NPSindex(i+1),1:3)');
end

fclose(fid);
