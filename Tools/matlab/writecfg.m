function writecfg(filename,np,s,pote,h,element,mass)

% WRITECFG    write .cfg type configuration files
%     based on np,s,h,pote
%
% Usage:    writecfg(cfgfilename)
%

% find number of atom
fid = fopen(filename,'w');
fprintf(fid,'Number of particles = %d\n',np);
fprintf(fid,'A = 1.0 Angstrom (basic length-scale)\n');
for i=1:3,
    for j=1:3,
        fprintf(fid,'H0(%d,%d) = %9.6f A\n',i,j,h(i,j));
    end
end

fprintf(fid,'.NO_VELOCITY.\n');
fprintf(fid,'entry_count = 4\n');
fprintf(fid,'auxiliary[0] = pote [eV]\n');

for i=1:np,
  fprintf(fid,'%8.6f\n',mass);
  fprintf(fid,'%s\n',element);
  fprintf(fid,' %8.6f %8.6f %8.6f %8.6f\n',s(i,1),s(i,2),s(i,3),pote(i));
end

fclose(fid);
