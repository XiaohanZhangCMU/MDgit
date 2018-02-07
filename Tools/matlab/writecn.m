function writecn(filename,np,s,h)

% WRITECN    write .cn type configuration files
%     based on np,s,h
%
% Usage:    writecfg(cfgfilename)
%

% write number of atoms
fid = fopen(filename,'w');
fprintf(fid,'%d\n',np);

% write scaled coordinates
for i=1:np,
  fprintf(fid,' %8.6f %8.6f %8.6f \n',s(i,1),s(i,2),s(i,3));
end

% write box matrix H
fprintf(fid,' %8.6f %8.6f %8.6f \n',h(1,1),h(1,2),h(1,3));
fprintf(fid,' %8.6f %8.6f %8.6f \n',h(2,1),h(2,2),h(2,3));
fprintf(fid,' %8.6f %8.6f %8.6f \n',h(3,1),h(3,2),h(3,3));

% write the species line at the end
fprintf(fid,'1 Mo\n');
 
% write the two zeros at the very end
fprintf(fid,'%f %f',0,0);

fclose(fid);
