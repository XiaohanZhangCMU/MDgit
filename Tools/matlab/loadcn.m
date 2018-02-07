function [np,s,h]=loadcn(cnfilename)

% LOADCN    Load .cn type configuration files
%     For a given configuration file *.cn, LOADCN reads from the file
%     and get np,s,h.
%
% Usage:    [np,s,h] = loadcn(cnfilename)
%


% (new code: every atom line can have arbitrary number of entries)
fid = fopen(cnfilename,'r');
myline = fgets(fid);
np=sscanf(myline,'%d',1);
s=zeros(np,3);  h=zeros(3,3); 
for i=1:np,  
    myline = fgets(fid);
    s(i,:)=(sscanf(myline,'%f',3))';
end;
for i=1:3,   
    myline = fgets(fid);
    h(i,:)=(sscanf(myline,'%f',3))'; 
end;
fclose(fid);

% (old code: require every atom line to have 3 entries)
%fopen(cnfilename,'r');
%np=fscanf(fid,'%d',1);
%s=zeros(np,3);  h=zeros(3,3); 
%for i=1:np,  s(i,:)=(fscanf(fid,'%f',3))'; end;
%for i=1:3,   h(i,:)=(fscanf(fid,'%f',3))'; end;
%fclose(fid);
