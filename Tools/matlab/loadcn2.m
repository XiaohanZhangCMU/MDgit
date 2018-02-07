function [np,s,h]=loadcn2(cnfilename)

% LOADCN    Load .cn type configuration files
%     For a given configuration file *.cn, LOADCN reads from the file
%     and get np,s,h.
%
% Usage:    [np,s,h] = loadcn(cnfilename)
%

fid = fopen(cnfilename,'r');
np=fscanf(fid,'%d',1);
s=zeros(np,10);  h=zeros(3,3); 
for i=1:np,  s(i,:)=(fscanf(fid,'%f',10))'; end;
for i=1:3,   h(i,:)=(fscanf(fid,'%f',3 ))'; end;
fclose(fid);
