function [np,data,h]=loadcfg(filename)

% LOADCN    Load .cfg type configuration files
%     For a given configuration file *.cfg, LOADCFG reads from the file
%     and get np,s,h,pote
%
% Usage:    [np,s,h] = loadcfg(cnfilename)
%

% find number of atom
fid = fopen(filename,'r');
line = fgets(fid); np = sscanf(line(23:end),'%d');
line = fgets(fid);

% find box size
h=zeros(3,3);
line = fgets(fid); h(1,1) = sscanf(line(11:end),'%f');
line = fgets(fid); h(1,2) = sscanf(line(11:end),'%f');
line = fgets(fid); h(1,3) = sscanf(line(11:end),'%f');
line = fgets(fid); h(2,1) = sscanf(line(11:end),'%f');
line = fgets(fid); h(2,2) = sscanf(line(11:end),'%f');
line = fgets(fid); h(2,3) = sscanf(line(11:end),'%f');
line = fgets(fid); h(3,1) = sscanf(line(11:end),'%f');
line = fgets(fid); h(3,2) = sscanf(line(11:end),'%f');
line = fgets(fid); h(3,3) = sscanf(line(11:end),'%f');

% find number of entries per atom
line = fgets(fid); %.NO_VELOCITY.
line = fgets(fid); nc = sscanf(line(15:end),'%d');
line = fgets(fid); %auxiliary[0] = pote [eV]

% atomic data
data = zeros(np,nc);
for i=1:np,
    line = fgets(fid);  %atomic  mass
    line = fgets(fid);  %element name
    line = fgets(fid);  data(i,:) = sscanf(line,'%f')'; 
end

fclose(fid);
