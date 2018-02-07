function[nNodes,boxsize, R] = loadCNfile(cnfile)
fid = fopen(cnfile,'r');
nNodes = fscanf(fid, '%d\n', 1);
R = zeros(nNodes,3);
boxsize =zeros(3,3);
for i = 1:nNodes
    for j = 1:3
        R(i,j) = fscanf(fid,'%f',1);
    end
end
for i = 1:3  
    for j = 1:3
        boxsize(i,j) = fscanf(fid,'%f',1);
    end
end
fclose(fid);
end