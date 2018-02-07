folder='~/Planet/Libs/MD++Fem/runs/x-test-fem-0-10-15/NEBinit';
nebchainfile = strcat(folder,'-neb.chain.500');
file = strcat(folder,num2str(0),'.cn');
chainpools = 50:50:1000;

[nNodes,box,R] = loadCNfile(file);
nchain = length(chainpools);
dlmwrite(nebchainfile,[nchain]);
dlmwrite(nebchainfile,[nNodes], '-append');
dlmwrite(nebchainfile,[0:nNodes-1]','-append', 'delimiter',' ');

neb1 = strcat(folder,'-1.cn');
neb2 = strcat(folder,'-2.cn');

file = strcat(neb1);
[nNodes,box,R] = loadCNfile(file);
dlmwrite(nebchainfile,R,'-append','delimiter',' ');

for iter = 1:nchain
chainno=chainpools(iter);
file = strcat(folder,num2str(chainno),'.cn');
[nNodes,box,R] = loadCNfile(file);
dlmwrite(nebchainfile,R,'-append','delimiter',' ');
end

file = strcat(neb2);
[nNodes,box,R] = loadCNfile(file);
dlmwrite(nebchainfile,R,'-append','delimiter',' ');
