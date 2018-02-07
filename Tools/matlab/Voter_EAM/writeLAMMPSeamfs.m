function writeLAMMPSeamfs(potfilename,griddata,elename,elemdata,rhordata,phirdata,Frhodata)


fid = fopen(potfilename,'w');

fprintf(fid,'EAM/FS potential for LAMMPS written by writeLAMMPSeamfs.m\n');
fprintf(fid,'                                                         \n');
fprintf(fid,'Produced on %s\n',date);

Nele = griddata.Nele;
tmp = sprintf('%d', Nele);
for i=1:Nele
    ename = elename{i};
    tmp = [tmp, sprintf('  %s',ename)];
end

% write number of elements and their names
disp(sprintf('%s',tmp));
fprintf(fid,'%s\n',tmp);


% test whether the elements are the same
same_element = 1;
for i=2:Nele,
    if ~strcmp(elename{i-1},elename{i})
        same_element = 0;
    end
end


rmin = griddata.rmin;
dr   = griddata.dr;
N0   = rmin/dr
Nrho = griddata.Nrho;
drho = griddata.drho;
Nr   = griddata.Nr + N0;
rcut = griddata.rcut;

% write grid info
disp(sprintf('%d  %26.16E  %d  %26.16E  %26.16E',Nrho,drho,Nr,dr,rcut));
fprintf(fid,'%d  %26.16E  %d  %26.16E  %26.16E\n',Nrho,drho,Nr,dr,rcut);

for i=1:Nele
    Z = elemdata{i}.Z;
    atommass = elemdata{i}.atommass;
    latt = elemdata{i}.latticeconst;
    latttype = elemdata{i}.latttype;
    
    % write element info
    disp(sprintf('%d  %8.4f  %8.4f  %s',Z,atommass,latt,latttype));
    fprintf(fid,'%d  %8.4f  %8.4f  %s\n',Z,atommass,latt,latttype);

    
    
    % write F(rho)
    Frho = Frhodata{i}; 
    disp(sprintf('writing F(rho) of element %d',i));
    for j=1:Nrho
        if j < 20
            disp(sprintf('%26.16E',Frho(j)));
        end
        fprintf(fid,'%26.16E\n', Frho(j));
    end
    disp(' ........... ');
    
    % write rho(r)
    rhordata_iele = rhordata{i};
    disp(sprintf('writing rho(r) of element %d',i));
    for j=1:Nele
        
        if same_element && (i ~= j)
            continue;
        end
        
        rhor = rhordata_iele{j};
        for k=1:Nr
            if k <= N0
                fprintf(fid,'%26.16E\n', rhor(1));
            else
                if k < (N0 + 20)
                    disp(sprintf('%26.16E',rhor(k-N0)));
                end
                fprintf(fid,'%26.16E\n', rhor(k-N0));
            
            end
        end   
    end
    disp(' ........... ');
    
end


for i=1:Nele
    % write phi(r)
    phirdata_iele = phirdata{i};
    for j=i:Nele
        if same_element && (i ~= j)
            continue;
        end        
        disp(sprintf('writing phi(r) between element %d and %d',i,j));
    
        phir = phirdata_iele{j};
        for k=1:Nr
            if k <= N0
                fprintf(fid,'%26.16E\n', phir(1)/rmin*k*dr);
            else
                if k < (N0 + 20)
                    disp(sprintf('%26.16E',phir(k-N0)));
                end
                fprintf(fid,'%26.16E\n', phir(k-N0));
            end 
        end        
    end
    
end
fprintf(fid,'\n');
fclose(fid);

