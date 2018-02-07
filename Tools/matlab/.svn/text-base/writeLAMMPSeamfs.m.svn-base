function writeLAMMPSeamfs(potfilename,griddata,elename,elemdata,rhordata,phirdata,Frhodata)


fid = fopen(potfilename,'w');

fprintf(fid,'EAM/FS Mo potential for LAMMPS written by writeLAMMPSeamfs.m\n');
fprintf(fid,'G. J. Ackland and R. Thetford, Phil. Mag. A. (1987) vol.56 pp.15\n');
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


Nrho = griddata.Nrho;
drho = griddata.drho;
Nr   = griddata.Nr;
dr   = griddata.dr;
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
    for j=1:5:Nrho
        if j < 20
            disp(sprintf('%26.16E %26.16E %26.16E %26.16E %26.16E',Frho(j:j+4)));
        end
        fprintf(fid,'%26.16E %26.16E %26.16E %26.16E %26.16E\n',...
                Frho(j:j+4));
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
        for k=1:5:Nr
        if k < 20
            disp(sprintf('%26.16E %26.16E %26.16E %26.16E %26.16E',rhor(k:k+4)));
        end
        fprintf(fid,'%26.16E %26.16E %26.16E %26.16E %26.16E\n',...
                rhor(k:k+4));
            
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
        for k=1:5:Nr
        if k < 20
            disp(sprintf('%26.16E %26.16E %26.16E %26.16E %26.16E',phir(k:k+4)));
        end
        fprintf(fid,'%26.16E %26.16E %26.16E %26.16E %26.16E\n',...
                phir(k:k+4));
            
        end        
    end
    
end
fprintf(fid,'\n');
fclose(fid);
