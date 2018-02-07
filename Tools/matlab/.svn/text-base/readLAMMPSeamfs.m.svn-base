function [griddata, elemdata, Frhodata, rhordata, phirdata] = readLAMMPSeamfs(potfilename)
%read LAMMPS potential file of the eam.fs type
%notice: phirdata stores r*phi (eV*Angstrom)

fid = fopen(potfilename,'r');
% skip the first three lines
myline = fgets(fid);
myline = fgets(fid);
myline = fgets(fid);

% read number of elements
myline = fgets(fid);
tmp = strread(myline,'%s');
Nele = sscanf(myline,'%d',1);
for i=1:Nele,
    elename{i} = tmp(i+1,:)
end

% test whether the elements are the same
same_element = 1;
for i=2:Nele,
    if ~strcmp(elename{i-1},elename{i})
        same_element = 0;
    end
end

% read grid info
myline = fgets(fid);
tmp = sscanf(myline,'%e',5);
Nrho = tmp(1); drho = tmp(2); Nr = tmp(3); dr = tmp(4); rcut = tmp(5);
griddata = struct('Nele',Nele,'Nrho',Nrho,'drho',drho,'Nr',Nr,'dr',dr,'rcut',rcut);

disp(sprintf('Nele = %d  Nr = %d',Nele,Nr));

% loop over each element
for iele = 1:Nele,
    % read element info
    myline = fgets(fid);
    tmp = sscanf(myline,'%e',3);
    Z=tmp(1); atommass=tmp(2); latticeconst=tmp(3);
    elemdata{iele} = struct('Z',Z,'atommass',atommass,'latticeconst',latticeconst);
    
    % read F(rho)
    Frho = zeros(Nrho,1);
    disp(sprintf('reading F(rho) of element %d',iele));
    j = 0;
    while j < Nrho,
        myline = fgets(fid);
        tmp = sscanf(myline,'%e',inf);
        Frho(j+1:j+length(tmp)) = tmp;
        j = j+length(tmp);
    end
    Frhodata{iele} = Frho;
    %
    figure(iele);
    subplot(1,3,1);
    plot([0:Nrho-1],Frho);
    
    % read rho(r)
    disp(sprintf('reading rho(r) of element %d',iele));
    rhor = zeros(Nr,1);
    for jele = 1:Nele,
        if same_element && (iele ~= jele)
            continue;
        end
        n = 0;
        while n < Nr,
            myline = fgets(fid);
            %disp(sprintf('%s',myline));
            tmp = sscanf(myline,'%e',inf);
            rhor(n+1:n+length(tmp)) = tmp;
            n = n+length(tmp);
        end
        rhordata_iele{jele} = rhor;
    end
    rhordata{iele} = rhordata_iele;
    
    figure(iele);
    subplot(1,3,2);
    plot([0:Nr-1],rhor);
    if Nele == 3,
        subplot(1,3,2);
        plot([0:Nr-1],rhordata_iele{1},'-', ...
             [0:Nr-1],rhordata_iele{2},'.', ...
             [0:Nr-1],rhordata_iele{3},'--');
    end
end    
    
for iele = 1:Nele,
    % read phi(r)
    phir = zeros(Nr,1);
    for jele = 1:iele,
        if same_element && (iele ~= jele)
            continue;
        end
        disp(sprintf('reading phi(r) between element %d and %d',iele,jele));
        n = 0;
        while n < Nr,
            myline = fgets(fid);
            tmp = sscanf(myline,'%e',inf);
            phir(n+1:n+length(tmp)) = tmp;
            n = n+length(tmp);
        end
        
        phirdata_iele{jele} = phir;
        
        figure(iele);
        subplot(1,3,3);
        plot([0:Nr-1],phir);
        ylim([-10 100]);
    end
    phirdata{iele} = phirdata_iele;
end

for iele = 1:Nele,
    for jele = iele+1:Nele,
        phirdata{iele}{jele} = phirdata{jele}{iele};
    end
    figure(iele);
    if Nele == 2 && ~same_element,
        subplot(1,3,3);
        hold on
        plot([0:Nr-1],phirdata{iele}{1},'-', ...
             [0:Nr-1],phirdata{iele}{2},'-');
        ylim([-10 100]);
        hold off
    end
    if Nele == 3 && ~same_element,
        subplot(1,3,3);
        hold on
        plot([0:Nr-1],phirdata{iele}{1},'-', ...
             [0:Nr-1],phirdata{iele}{2},'-', ...
             [0:Nr-1],phirdata{iele}{3},'-');
        ylim([-10 100]);
        hold off
    end
end

fclose(fid);