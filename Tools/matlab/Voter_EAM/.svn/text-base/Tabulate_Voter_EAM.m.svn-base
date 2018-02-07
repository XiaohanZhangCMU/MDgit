%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tabulate Voter EAM potential in MD++ EAM format                         %
%                                                                         %
% Keonwook Kang, kwkang@gmail.com                                         %
% Feb 10, 2011                                                            %
%                                                                         %
% Ref. 1.                                                                 %
% Arthur F. Voter, Los Alamos Unclassfied Technical Report                %
% # LA-UR 93-3901 (1993)                                                  %
% Ref. 2.                                                                 %
% Y. Mishin et al, PRB (2001) vol. 63, 224106                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In Ref. 1., A table of fitting parameters for Ni, Pd, Pt, Cu, Ag, Au, and
% Al is given.

%elename{1} = 'Ag';
elename{1} = 'Cu';

if (strcmp(elename{1}, 'Ag'))
    Nele = 1; Nrho = 5000; Nr = 5000;
    Z = 47; % atomic number
    atommass = 107.8682; % atomic mass
    potfilename = 'Ag_Voter.eam';
    latttype = 'FCC'; % lattice type
    outpotfile = 'eamdata.Ag.Voter';

    % Fitting Parameters for Ag
    D_M = 0.6721;           % in eV
    R_M = 2.5700;           % in Angstrom
    alpha = -1.8260;        % in 1/Angstrom
    beta = 3.906;           % in 1/Angstrom
    r_c = 5.5420;           % cut-off radius in Angstrom

    % Properties used in fit
    latticeconst = 4.09;        % in Angstrom
    E_coh = 2.85;               % cohesive energy in eV
    B = 1.04e12;                % bulk modulus in erg/cm^3
elseif (strcmp(elename{1}, 'Cu'))
    Nele = 1; Nrho = 5000; Nr = 5000;
    Z = 29; % atomic number
    atommass = 63.546; % atomic mass
    latttype = 'FCC'; % lattice type
    potfilename = 'Cu_Voter.eam';
    outpotfile = 'eamdata.Cu.Voter';

    % Fitting Parameters for Cu
    D_M = 0.7366;           % in eV
    R_M = 2.3250;           % in Angstrom
    alpha = -1.9190;        % in 1/Angstrom
    beta = 4.043;           % in 1/Angstrom
    r_c = 4.9610;           % cut-off radius in Angstrom

    % Properties used in fit
    latticeconst = 3.615;        % in Angstrom
    E_coh = 3.54;               % cohesive energy in eV
    B = 1.42e12;                % bulk modulus in erg/cm^3
else
    disp('unknown elename! do nothing.');
    return;
end

Omega = latticeconst^3/4;   % equilibrium atomic volume in Angstrom^3
% convert B from erg/cm^3 to eV/Angstrom^3
B = B*6.24150974/1e13;      % bulk modulus in eV/Angstrom^3

rmax    = round(r_c*1.1);
dr   = rmax/Nr;
r    = [0:Nr-1 ]*dr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Density function rho(r)
% Smoothly cut-off rho(r)
rhor_sm = rhor_sm_Voter(r,r_c,beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust the range of r.
% The values of r < r_min are forbidden in simulations.
% If r < r_min, rho(r) is not a monotonically decreasing function.
% Voter's rho(r) is not a monotonically decreasing fucntion as r increases.
% Thus, r should be larger than r_min at which rho(r) is max.
[rhor_sm_max, idx] = max(rhor_sm);
r_min = r(idx); r = r_min:dr:(r_min+dr*(Nr-1));

% reevaluate in the adjusted r
rhor  = rhor_Voter(r,beta);         % rho(r)            1st neighbor
rhor2 = rhor_Voter(sqrt(2)*r,beta); % rho(sqrt(2)*r)    2nd neighbor
rhor3 = rhor_Voter(sqrt(3)*r,beta); % rho(sqrt(3)*r)    3rd neighbor
rhor_rc = rhor_Voter(r_c,beta);     % rho(r_c)

% Smoothly cut-off rho(r)
rhor_sm = rhor_sm_Voter(r,r_c,beta);
rhor2_sm = rhor_sm_Voter(sqrt(2)*r,r_c,beta);
rhor3_sm = rhor_sm_Voter(sqrt(3)*r,r_c,beta);

b0 = latticeconst/sqrt(2);      % 1st neighbor
b1 = latticeconst;              % 2nd neighbor
b2 = latticeconst*sqrt(1.5);    % 3rd neighbor

% Normalize rhor_sm such that /bar{rho} to be 1 in the equilibrium FCC
% crystal or when r = latticeconst/sqrt(2).
rhor_sm_1 = rhor_sm_Voter(b0, r_c, beta);
rhor_sm_2 = rhor_sm_Voter(b1, r_c, beta);
rhor_sm_3 = rhor_sm_Voter(b2, r_c,beta);
rhobar_eq = 12*rhor_sm_1 + 6*rhor_sm_2 + 24*rhor_sm_3;

rhor_sm = (r<r_c).*rhor_sm/rhobar_eq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Morse type pair potential, phi(r)
% Smoothly cut-off phi(r)
% D_M : depth of the minimum
% R_M : position of the minimum
% alpha_M : a measure of the curvature at the minimum
phir_sm = phir_sm_Voter(r, r_c, alpha, R_M, D_M).*(r<r_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rose's universal energy, E_rose
a = 2/sqrt(2)*r; a_c = 2/sqrt(2)*r_c;
% reduced lattice constant a*
a_star = (a/latticeconst - 1)/(E_coh/(9*B*Omega))^(1/2);
acut_star = (a_c/latticeconst - 1)/(E_coh/(9*B*Omega))^(1/2);
% original Rose function
E_rose = -E_coh * Erose_fcn(a_star); 
E_rose_ac = -E_coh * Erose_fcn(acut_star);

% Iterative correction of Erose_fcn(a*) to satisfy Erose_fcn(acut*) = 0
% initial value of q and epsilon
q = acut_star; epsilon = Erose_fcn(q); epsilon_c = epsilon;
f_mod_acut_star = Erose_fcn_mod(acut_star, acut_star, epsilon);
disp('Start iteratively solving (1-epsilon)^.5*acut_star = q and f(q) = epsilon.')
for i=1:100
    if (mod(i-1,10)==0)
        disp(sprintf('i=%4d, q=%20.14e, epsilon = %20.14e, f_mod(a_c*) = %20.14e', ...
             i-1, q, epsilon, f_mod_acut_star))
    end
    
    % update of q and epsilon
    %q = .5*(q + (1-f_mod)^.5*acut_star);
    q = .5*(q + (1-epsilon)^.5*acut_star);
    epsilon = Erose_fcn(q);
    
    f_mod_acut_star = Erose_fcn_mod(acut_star, acut_star, epsilon);
end
disp(sprintf('i=%4d, q=%20.14e, epsilon = %20.14e, f_mod(a_c*) = %20.14e', ...
              i, q, epsilon, f_mod_acut_star))
disp('End')
E_rose_mod = -E_coh * Erose_fcn_mod(a_star, acut_star, epsilon);
E_rose_mod = (a<a_c).*E_rose_mod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Embedding function F(rho)
                     
% "interatomic_dst.fcc" file created by running cal_neighbor_dist_num.tcl
interatomic_dist = load('interatomic_dist.fcc'); 
interatomic_dist = interatomic_dist * latticeconst;
xin = unique(interatomic_dist);
[nx, xout ] = hist(interatomic_dist, xin);

idx = find((xout/(latticeconst/sqrt(2))*r_min) > r_c, 1, 'first');
phir_sum = 0; rhobar = 0;
disp('Start summing phi(r_i)*n_i over i, where n_i is the number of i-th closest neighbors') 
for i=1:idx-1
    coeff = xout(i)/(latticeconst/sqrt(2));
    disp(sprintf('i = %d, n_i = %d, r_i/r_1 = %f', i, nx(i), coeff))
    phir_sum = phir_sum + nx(i)*phir_sm_Voter(coeff*r, r_c, alpha, R_M, D_M)...
                          .*(coeff*r < r_c); 
    rhobar = rhobar + nx(i)*rhor_sm_Voter(coeff*r, r_c, beta)...
                          .*(coeff*r < r_c); 
end
disp('End')
Frho = E_rose_mod - 0.5*phir_sum;
rhobar = rhobar/rhobar_eq;

rho_max = max(round((rhobar_eq)* 50), round(max(rhobar)/10));
drho = rho_max/Nrho;
rho  = [0:Nrho-1]*drho;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform phi(r) and F(rho) in effective pair type
% See "Handbook of materials modeling" by S. Yip
% Ch. 2.2 Interatomic potentials for metals by Y. Mishin p.461-462
idx = find(abs(diff(rhobar))==0,1,'first');
Frho_interp = spline(rhobar(1:idx), Frho(1:idx), rho);

cs4 = spline(rho,Frho_interp);
cs4p=cs4; cs4p.order=3; cs4p.coefs=[cs4.coefs(:,1)*3,cs4.coefs(:,2)*2,cs4.coefs(:,3)];
%cs4pp=cs4; cs4pp.order=2; cs4pp.coefs=[cs4p.coefs(:,1)*2,cs4p.coefs(:,2)];
dFdrho = ppval(cs4p,1); %d2Fdrho2 = ppval(cs4pp,1);

phir_sm = phir_sm + 2*dFdrho*rhor_sm;
Frho_interp = Frho_interp - dFdrho*rho;


% In order to compare with Mishin Cu,
if (strcmp(elename{1}, 'Cu'))
    % Read Mishin Cu EAM table
    [r_mishin_cu, rho_mishin_cu, rhor_mishin_cu, phir_mishin_cu, Frho_mishin_cu ] = ...
        read_mishin_eam('./eamdata.Cu.Mishin', 5000, 5000);
    
    % Read Voter Cu EAM table
    fid = fopen('./eamdata.Cu.Voter.rhor');
    comment = fgets(fid); comment = fgets(fid);
    rhor_Cu_Voter = fscanf(fid,'%f %f',[2, 3000]); rhor_Cu_Voter = rhor_Cu_Voter';
    fclose(fid);
    
    fid = fopen('./eamdata.Cu.Voter.phir');
    comment = fgets(fid); comment = fgets(fid);
    phir_Cu_Voter = fscanf(fid,'%f %f',[2, 3000]); phi_Cu_Voter = phir_Cu_Voter';
    fclose(fid);
    
    fid = fopen('./eamdata.Cu.Voter.Frho');
    comment = fgets(fid); comment = fgets(fid);
    Frho_Cu_Voter = fscanf(fid,'%f %f',[2, 3000]); Frho_Cu_Voter = Frho_Cu_Voter';
    fclose(fid);
end

figure(1), set(gca, 'fontsize', 14)
if (strcmp(elename{1}, 'Cu'))
    plot(r, rhor_sm,'-',rhor_Cu_Voter(:,1),rhor_Cu_Voter(:,2),'r-.',...
         r_mishin_cu, rhor_mishin_cu)
    legend('Voter EAM Cu - analytic form', 'Voter EAM Cu - tabulated form', 'Mishin EAM Cu')    
    axis([2.2 6.0 0.0 0.16])
    line([b0; b0],[0; 0.16],'color','k'), text(b0+.05, 0.10, '1st neighbor')
    line([b1; b1],[0; 0.16],'color','k'), text(b1+.05, 0.09, '2nd neighbor')
    line([b2; b2],[0; 0.16],'color','k'), text(b2+.05, 0.08, '3rd neighbor')
    line([r_c; r_c],[0; 0.16],'color','m'), text(r_c+.05, 0.06, 'cut-off radius')
else
    plot(r, rhor_sm,'-')
end
xlabel('r'), ylabel('\rho(r)')

figure(2), set(gca, 'fontsize', 14)
if (strcmp(elename{1}, 'Cu'))
    plot(r, phir_sm, phi_Cu_Voter(:,1), phi_Cu_Voter(:,2), '-.', ...
         r_mishin_cu, phir_mishin_cu)     
    legend('Voter EAM Cu - analytic form', 'Voter EAM Cu - tabulated form', 'Mishin EAM Cu',4)    
    axis([2.2 6.0 -0.2 0.2])
    line([b0; b0],[-0.2; 0.2],'color','k'), text(b0+.05, 0.12, '1st neighbor')
    line([b1; b1],[-0.2; 0.2],'color','k'), text(b1+.05, 0.11, '2nd neighbor')
    line([b2; b2],[-0.2; 0.2],'color','k'), text(b2+.05, 0.10, '3rd neighbor')
    line([r_c; r_c],[-0.2; 0.2],'color','m'), text(r_c+.05, 0.09, 'cut-off radius')
else
    plot(r, phir_sm,'-')
end
xlabel('r'), ylabel('\phi(r)')

figure(3)
plot(a_star, E_rose, a_star, E_rose_mod), ylim([E_coh*(-2) E_coh*10])
xlabel('scaled seperation, a*'), ylabel('universal binding energy (eV)')
legend('original','modified')

figure(4), set(gca, 'fontsize', 14)
if (strcmp(elename{1}, 'Cu'))
    plot(rho, Frho_interp, Frho_Cu_Voter(:,1), Frho_Cu_Voter(:,2),'-.',...
         rho_mishin_cu, Frho_mishin_cu)
    legend('Voter EAM Cu - analytic form', 'Voter EAM Cu - tabulated form', 'Mishin EAM Cu')    
    axis([0.0 2.5 -2.5 0.0])
else
    plot(rho, Frho_interp)
end
xlabel('\rho'), ylabel('Embedding energy F(\rho) (eV)')

%figure(5),
%plot(r, rhobar)

griddata = struct('Nele',Nele,'Nrho',Nrho,'drho',drho,'Nr',Nr,'dr',dr,'rcut',r_c,'rmin',r_min);
elemdata{1} = struct('Z',Z,'atommass',atommass,'latticeconst',latticeconst,'latttype',latttype);
Frhodata{1} = Frho_interp;
rhordata_iele{1} = rhor_sm; rhordata{1} = rhordata_iele;
phirdata_iele{1} = phir_sm.*r; phirdata{1} = phirdata_iele;

% write EAM potential in MD++ format
writeMDPPeam(outpotfile,griddata,elemdata{1},rhordata{1}{1},phirdata{1}{1},Frhodata{1});

% write FS potential in LAMMPS format
writeLAMMPSeamfs(potfilename,griddata,elename,elemdata,rhordata,phirdata,Frhodata)
