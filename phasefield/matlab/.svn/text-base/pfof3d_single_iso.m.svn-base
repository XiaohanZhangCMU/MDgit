%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File: pfof3d_single_iso.m
%
%Purpose: 3D isotropic phase field orientational field model
%
% Yanming Wang (yanmingw@stanford.edu)
% Wei Cai (caiwei@stanford.edu)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D single phase field model with isotropic surface energy
gradq_eps = 1e-4;

if(~exist('N_grid'))
    N_grid = [200 200 2];
end

if(~exist('model_param'))
    model_param = [0.1 1.0 1.0 0.01 0.1 0.05 1e-2]; 
end

h = model_param(1);      % grid size
M = model_param(2);      % kinetic coefficient
U = model_param(3);      % energy parameters
vareps = model_param(4); % energy parameters
mu = model_param(5);     % chemical potential
H = model_param(6);      % grain boundary energy
Mq = model_param(7);     % mobility for orientation field


if(~exist('Niter'))
    Niter = 10000; % total number of simulation steps
end

if(~exist('dt'))
    dt = 1e-3; % time step
end

if(~exist('equation_type'))
    equation_type = 0; %0: Ginzburg-Landau, 1: Cahn-Hillard
end

if(~exist('plotfreq'))
    plotfreq = 200;
end

if(~exist('savefreq'))
   savefreq = 5000;
end

if(~exist('printfreq'))
    printfreq = 50;
end

cell_size = N_grid * h;
x = [0:N_grid(1)-1] * h;
y = [0:N_grid(2)-1] * h;
z = [0:N_grid(3)-1] * h;
Fdata = zeros(Niter,1);
Gdata = zeros(Niter,1);

if(~exist('phi'))
    phi = rand(N_grid);         
end

%%%%%%%%%%%%%%%
% initialize the orientation field
if (~exist('eta'))
    eta = rand(N_grid)*pi;  
end

if (~exist('c'))
    c = [1 1 1]/sqrt(3);
end

if(~exist('q'))
    q{1} = cos(eta/2);
    q{2} = c(1)*sin(eta/2);
    q{3} = c(2)*sin(eta/2);
    q{4} = c(3)*sin(eta/2);
end

% Choice 1 tanh function
% vol_func = inline('tanh((x-0.5)/0.10)/2+0.5', 'x');
% grad_vol_func = inline('(1-tanh((x-0.5)/0.10).^2)/2/0.10', 'x');
% Choice 2 polynomial function
vol_func = inline('x.^3.*(10-15*x+6*x.^2)', 'x');
grad_vol_func = inline('30*x.^4-60*x.^3+30*x.^2', 'x');

pxind=[2:N_grid(1),1]; nxind=[N_grid(1),1:N_grid(1)-1];
pyind=[2:N_grid(2),1]; nyind=[N_grid(2),1:N_grid(2)-1];
pzind=[2:N_grid(3),1]; nzind=[N_grid(3),1:N_grid(3)-1];

%iteration begins
for iter = 1:Niter,
    dphidx = (phi(pxind,:,:)-phi(nxind,:,:))/(2*h);
    dphidy = (phi(:,pyind,:)-phi(:,nyind,:))/(2*h);
    dphidz = (phi(:,:,pzind)-phi(:,:,nzind))/(2*h);
    d2phi  = (phi(pxind,:,:)+phi(nxind,:,:)+phi(:,pyind,:)+phi(:,nyind,:) ...
           + phi(:,:,pzind)+phi(:,:,nzind)-6*phi)/(h^2);
         
    for i = 1:4,         
        dqdx{i} = (q{i}(pxind,:,:)-q{i}(nxind,:,:))/(2*h);
        dqdy{i} = (q{i}(:,pyind,:)-q{i}(:,nyind,:))/(2*h);
        dqdz{i} = (q{i}(:,:,pzind)-q{i}(:,:,nzind))/(2*h);
        d2q{i}  = (q{i}(pxind,:,:)+q{i}(nxind,:,:)+q{i}(:,pyind,:)+q{i}(:,nyind,:) ...
                + q{i}(:,:,pzind)+q{i}(:,:,nzind)-6*q{i})/(h^2);
    end      

    vol_phi = vol_func(phi);
    dvol_phidx = (vol_phi(pxind,:,:)-vol_phi(nxind,:,:))/(2*h);
    dvol_phidy = (vol_phi(:,pyind,:)-vol_phi(:,nyind,:))/(2*h);
    dvol_phidz = (vol_phi(:,:,pzind)-vol_phi(:,:,nzind))/(2*h);
    
    
    F = (U * sum(sum(sum( phi.^2.*(1-phi).^2 ))) ...
        + vareps * sum(sum(sum( dphidx.^2 + dphidy.^2 + dphidz.^2))) ...
        + mu*sum(sum(sum(vol_phi)))) * (h^3);
    
    gradq_normsqr = zeros(N_grid);
    for i = 1:4,
        grad_qsqr{i} = dqdx{i}.^2 + dqdy{i}.^2 + dqdz{i}.^2;
        gradq_normsqr = gradq_normsqr + grad_qsqr{i};
    end
   
    gradq_normsqr = gradq_normsqr + gradq_eps^2; 
    gradq_norm = sqrt(gradq_normsqr);
    
    F_ori = H*sum(sum(sum(vol_phi.*gradq_norm))) * (h^3);
        
    F = F + F_ori;
    
    Fdata(iter)=F;
    
    dFdphi = 2*U*phi.*(1-phi).*(1-2*phi) - 2*vareps * d2phi ...
           + mu*grad_vol_func(phi) ...
           + H*grad_vol_func(phi).*gradq_norm;
    
    % implementation 1        
    for i = 1:4,
        dFdq{i} = -H*(dvol_phidx.*dqdx{i}+dvol_phidy.*dqdy{i}+dvol_phidz.*dqdz{i})./gradq_norm ...
                - H*vol_phi./gradq_norm.*(d2q{i} - grad_qsqr{i}./gradq_normsqr);
    end
    
    dFdq_dot_q = zeros(N_grid);
    
    for i = 1:4,
        dFdq_dot_q = dFdq_dot_q + dFdq{i}.*q{i};
    end
    
    if (equation_type == 0)
        % Ginzburg-Landau
        dphidt = -M*dFdphi;        
        for i = 1:4,
            dqdt{i} = -Mq*(dFdq{i} - dFdq_dot_q.*q{i});
        end
        
    elseif (equation_type == 1)
        
        d2dFdphi = (dFdphi(pxind,:,:)+dFdphi(nxind,:,:) ...
                   +dFdphi(:,pyind,:)+dFdphi(:,nyind,:) ...
                   +dFdphi(:,:,pzind)+dFdphi(:,:,nzind)-6*dFdphi)/(h^2);
        dphidt = M*d2dFdphi;
        for i = 1:4,
            dqdt{i} = -Mq*(dFdq{i} - dFdq_dot_q.*q{i});
        end
    else
        disp(sprintf('Unknown equation_type %d',equation_type));
        return
    end
    if mod(iter,2) == 0 || iter == 1,
        phi = phi + dt*dphidt;
    end
        
    for i = 1:4,
        q{i} = q{i} + dt*dqdt{i};
    end
    
    qnorm = sqrt(q{1}.^2 + q{2}.^2 + q{3}.^2 + q{4}.^2);
    for i = 1:4,
        q{i} = q{i}./qnorm;
    end
    
    Gdata(iter) = max(max(max(abs(dphidt))));
    
    if(mod(iter,printfreq)==0)
        disp(sprintf('iter=%d F=%14.6e G=%14.6e',iter,F,Gdata(iter)));
    end
    
    if(mod(iter,plotfreq)==0)
        
        % draw the phase field
        figure(1);
        clf;
        mesh(x,y,phi(:,:,1));
        zlim([0 1.0]);
        if mod(iter,savefreq) == 0,
            figname1 = sprintf('phase%d',iter);
            print(figname1,'-djpeg');
        end
       
        % draw orientation angle
        eta = 2*acos(q{1});
        figure(2);
        clf;
        mesh(x,y,eta(:,:,1));
        colorbar;
        drawnow;
        if mod(iter,savefreq) == 0,
            figname1 = sprintf('angle%d',iter);
            print(figname1,'-djpeg');
        end
       
        % plot grain color
        Graincolor = zeros(N_grid(1),N_grid(2),3);
        Graincolor(:,:,1) = (q{1}(:,:,1)+1)/2*255;
        Graincolor(:,:,2) = (q{2}(:,:,1)+1)/2*255;
        Graincolor(:,:,3) = (q{3}(:,:,1)+1)/2*255;
        Graingray = 0.299*Graincolor(:,:,1) + 0.587*Graincolor(:,:,2) + 0.114*Graincolor(:,:,3);
        
        figure(3);
        clf;
        contour(x,y,Graingray);
            
        figure(4);
        subplot(2,1,1); 
        plot([0:iter-1]*dt,Fdata(1:iter));
        xlabel('t');
        ylabel('F');
        subplot(2,1,2);
        plot([0:iter-1]*dt,Gdata(1:iter));
        xlabel('t');
        ylabel('max(G)');        
        drawnow;
        if mod(iter,Niter) == 0,
            figname1 = sprintf('F_G');
            print(figname1,'-djpeg');
        end
    end
end






