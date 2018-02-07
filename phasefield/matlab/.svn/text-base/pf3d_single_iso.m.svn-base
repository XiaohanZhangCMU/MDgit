%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File: pf3d_single_iso.m
%
%Purpose: 3-dimensional Phase Field model 
%         solving Ginzburg-Landau (or Cahn-Hillard) equation
%         by finite difference method
%        
% Wei Cai (caiwei@stanford.edu)
% Seunghwa Ryu (shryu@stanford.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D single phase field model with isotropic surface energy

if(~exist('N_grid'))
    N_grid = [50 50 50];
end

if(~exist('model_param'))
    model_param = [0.4 1.0 1.0 1.0   0.0 ]; 
end

h = model_param(1);      % grid size
M = model_param(2);      % kinetic coefficient
U = model_param(3);      % energy parameters
vareps = model_param(4); % energy parameters
m = model_param(5);      % chemical potential difference between two phases

if(~exist('Niter'))
    Niter = 5000; % total number of simulation steps
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

constr_func = inline('x.^3/3 + (1-x).*x.^2 - (1-x).^3/3 - x.*(1-x).^2', 'x');
grad_constr_func = inline('4*x.*(1-x)', 'x');

if N_grid(1)>1
    pxind=[2:N_grid(1),1]; nxind=[N_grid(1),1:N_grid(1)-1];
else
    pxind=1; nxind=1;
end

if N_grid(2)>1
    pyind=[2:N_grid(2),1]; nyind=[N_grid(2),1:N_grid(2)-1];
else
    pyind=1; nyind=1;
end

if N_grid(3)>1
    pzind=[2:N_grid(3),1]; nzind=[N_grid(3),1:N_grid(3)-1];
else
    pzind=1; nzind=1;
end

%iteration begins
for iter = 1:Niter,
    dphidx = (phi(pxind,:,:)-phi(nxind,:,:))/(2*h);
    dphidy = (phi(:,pyind,:)-phi(:,nyind,:))/(2*h);
    dphidz = (phi(:,:,pzind)-phi(:,:,nzind))/(2*h);
    d2phi  = (phi(pxind,:,:)+phi(nxind,:,:)+phi(:,pyind,:)+phi(:,nyind,:) ...
             +phi(:,:,pzind)+phi(:,:,nzind)-6*phi)/(h^2);

    F = (U * sum(sum(sum( phi.^2.*(1-phi).^2 ))) ...
        - m * sum(sum(sum( constr_func(phi) ))) ...
        + vareps * sum(sum(sum( dphidx.^2 + dphidy.^2 + dphidz.^2)))) * (h^3);
    Fdata(iter)=F;
   
    % (not fully self-consistent, but numerically stable)
    dFdphi = 2*U*phi.*(1-phi).*(1-2*phi) - m*grad_constr_func(phi) - 2*vareps * d2phi;
    
    % (fully self-consistent, but numerically unstable)
    %tmp_x = 2*vareps.*dphidx ;
    %tmp_y = 2*vareps.*dphidy ;
    %tmp_z = 2*vareps.*dphidz ;    
    %
    %dFdphi = 2*U*phi.*(1-phi).*(1-2*phi) ...
    %       + 4*m*phi.*(phi-1) ...
    %       - (tmp_x(pxind,:,:)-tmp_x(nxind,:,:))/(2*h) ...
    %       - (tmp_y(:,pyind,:)-tmp_y(:,nyind,:))/(2*h) ...
    %       - (tmp_z(:,:,pzind)-tmp_z(:,:,nzind))/(2*h) ;
 
    if (equation_type == 0)
        % Ginzburg-Landau
        dphidt = -M*dFdphi;        
    elseif (equation_type == 1)
        d2dFdphi = (dFdphi(pxind,:,:)+dFdphi(nxind,:,:) ...
                   +dFdphi(:,pyind,:)+dFdphi(:,nyind,:) ...
                   +dFdphi(:,:,pzind)+dFdphi(:,:,nzind)-6*dFdphi)/(h^2);
        dphidt = M*d2dFdphi;
    elseif (equation_type == 2)
        % constraint: int phi(x) dx = const
        %dphidt = -M*(dFdphi - mean(mean(mean(dFdphi))));
        % constraint: int m(phi(x)) dx = const
        m_prime = grad_constr_func(phi);
        dphidt = -M*(dFdphi - sum(sum(sum(m_prime.*dFdphi)))/sum(sum(sum(m_prime.^2)))*m_prime);
    else
        disp(sprintf('Unknown equation_type %d',equation_type));
        return
    end
    
    phi = phi + dt* dphidt;
    Gdata(iter) = max(max(max(abs(dphidt))));
    
    if(mod(iter,printfreq)==0 || iter==1)
        disp(sprintf('iter=%d F=%14.6e G=%14.6e rel. vol. = (%.1f %%)', ...
            iter,F,Gdata(iter),sum(sum(sum((constr_func(phi)+1/3)*(3/2))))/prod(N_grid)*100));
    end
    
    if(mod(iter,plotfreq)==0 || iter==1)
        figure(1);
        clf
        p_surf = patch(isosurface(y,x,z,phi,0.5));
        isonormals(y,x,z,phi,p_surf)
        set(p_surf,'FaceColor','cyan','EdgeColor','none');
        daspect([1 1 1])
        view(3); axis tight
        camlight 
        lighting gouraud
        xlim([min(y) max(y)]); ylim([min(x) max(x)]); zlim([min(z) max(z)]);
        alpha(0.7);
        
        figure(2);
        subplot(2,1,1); 
        plot([0:iter-1]*dt,Fdata(1:iter));
        xlabel('t');
        ylabel('F');
        subplot(2,1,2);
        plot([0:iter-1]*dt,Gdata(1:iter));
        xlabel('t');
        ylabel('max(G)');        
    end
end

