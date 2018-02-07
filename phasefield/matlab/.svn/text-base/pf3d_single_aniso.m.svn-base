%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File: pf3d_single_aniso.m
%
%Purpose: 3-dimensional Phase Field model 
%         solving Ginzburg-Landau (or Cahn-Hillard) equation
%         by finite difference method
%         Anisotropic surface energy (cubic symmetry)
%        
% Seunghwa Ryu (shryu@stanford.edu)
% Wei Cai (caiwei@stanford.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D single phase field model with anisotropic surface energy

if(~exist('N_grid'))
    N_grid = [50 50 50];
end

if(~exist('model_param'))
    % surface energy parameters cooresponding to anisotropy of Si 
    model_param = [0.4 1.0 2.0   0.0   1.0 0.6993 -3.6595 -1.7384];
end

h = model_param(1);      % grid size
M = model_param(2);      % kinetic coefficient
U = model_param(3);      % energy parameters
m = model_param(4);      % chemical potential difference between two phases
eps0 = model_param(5);   % energy parameters
eps1 = model_param(6);   % energy parameters
eps2 = model_param(7);   % energy parameters
eps3 = model_param(8);   % energy parameters

if(~exist('rot_matrix'))
    % contains orientations of x, y, z axes in crystal coordinates
    %rot_matrix = eye(3);
    %rot_matrix = [ [1 1 1]'/sqrt(3), [1 -1 0]'/sqrt(2), [1 1 -2]'/sqrt(6)];
    %rot_matrix = [ [1 1 -2]'/sqrt(6), [1 1 1]'/sqrt(3), [1 -1 0]'/sqrt(2)];
    rot_matrix = [ [1 -1 0]'/sqrt(2), [1 1 -2]'/sqrt(6), [1 1 1]'/sqrt(3)];
end

if(~exist('Niter'))
    Niter = 2000; % total number of simulation steps
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
    phi = zeros(N_grid);
    for i=1:N_grid(1), for j=1:N_grid(2), for k=1:N_grid(3),
       if (i-N_grid(1)/2)^2+(j-N_grid(2)/2)^2+(k-N_grid(3)/2)^2 < (max(N_grid(1:2))*0.3)^2
           phi(i,j,k) = 1.0 - 1e-6;
       end
    end; end; end
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
    % spatial derivatives
    dphidx = (phi(pxind,:,:)-phi(nxind,:,:))/(2*h);
    dphidy = (phi(:,pyind,:)-phi(:,nyind,:))/(2*h);
    dphidz = (phi(:,:,pzind)-phi(:,:,nzind))/(2*h);
    d2phi  = (phi(pxind,:,:)+phi(nxind,:,:)+phi(:,pyind,:)+phi(:,nyind,:) ...
             +phi(:,:,pzind)+phi(:,:,nzind)-6*phi)/(h^2);

    dphi_SQR = (dphidx.*dphidx) + (dphidy.*dphidy) + (dphidz.*dphidz);
    absdphi = sqrt(dphi_SQR + 1e-24);
    dphi_CUBE = absdphi.^3;
    
    nx_lab = dphidx ./ absdphi;
    ny_lab = dphidy ./ absdphi;
    nz_lab = dphidz ./ absdphi;
    
    R = rot_matrix;
    % crystal <- lab
    nx_cryst = R(1,1)*nx_lab + R(1,2)*ny_lab + R(1,3)*nz_lab;
    ny_cryst = R(2,1)*nx_lab + R(2,2)*ny_lab + R(2,3)*nz_lab;
    nz_cryst = R(3,1)*nx_lab + R(3,2)*ny_lab + R(3,3)*nz_lab;
    
    triplesum = (nx_cryst.^2).*(ny_cryst.^2) ...
              + (ny_cryst.^2).*(nz_cryst.^2) ...
              + (nz_cryst.^2).*(nx_cryst.^2);
    EPS_loc   = eps0 + eps1.*triplesum + eps2.*(nx_cryst.*ny_cryst.*nz_cryst).^2 ...
                     + eps3*(triplesum.^2);
                     
    % vareps in 'iso.m' now becomes (0.5 EPS_loc.^2)
    F = (U * sum(sum(sum( phi.^2.*(1-phi).^2 ))) ...
        - m * sum(sum(sum( constr_func(phi) ))) ...
        + 0.5 * sum(sum(sum( EPS_loc.^2 .* dphi_SQR)))) * (h^3);
    Fdata(iter)=F;
   
    dnx_lab_dphidx = 1./absdphi - dphidx .* dphidx ./ dphi_CUBE;
    dny_lab_dphidx =            - dphidy .* dphidx ./ dphi_CUBE;
    dnz_lab_dphidx =            - dphidz .* dphidx ./ dphi_CUBE;

    dnx_lab_dphidy =            - dphidx .* dphidy ./ dphi_CUBE;
    dny_lab_dphidy = 1./absdphi - dphidy .* dphidy ./ dphi_CUBE;
    dnz_lab_dphidy =            - dphidz .* dphidy ./ dphi_CUBE;

    dnx_lab_dphidz =            - dphidx .* dphidz ./ dphi_CUBE;
    dny_lab_dphidz =            - dphidy .* dphidz ./ dphi_CUBE;
    dnz_lab_dphidz = 1./absdphi - dphidz .* dphidz ./ dphi_CUBE;

    dEPS_dnx_cryst = 2.*eps1*nx_cryst.*((ny_cryst.^2)+(nz_cryst.^2)) ...
             + 2.*eps2*nx_cryst.*(ny_cryst.^2).*(nz_cryst.^2)  ...
             + 4.*eps3*nx_cryst.*((ny_cryst.^2)+(nz_cryst.^2)).*triplesum;
    dEPS_dny_cryst = 2.*eps1*ny_cryst.*((nz_cryst.^2)+(nx_cryst.^2)) ...
             + 2.*eps2*ny_cryst.*(nz_cryst.^2).*(nx_cryst.^2)  ...
             + 4.*eps3*ny_cryst.*((nz_cryst.^2)+(nx_cryst.^2)).*triplesum;
    dEPS_dnz_cryst = 2.*eps1*nz_cryst.*((nx_cryst.^2)+(ny_cryst.^2)) ...
             + 2.*eps2*nz_cryst.*(nx_cryst.^2).*(ny_cryst.^2)   ...
             + 4.*eps3*nz_cryst.*((nx_cryst.^2)+(ny_cryst.^2)).*triplesum;

    % lab <- crystal
    dEPS_dnx_lab = R(1,1)*dEPS_dnx_cryst + R(2,1)*dEPS_dny_cryst + R(3,1)*dEPS_dnz_cryst;
    dEPS_dny_lab = R(1,2)*dEPS_dnx_cryst + R(2,2)*dEPS_dny_cryst + R(3,2)*dEPS_dnz_cryst;
    dEPS_dnz_lab = R(1,3)*dEPS_dnx_cryst + R(2,3)*dEPS_dny_cryst + R(3,3)*dEPS_dnz_cryst;
    
    dEPS_dphidx = dEPS_dnx_lab.*dnx_lab_dphidx + dEPS_dny_lab.*dny_lab_dphidx + dEPS_dnz_lab.*dnz_lab_dphidx;
    dEPS_dphidy = dEPS_dnx_lab.*dnx_lab_dphidy + dEPS_dny_lab.*dny_lab_dphidy + dEPS_dnz_lab.*dnz_lab_dphidy;
    dEPS_dphidz = dEPS_dnx_lab.*dnx_lab_dphidz + dEPS_dny_lab.*dny_lab_dphidz + dEPS_dnz_lab.*dnz_lab_dphidz;

%     tmp_x = (EPS_loc.^2).*dphidx + dphi_SQR.*EPS_loc.*dEPS_dphidx;
%     tmp_y = (EPS_loc.^2).*dphidy + dphi_SQR.*EPS_loc.*dEPS_dphidy;
%     tmp_z = (EPS_loc.^2).*dphidz + dphi_SQR.*EPS_loc.*dEPS_dphidz;    
%     
%     dFdphi = 2*U*phi.*(1-phi).*(1-2*phi) ...
%            - m*grad_constr_func(phi) ...
%            - (tmp_x(pxind,:,:)-tmp_x(nxind,:,:))/(2*h) ...
%            - (tmp_y(:,pyind,:)-tmp_y(:,nyind,:))/(2*h) ...
%            - (tmp_z(:,:,pzind)-tmp_z(:,:,nzind))/(2*h) ;

    tmp_x = dphi_SQR.*EPS_loc.*dEPS_dphidx;
    tmp_y = dphi_SQR.*EPS_loc.*dEPS_dphidy;
    tmp_z = dphi_SQR.*EPS_loc.*dEPS_dphidz;    
    
    dFdphi = 2*U*phi.*(1-phi).*(1-2*phi) ...
           - m*grad_constr_func(phi) ...
           - (EPS_loc.^2).*d2phi ...
           - 2*EPS_loc.*(EPS_loc(pxind,:,:)-EPS_loc(nxind,:,:))/(2*h).*dphidx ...
           - 2*EPS_loc.*(EPS_loc(:,pyind,:)-EPS_loc(:,nyind,:))/(2*h).*dphidy ...
           - 2*EPS_loc.*(EPS_loc(:,:,pzind)-EPS_loc(:,:,nzind))/(2*h).*dphidz ...
           - (tmp_x(pxind,:,:)-tmp_x(nxind,:,:))/(2*h) ...
           - (tmp_y(:,pyind,:)-tmp_y(:,nyind,:))/(2*h) ...
           - (tmp_z(:,:,pzind)-tmp_z(:,:,nzind))/(2*h) ;
       
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
        xlabel('y'); ylabel('x'); zlabel('z');
        
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

