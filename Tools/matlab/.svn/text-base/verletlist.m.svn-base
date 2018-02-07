function [nn,index] = verletlist(r,h,rv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: verletlist.m
%%
%% Purpose: Construct Verlet List (neighbor list) in 3D
%%          uses cell list to achieve O(N)
%%
%% Input:
%%    r  - (real) coordinate of atoms
%%    h  - PBC box h = (c1|c2|c3)
%%    rv - Verlet cut-off radius
%%
%% Output:
%%    nn - 1D array, nn(i) is the number of neighbors for atom i
%%    index - 2D array, index(i,j) is the index of jth neighbor
%%            of atom i, j = 1,...,nn(i)
%%
%% ME346 Introduction to Molecular Simulations (Stanford)
%% Wei Cai (caiwei@stanford.edu)
%% generalized from Solution to Homework 4.3
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first determine size of the cell list
c1 = h(:,1);
c2 = h(:,2);
c3 = h(:,3);
V = abs(det(h));
hx = abs( V / norm(cross(c2,c3)) );
hy = abs( V / norm(cross(c3,c1)) );
hz = abs( V / norm(cross(c1,c2)) );

Nx = floor(hx/rv);
Ny = floor(hy/rv);
Nz = floor(hz/rv);

NMAX = 50;
CMAX = 20;

% inverse of h matrix
hinv = inv(h);

if (Nx<3)|(Ny<3)|(Nz<3)
    disp(sprintf('reconstruct cell list Nx=%d Ny=%d Nz=%d',Nx,Ny,Nz));
end
[cell,cellid] = celllist(r,h,Nx,Ny,Nz,CMAX);

% find number of atoms
NP = length(r(:,1));

% initialize Verlet list
nn = zeros(NP,1);
index = zeros(NP,NMAX);

for i=1:NP,
   % position of atom i
   ri = r(i,:);
   
   % find which cell (ix,iy,iz) that atom i belongs to
   ix = cellid(i,1);
   iy = cellid(i,2);
   iz = cellid(i,3);
   
   % go through all neighboring cells
   for nx = ix-1:ix+1,
     for ny = iy-1:iy+1,
       for nz = iz-1:iz+1,
         nnx = nx; nny = ny; nnz = nz;
         
         % apply periodic boundary condition
         if(nnx<1) nnx=Nx; end
         if(nnx>Nx) nnx=1; end
         if(nny<1) nny=Ny; end
         if(nny>Ny) nny=1; end
         if(nnz<1) nnz=Nz; end
         if(nnz>Nz) nnz=1; end
         
         % extract atom id in this cell
         nc = cell(nnx,nny,nnz,1);
         ind = zeros(nc,1);
         ind(:) = cell(nnx,nny,nnz,2:nc+1);
         
         for k = 1:nc,
            j = ind(k);
            if (i==j)
            else
               rj = r(j,:);
               drij = ri - rj;
               dsij = (hinv * drij')';
               dsij = dsij - round(dsij);
               drij = (h * dsij')';
               if(norm(drij)<rv)
                  nn(i) = nn(i) + 1;
                  if(nn(i)>NMAX)
                     % extend memory
                     index2 = index;
                     index = zeros(NP,NMAX*2);
                     index(:,1:NMAX) = index2(:,:)
                     NMAX = NMAX * 2;
                  end               
                  index(i,nn(i)) = j;
               end
            end            
         end
        end
     end
   end
end
