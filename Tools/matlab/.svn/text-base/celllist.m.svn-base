function [cell,cellid] = celllist(r,h,Nx,Ny,Nz,CMAX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: celllist.m
%%
%% Purpose: Construct cell list in 3D , i.e.,
%%          grouping atoms into Nx by Ny by Nz cells
%%
%% Input:
%%    r  - (real) coordinate of atoms
%%    h  - PBC box h = (c1|c2|c3)
%%    Nx - number of cells in x direction
%%    Ny - number of cells in y direction
%%    Nz - number of cells in z direction
%%    CMAX - (suggested) maximum number of atoms in each cell 
%%           for initial memory allocation
%%
%% Output:
%%    cell - the cell list nc=cell(i,j,k,1) is number of atoms in cell(i,j,k)
%%           cell(i,j,k,2:nc+1) are the index of atoms in cell(i,j,k)
%%    cellid - 2D array ix=cellid(i,1), iy=cellid(i,2), iz=cellid(i,3)
%%             atom i belong to cell(ix,iy,iz)
%%
%% ME346 Introduction to Molecular Simulations (Stanford)
%% Wei Cai (caiwei@stanford.edu)
%% generalized from Solution to Homework 4.3
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allocate cell memory
cell = zeros(Nx,Ny,Nz,CMAX+1);

% find number of atoms
NP = length(r(:,1));
cellid = zeros(NP,3);

% find reduced coordinates of all atoms
s = (inv(h) * r')';

% fold reduced coordinates into [0, 1)
s = s - floor(s);

for i=1:NP,
   ix = mod(floor(s(i,1)*Nx),Nx) + 1;
   iy = mod(floor(s(i,2)*Ny),Ny) + 1;
   iz = mod(floor(s(i,3)*Nz),Nz) + 1;
   cell(ix,iy,iz,1) = cell(ix,iy,iz,1) + 1;
   ind = cell(ix,iy,iz,1) + 1;
   if(ind>CMAX+1)
      % extend memory
      cell2 = cell;
      cell = zeros(Nx,Ny,Nz,CMAX*2+1);
      cell(:,:,:,1:CMAX+1)=cell2(:,:,:,:);
      CMAX = CMAX * 2;
   end
   
   % adding atom i to cell list
   cell(ix,iy,iz,ind) = i;
   
   % associate cellid (ix,iy,iz) to atom i
   cellid(i,:)=[ix,iy,iz];
end