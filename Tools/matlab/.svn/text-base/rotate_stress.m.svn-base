% rotate_stress.m
%
% rotate 2nd order tensor given rotation matrix T
%  Sig'_ij = Q_ip Q_jp Sig_pq (sum over p,q)
%   where Q_ip = dot(e'_i, e_p)
%
% Wei Cai
% ME340 Elasticity of Microscopic Structures, Stanford University, Winter 2006
%

% stress in old coordinate system
sig = [1 0 0
       0 0 0
       0 0 0 ];
   
% we want to find stress in the new coordinate system
sigp = zeros(3,3);

% new coordinate system
e1p = [1 1 -2];
e2p = [-1 1 0];
e3p = [1 1  1];

e1p = e1p/norm(e1p);
e2p = e2p/norm(e2p);
e3p = e3p/norm(e3p);

% original coordinate system
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

Q = [dot(e1p,e1) dot(e1p,e2) dot(e1p,e3)
     dot(e2p,e1) dot(e2p,e2) dot(e2p,e3)
     dot(e3p,e1) dot(e3p,e2) dot(e3p,e3)];
 
for i=1:3,
 for j=1:3,
  for p=1:3,
   for q=1:3,
       sigp(i,j) = sigp(i,j) + Q(i,p)*Q(j,q)*sig(p,q);
   end
  end
 end
end

sig,sigp