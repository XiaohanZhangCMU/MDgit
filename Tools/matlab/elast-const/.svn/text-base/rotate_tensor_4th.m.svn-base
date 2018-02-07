% rotate_elast_stiff.m
%
% rotate 4th order tensor given rotation matrix T
%  Cp_ijkl = T_mi T_nj T_sk T_tl C_mnst (sum over m,n,s,t)
%
% Wei Cai
% ME340 Elasticity of Microscopic Structures, Stanford University, Winter 2006
%
function Cp = rotate_elast_stiff(C,T)

Cp = zeros(size(C));

for i=1:3,
 for j=1:3,
  for k=1:3,
   for l=1:3,
       for m=1:3,
        for n=1:3,
         for s=1:3,
          for t=1:3,
              Cp(i,j,k,l)=Cp(i,j,k,l)+T(m,i)*T(n,j)*T(s,k)*T(t,l)*C(m,n,s,t);
          end
         end
        end
       end
   end
  end
 end
end
