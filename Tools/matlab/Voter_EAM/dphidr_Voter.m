function f = dphidr_Voter(r, alpha, r0, E1)

f = 2*E1*alpha*( exp(2*alpha*(r-r0)) - exp(alpha*(r-r0)) );

end