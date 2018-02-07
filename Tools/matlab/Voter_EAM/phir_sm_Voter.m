function f = phir_sm_Voter(r, r_c, alpha, r0, E1)

m = 20;
dphidr_rc = dphidr_Voter(r_c, alpha, r0, E1);
f = phir_Voter(r,alpha,r0,E1) - phir_Voter(r_c,alpha,r0,E1) + (r_c/m)*(1-(r/r_c).^m)*dphidr_rc;

end