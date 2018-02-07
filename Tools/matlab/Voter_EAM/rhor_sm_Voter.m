function f = rhor_sm_Voter(r, r_c, beta)

m = 20;
drhodr_rc = drhodr_Voter(r_c, beta);
f = rhor_Voter(r,beta) - rhor_Voter(r_c,beta) + (r_c/m)*(1-(r/r_c).^m)*drhodr_rc;

end