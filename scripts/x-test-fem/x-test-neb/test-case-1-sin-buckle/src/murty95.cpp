/*
  murty95.h
  by Keonwook Kang kwkang75@yonsei.ac.kr
  Last Modified : 

  FUNCTION  :  
  Reference:   M. V. Ramana Murty and Harry A. Atwater
               "Empirical interatomic potential for Si-H interactions"
			   Phys. Rev. B, 51, 4889 (1995)
			   
			   Si-Si interaction is based on 
               J. Tersoff, Empirical interatomic potential for silicon 
               with improved elastic properties
               Phys. Rev. B, 38, 9902 (1988), Si(C)

			   H-H interaction is based on
               D. W. Brenner, "Empirical potential for hydrogen for use
			   in simulating the chemical vapor deposition of diamond films"
               Phys. Rev. B, 42, 9458 (1990)			   
			   
  Note:        Implemented by modifying tersoff.cpp
               The cut-off fuction f_c in eqn.(4) is identical expression 
               with the cut-off function in Tersoff PRB 39 5566(1989).
               Just, R_{ij} = R-D and S_{ij} = R+D or
                     R = (R_{ij} + S_{ij})/2, D = - (R_{ij} - S_{ij})/2
*/

#include "murty95.h"

inline double spline(double spline_coeff[4], double func[], int ind, double x)
{
    /* cubic spline based on Lagrange polynomials 
  
      f(x) = a_i/(6*h_i)*(x_{i+1} - x)^3 + a_{i+1}/(6*h_i)*(x-x_i)^3
             + (y_i/h_i-a_i*h_i/6)*(x_{i+1}-x)
             + (y_{i+1}/h_i-a_{i+1}*h_i/6)*(x-x_i) 
      where h_i = x_{i+1}-x_i
    */
    int xi, xip1;
    double a, b, fi, fip1, f;
    double pp, pp2, pp3;
    double qq, qq2, qq3;
  
    xi = ind; xip1 = ind+1; 
    a = spline_coeff[ind];
    b = spline_coeff[ind+1];
    fi = func[ind]; fip1 = func[ind+1];
    pp = x-xi; pp2 = pp*pp; pp3=pp2*pp;
    qq = xip1-x; qq2=qq*qq; qq3=qq2*qq;
    f=a/6.0*qq3+b/6.0*pp3+(fi-a/6.0)*qq+(fip1-b/6.0)*pp;
    return f;
}

void compute_spline_coeff(int ngrd, double func[],double h, double spline_coeff[])
{
    int i,ind,n=ngrd-2;
    if (n<1) ERROR("Increase ngrd in murty95.h"); 
	double dia[n], adia[n-1], bdia[n-1], B[n], x[n];

    /* Thomas algorithm to solve Ax=B
       A = [2h  h  0  0 ... 
            h  2h  h  0 ...
            0   h 2h  h ...
                    ...
                    h 2h ] */ 
    for (i=0;i<n;i++) 
    {
        dia[i]=2.0*h; B[i] = 6.0*(func[i+2]-func[i])/h;
    }
    for (i=0;i<(n-1);i++) 
    {
        adia[i]=h; // above diagonal
        bdia[i]=h; // below diagonal
    }
    adia[0]=adia[0]/dia[0];
    B[0]=B[0]/dia[0];
	for (i=1;i<(n-1);i++)
    {
        adia[i]=adia[i]/(dia[i]-bdia[i-1]*adia[i-1]);
        B[i]=(B[i]-bdia[i-1]*B[i-1])/(dia[i]-bdia[i-1]*adia[i-1]);
    }
    B[n-1]=(B[n-1]-bdia[n-2]*B[n-2])/(dia[n-1]-bdia[n-2]*adia[n-2]);
    x[n-1]=B[n-1];
    for (i=n-2;i>=0;i--)
        x[i]=B[i]-adia[i]*x[i+1];
    
    /* natural cubic spline : 
             the 2nd derivatives at the endpoints are zero */
    spline_coeff[0] = 0.0;
    for(ind=1;ind<ngrd-1;ind++)
    {
        spline_coeff[ind] = x[ind-1];
    }
    spline_coeff[ngrd-1] = 0.0;
}

void Murty95Frame::initvars()
{
#ifdef _MA_ORIG
    /* Fitting parameters for Si-Si, PRB 51 4889 (1995) */
    A_SiSi             = 1.8308e3;          // A (eV)
    B_SiSi             = 4.7118e2;          // B (eV)
    Lam_SiSi           = 2.4799;            // lambda_1 (1/\AA)
    Mu_SiSi            = 1.7322;            // lambda_2 (1/\AA)
    Alpha_SiSi         = 5.1975;            // 
    Beta_SiSi          = 3.0;               // beta
    R0_SiSi            = 2.35;              // equilibrium distance to 1nn (\AA)
    R_SiSi             = 2.85;              // R (\AA)
    D_SiSi             = 0.15;              // D (\AA)
    cTf_SiSi           = 0.0;               // c
    dTf_SiSi           = 0.160;             // d
    hTf_SiSi           =-5.9826e-1;         // h
    eta_SiSi           = 7.8734e-1;         
    delta_SiSi         = 0.635;
    F1_SiSi            = 1.0;
    F2_SiSi            = 1.0;
	
    /* Fitting parameters for H-H, PRB 51 4889 (1995)  */
    A_HH               = 80.07;
    B_HH               = 31.38;
    Lam_HH             = 4.2075;
    Mu_HH              = 1.7956;
    Alpha_HH           = 3.0;             
    Beta_HH            = 1.0;               // beta
    R0_HH              = 0.74;    
    R_HH               = 1.40;    
    D_HH               = 0.30;
    cTf_HH             = 4.0;               // c
    dTf_HH             = 0.0;               // d
    hTf_HH             = 0.0;               // h
    eta_HH             = 1.0;
    delta_HH           = 0.80469;
    F1_HH              = 1.0;
    F2_HH              = 1.0;
#else
#if 1
    /* Fitting parameters for Si-Si, 
       M.S. Kim and B.C Lee EPL 102 (2013) 33002 */
    A_SiSi             = 2.0119e3;          // A (eV)
    B_SiSi             = 7.1682e1;          // B (eV)
    Lam_SiSi           = 3.0145;            // lambda_1 (1/\AA)
    Mu_SiSi            = 1.1247;            // lambda_2 (1/\AA)
    Alpha_SiSi         = 3.27;              // 
    Beta_SiSi          = 2.0;               // beta
    R0_SiSi            = 2.37;              // equilibrium distance to 1nn (\AA)
    R_SiSi             = 3.33;              // R (\AA)
    D_SiSi             = 0.15;              // D (\AA)
    cTf_SiSi           = 0.01123;           // c
    dTf_SiSi           = 0.18220;           // d
    hTf_SiSi           =-3.3333e-1;         // h
    eta_SiSi           = 5.2772e-1;
    delta_SiSi         = 1.0;
    F1_SiSi            = 1.0;
    F2_SiSi            = 1.0;
#else
    /* new fitting Nov 28, 2014 */
    A_SiSi             = 8025.195302116495; // A (eV)
    B_SiSi             = 40.701115782280240;// B (eV)
    Lam_SiSi           = 3.909594726562502; // lambda_1 (1/\AA)
    Mu_SiSi            = 0.867186853088723; // lambda_2 (1/\AA)
    Alpha_SiSi         = 0.001;             // 
    Beta_SiSi          = 2.0;               // beta
    R0_SiSi            = 2.37;              // equilibrium distance to 1nn (\AA)
    R_SiSi             = 3.33;              // R (\AA)
    D_SiSi             = 0.15;              // D (\AA)
    cTf_SiSi           = 0.233029307210524; // c
    dTf_SiSi           = 1.072002775222063; // d
    hTf_SiSi           =-0.448928185468685; // h
    eta_SiSi           = 1.08671875;
    delta_SiSi         = 0.761787281185388;
    F1_SiSi            = 1.0;
    F2_SiSi            = 1.0;
#endif
    /* Fitting parameters for H-H,        
       M.S. Kim and B.C. Lee EPL 102 (2013) 33002 */
    A_HH               = 66.8276;
    B_HH               = 21.5038;
    Lam_HH             = 4.4703;
    Mu_HH              = 1.1616;
    Alpha_HH           = 3.0;
    Beta_HH            = 1.0;               // beta
    R0_HH              = 0.75;
    R_HH               = 3.50;
    D_HH               = 0.50;
    cTf_HH             = 0.9481;            // c
    dTf_HH             = 0.0;               // d
    hTf_HH             = 0.0;               // h
    eta_HH             = 1.6183;
    delta_HH           = 0.8684;
    F1_HH              = 1.0;
    F2_HH              = 1.0;

#endif

    /* Fitting parameters for Si-H, PRB 51 4889 (1995)  */
    A_SiH              = 323.54;
    B_SiH              = 84.18;
    Lam_SiH            = 2.9595;
    Mu_SiH             = 1.6158;
    Alpha_SiH          = 4.0;             Alpha_SiH2 = 0.00;
    Beta_SiH           = 3.0;             Beta_SiH2 = 0.00;
    R0_SiH             = 1.475;    
    R_SiH              = 1.85;    
    D_SiH              = 0.15;
    cTf_SiH            = 0.0216;          cTf_SiH2 = 0.70;
    dTf_SiH            = 0.27;            dTf_SiH2 = 1.00;
    eta_SiH            = 1.0;
    delta_SiH          = 0.80469;
    hTf_SiH[0]         = -0.040;          hTf_SiH2 =-1.00;
    hTf_SiH[1]         = -0.040;
    hTf_SiH[2]         = -0.040;
    hTf_SiH[3]         = -0.276;
    hTf_SiH[4]         = -0.47;
    F1_SiH[0]          = 0.0;
    F1_SiH[1]          = 1.005;
    F1_SiH[2]          = 1.109;
    F1_SiH[3]          = 0.953;
    F1_SiH[4]          = 1.0;
    F2_SiH[0]          = 0.0;
    F2_SiH[1]          = 0.930;
    F2_SiH[2]          = 1.035;
    F2_SiH[3]          = 0.934;
    F2_SiH[4]          = 1.0;
	
    acut=max(R_SiSi+D_SiSi,R_HH+D_HH);
	acut=max(acut,R_SiH+D_SiH)*1.0;
    _RLIST=acut*1.1;
    _SKIN=_RLIST-acut;
    acutsq=acut*acut;

    strcpy(incnfile,"../si.cn");
    DUMP(HIG"Murty95Frame initvars"NOR);
    MDFrame::initvars();
}



void Murty95Frame::murty()
{
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    /* Multi process function */
    int i,j,k,ipt,jpt,kpt,neg,n1nn=0;
	int species_ipt,species_jpt,species_kpt;
    Vector3 sij,rij,f3i,fZ,fTheta;
    double fC, fA, fR, dfCdr, dfAdr, dfRdr;
    double Arg1, Arg2, temp, Z, CosTheta, dCosTheta, gTheta, gR, dgRdr, Zeta, Zeta1, bij;
    double r2ij, rrij, r, neff1nn;
    double V, dVdr, dVdZ;
    double rr[NNM], Re[NNM], Neff1NN[_NP];
    //double Asave[NNM], Bsave[NNM], Musave[NNM], Lamsave[NNM];
    //double Rsave[NNM], Ssave[NNM];
    double fCsave[NNM], dfCdrsave[NNM], fAsave[NNM], dfAdrsave[NNM];
    double fRsave[NNM], dfRdrsave[NNM], dZdCosTheta[NNM];
    double dZdrij, dZdrik[NNM];
    Vector3 rg[NNM], rt[NNM];
    int list2[NNM];
    int n0, n1;

	DUMP("Murty & Atwater 95");
    refreshneighborlist();

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    n0=0;
    n1=_NP;

    compute_spline_coeff(NGRD,F1_SiH,1.0,F1_csp);
    compute_spline_coeff(NGRD,F2_SiH,1.0,F2_csp);
    compute_spline_coeff(NGRD,hTf_SiH,1.0,H_csp);

    memset(Neff1NN,0,sizeof(double)*_NP);
    /* 1st iteration to calculate coordination of Si atoms */
    for(ipt=n0;ipt<n1;ipt++)
    {
        if(fixed[ipt]==-1) continue; /* ignore atom if -1 */
        species_ipt=species[ipt];
        //if(species_ipt==1) continue; /* ignore if atom is H */
		
        neff1nn = 0; 
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
            
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            if(r2ij>acutsq) continue;
            rrij=sqrt(r2ij);

            r = rrij;
            species_jpt=species[jpt];
            /* set potential parameters */
            if (species_ipt==0 && species_jpt==0)
            { // Si-Si interaction
                R = R_SiSi; D = D_SiSi;
		    }
			else if (species_ipt==1 && species_jpt==1) 
			{ // H-H interaction
				R = R_HH; D = D_HH;
			}
			else if ((species_ipt==0&&species_jpt==1)
			         || (species_ipt==1&&species_jpt==0))
			{ // Si-H interaction
                R = R_SiH; D = D_SiH;
			}
			
            if (r>=(R+D)) continue;
            if(r<=(R-D))
                fC = 1;
            else if (r<(R+D))
            {
                Arg1 = M_PI*(r-R)/(2*D);
                Arg2 = 3.0*Arg1;
                fC = 0.5-9.0/16.0*sin(Arg1)-1/16.0*sin(Arg2);
            }
            else
                fC = 0;
            neff1nn+=fC;
        }
        Neff1NN[ipt]=neff1nn;
        //INFO("neff1nn="<<neff1nn);
    }    
    
    /* 2nd iteration to calculate individual energy and force */
    for(ipt=n0;ipt<n1;ipt++)
    {
        if(fixed[ipt]==-1) continue; /* ignore atom if -1 */
        species_ipt=species[ipt];

        memset(fCsave,0,sizeof(double)*NNM);
        memset(dfCdrsave,0,sizeof(double)*NNM);
        memset(fAsave,0,sizeof(double)*NNM);
        memset(dfAdrsave,0,sizeof(double)*NNM);
        memset(fRsave,0,sizeof(double)*NNM);
        memset(dfRdrsave,0,sizeof(double)*NNM);
        memset(dZdrik,0,sizeof(double)*NNM);
        memset(dZdCosTheta,0,sizeof(double)*NNM);
        memset(rr,0,sizeof(double)*NNM);
        memset(Re,0,sizeof(double)*NNM);
        memset(list2,0,sizeof(int)*NNM);
        for(j=0;j<NNM;j++) 
        {
            rg[j].clear(); rt[j].clear();
        }
        neg = 0;        
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            if(fixed[jpt]==-1) continue; /* ignore atom if -1 */
            
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            if(r2ij>acutsq) continue;
            rrij=sqrt(r2ij);

            r = rrij;
			species_jpt=species[jpt];
            /* set potential parameters */
			if (species_ipt==0 && species_jpt==0) 
			{ // Si-Si interaction
                A = A_SiSi; B = B_SiSi;
				Lam = Lam_SiSi; Mu = Mu_SiSi;
				R0 = R0_SiSi; R = R_SiSi; D = D_SiSi;
		    }
			else if (species_ipt==1 && species_jpt==1) 
			{ // H-H interaction
                A = A_HH; B = B_HH;
				Lam = Lam_HH; Mu = Mu_HH;
				R0 = R0_HH; R = R_HH; D = D_HH;
			}
			else if ((species_ipt==0&&species_jpt==1)
			         || (species_ipt==1&&species_jpt==0))
			{
                A = A_SiH; B = B_SiH;
				Lam = Lam_SiH; Mu = Mu_SiH;
				R0 = R0_SiH; R = R_SiH; D = D_SiH;
			}
			
            if (r>=(R+D)) continue;
            if(r<=(R-D))
            {
                fC = 1;
                dfCdr = 0;
            }
            else if (r<(R+D))
            {
                Arg1 = M_PI*(r-R)/(2.0*D);
                Arg2 = 3.0*Arg1;
                fC = 0.5-9.0/16.0*sin(Arg1)-1.0/16.0*sin(Arg2);
                dfCdr = -3.0/16.0*(3.0*cos(Arg1)+cos(Arg2))*M_PI/(2.0*D);
            }
            else
            {
                fC = 0;
                dfCdr = 0;
            }
			
            /* Repulsion, fR = A*F1*exp(-Lam*r) */
            fR = A*exp(-Lam*r);
            dfRdr = -Lam*fR;
            /* Attraction, fA = -B*F2*exp(-Mu*r) */
            fA = -B*exp(-Mu*r);
            dfAdr = -Mu*fA;

            /* save all neighbors of i */
            rg[neg]=rij;                /* displacement vector */
            rt[neg]=rij; rt[neg]/=rrij; /* unit vector */
            rr[neg]=rrij;               /* distance */
            list2[neg]=jpt;             /* neighbor index */
			Re[neg]=R0;                 /* equilibrium distance between i and j */
        
            fCsave[neg]=fC;
            dfCdrsave[neg]=dfCdr;
            fAsave[neg]=fA;
            dfAdrsave[neg]=dfAdr;
            fRsave[neg]=fR;
            dfRdrsave[neg]=dfRdr;
            
            neg++;
        }
        
        //INFO_Printf("neg=%d\n",neg);
        /* going through the neighbors of atom i again */
        for(j=0;j<neg;j++)
        {
            jpt=list2[j];
            if(fixed[jpt]==-1) continue; /* ignore atom if -1 */

			species_jpt=species[jpt];
            /* set potential parameters */
			if (species_ipt==0 && species_jpt==0 )
			{ // Si-Si interaction
                A = A_SiSi; B = B_SiSi;
				Lam = Lam_SiSi; Mu = Mu_SiSi;
                Alpha = Alpha_SiSi; Beta = Beta_SiSi;
                cTf = cTf_SiSi; dTf = dTf_SiSi; hTf = hTf_SiSi; 
				R0 = R0_SiSi; R = R_SiSi; D = D_SiSi;
				eta = eta_SiSi; delta = delta_SiSi;
				F1 = F1_SiSi; F2 = F2_SiSi; 
		    }
			else if (species_ipt==1 && species_jpt==1 ) 
			{ // H-H interaction
                A = A_HH; B = B_HH;
				Lam = Lam_HH; Mu = Mu_HH;
				Alpha = Alpha_HH; Beta = Beta_HH;
				R0 = R0_HH; cTf = cTf_HH; dTf = dTf_HH;
				hTf = hTf_HH; R = R_HH; D = D_HH;
				eta = eta_HH; delta = delta_HH;
				F1 = F1_HH; F2 = F2_HH; 
			}
			else if ((species_ipt==0&&species_jpt==1) 
			         ||(species_ipt==1&&species_jpt==0))
			{ // Si-H interaction
                A = A_SiH; B = B_SiH;
				Lam = Lam_SiH; Mu = Mu_SiH;
				Alpha = Alpha_SiH; Beta = Beta_SiH;
				R0 = R0_SiH; R = R_SiH; D = D_SiH;
                cTf = cTf_SiH; dTf = dTf_SiH;
				eta = eta_SiH; delta = delta_SiH;

                /* Assign the effective coordination of Si */
                if (species_ipt==0) neff1nn=Neff1NN[ipt];
                if (species_ipt==1) neff1nn=Neff1NN[jpt];
                if (neff1nn>=4.0)
                {
                    hTf = hTf_SiH[NGRD-1]; F1 = F1_SiH[NGRD-1]; F2 = F2_SiH[NGRD-1];
                }
                else if (neff1nn<4.0) 
                {
                    n1nn = (int)neff1nn;
                    F1 = spline(F1_csp,F1_SiH,n1nn,neff1nn);
                    F2 = spline(F2_csp,F2_SiH,n1nn,neff1nn);
                    /* hTf will be obtained in the three-body loop. */
                }
            }
            //INFO("neff1nn = "<<neff1nn);
            //INFO("F1 = "<<F1);
            //INFO("F2 = "<<F2);
            //INFO("H = "<<hTf);
            
            dZdrij=0.0;
            Z = 0.0; /* bond order for attractive term */
            /* three body loop */
            for(k=0;k<neg;k++)
            {
                if(k==j) continue; /* k and j are different neighbors */
                kpt=list2[k];
                /* if fixed == -1, simply ignore this atom */
                if(fixed[kpt]==-1) continue; /* ignore atom if -1 */
                species_kpt=species[kpt];
				/* Calculate F1(N), F2(N) and hTf using cubic spline */
				
                /* consider triplet (i,j,k)
                       j   k
                         V   
                         i                */
                if(species_ipt==0&&species_jpt==0&&species_kpt==0)
                { //(i,j,k) = (Si,Si,Si) I
                    Alpha = Alpha_SiSi; Beta = Beta_SiSi;
                    cTf = cTf_SiSi; dTf = dTf_SiSi; hTf = hTf_SiSi; 
                }
                if(species_ipt==0&&species_jpt==0&&species_kpt==1)
                { //(Si,Si,H) IIa
                    Alpha = Alpha_SiH; Beta = Beta_SiH;
                    cTf = cTf_SiH; dTf = dTf_SiH;
                    hTf = spline(H_csp,hTf_SiH,(int)Neff1NN[ipt],Neff1NN[ipt]);
                }
                if(species_ipt==0&&species_jpt==1&&species_kpt==0)
                { //(Si,H,Si) IIa
                    Alpha = Alpha_SiH; Beta = Beta_SiH;
                    cTf = cTf_SiH; dTf = dTf_SiH;
                    hTf = spline(H_csp,hTf_SiH,(int)Neff1NN[ipt],Neff1NN[ipt]);
                }
                if(species_ipt==0&&species_jpt==1&&species_kpt==1)
                { //(Si,H,H) IIa
                    Alpha = Alpha_SiH; Beta = Beta_SiH;
                    cTf = cTf_SiH; dTf = dTf_SiH;
                    hTf = spline(H_csp,hTf_SiH,(int)Neff1NN[ipt],Neff1NN[ipt]);
                }
                if(species_ipt==1&&species_jpt==0&&species_kpt==0)
                { //(H,Si,Si) IIb
                    Alpha = Alpha_SiH2; Beta = Beta_SiH2;
                    cTf = cTf_SiH2; dTf = dTf_SiH2; hTf = hTf_SiH2;
                }
                if(species_ipt==1&&species_jpt==0&&species_kpt==1)
                { //(H,Si,H) III
                    Alpha = Alpha_HH; Beta = Beta_HH;
                    cTf = cTf_HH; dTf = dTf_HH; hTf = hTf_HH;
                }
                if(species_ipt==1&&species_jpt==1&&species_kpt==0)
                { //(H,H,Si) III
                    Alpha = Alpha_HH; Beta = Beta_HH;
                    cTf = cTf_HH; dTf = dTf_HH; hTf = hTf_HH;
                }
                if(species_ipt==1&&species_jpt==1&&species_kpt==1)
                { //(H,H,H) III
                    Alpha = Alpha_HH; Beta = Beta_HH;
                    cTf = cTf_HH; dTf = dTf_HH; hTf = hTf_HH;
                }
                //INFO("H = "<<hTf);
                
                
                /* cos(theta_ijk): angle between bond ij and bond ik */
                CosTheta=rt[j]*rt[k];
                /* h-cos(theta_ijk) */
                dCosTheta = (hTf-CosTheta);
                /* g(theta)=c+d*(h-cos(theta_ijk))^2 */
                gTheta = cTf+dTf*dCosTheta*dCosTheta;
                /* g(r) = exp[alpha*{(r_ij-R0_ij)-(r_ik-R0_ik)}^beta] */
                gR = exp(Alpha*pow((rr[j]-Re[j]-rr[k]+Re[k]),Beta));
                Z += fCsave[k]*gTheta*gR;
                /* d(gR)/d(r_ik) */
                //dgRdr = gR*(-1.0)*Alpha*Beta*pow((rr[j]-Re[j]-rr[k]+Re[k]),Beta-1);
                /* d(gR)/d(r_ij) */
                dgRdr = gR*Alpha*Beta*pow((rr[j]-Re[j]-rr[k]+Re[k]),Beta-1);
                dZdrij+=fCsave[k]*gTheta*dgRdr;
                dZdrik[k] = dfCdrsave[k]*gTheta*gR-fCsave[k]*gTheta*dgRdr;
                dZdCosTheta[k] = (fCsave[k]*gR)*(-2.0)*dTf*(hTf-CosTheta);
                //INFO("g(theta)="<<gTheta<<", g(r)="<<gR);
            }
            
            fC = fCsave[j];
            dfCdr = dfCdrsave[j];
            fA = fAsave[j];
            dfAdr = dfAdrsave[j];
            fR = fRsave[j];
            dfRdr = dfRdrsave[j];
            
            Zeta = pow(Z,eta);
            Zeta1 = pow(Z,eta-1);
            if(neg==1) Zeta1=0.0;

            bij = pow((1.0+Zeta),(-1.0)*delta);
            V = 0.5*fC*(F1*fR+bij*F2*fA);
            //INFO("Z="<<Z<<", bij = "<<bij);
            //INFO("V="<<V<<" (eV)");
    
#if 0
// Debug
            FILE *fp;
            if(ipt==0)
            {
            INFO("write energy.dat");
            if(j==0) fp=fopen("energy.dat","w");
            else     fp=fopen("energy.dat","a");
            if(fp==NULL) FATAL("open file failure");
            fprintf(fp, "%d %d %12.10e %12.10e %12.10e %12.10e %12.10e %12.10e\n", ipt,jpt,fC, F1, F2, fR,fA,bij);
            fclose(fp);
            }
#endif

            _EPOT += V;
            _EPOT_IND[ipt] += V;

            dVdZ = 0.5*fC*F2*fA*(-bij*(delta*eta*Zeta1)/(1+Zeta));
            dVdr = 0.5*dfCdr*(F1*fR+bij*F2*fA) + 0.5*fC*(F1*dfRdr+bij*F2*dfAdr) + dVdZ*dZdrij;

            f3i=rt[j]*dVdr;
            _F[ipt]+=f3i;
            _F[jpt]-=f3i;
            _VIRIAL.addnvv(-1.,f3i,rg[j]);
            if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f3i,rg[j],-1.0);

//            INFO_Printf(" bij = %e V = %e dVdr = %e f3i=(%e,%e,%e)\n",bij,V,dVdr,f3i.x,f3i.y,f3i.z);            
//            INFO_Printf(" bij = %e rt[j]=(%e,%e,%e)\n",bij,rt[j].x,rt[j].y,rt[j].z);            

            
            /* three body loop */
            for(k=0;k<neg;k++)
            {
                if(k==j) continue; /* k and j are different neighbors */
                kpt=list2[k];
                if(fixed[kpt]==-1) continue; /* ignore atom if -1 */

                /* Forces due to derivatives with respect to Z, then rik */
                temp = dVdZ*dZdrik[k];

                fZ = rt[k]*temp;
                _F[ipt]+=fZ;
                _F[kpt]-=fZ;
                _VIRIAL.addnvv(-1.,fZ,rg[k]);
                if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],fZ,rg[k],-1.0);

                /* cos(theta_ijk): angle between bond ij and bond ik */
                CosTheta=rt[j]*rt[k];
                
                /* Forces due to derivatives with respect to cosTheta */
                temp = dVdZ*dZdCosTheta[k];

                fTheta = (rt[j]*CosTheta-rt[k])*(-temp/rr[j]);
                _F[ipt]+=fTheta;
                _F[jpt]-=fTheta;
                _VIRIAL.addnvv(-1.,fTheta,rg[j]);
                if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],fTheta,rg[j],-1.0);


                fTheta = (rt[k]*CosTheta-rt[j])*(-temp/rr[k]);
                _F[ipt]+=fTheta;
                _F[kpt]-=fTheta;
                _VIRIAL.addnvv(-1.,fTheta,rg[k]);
                if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],fTheta,rg[k],-1.0);

            }
        }          
    }
}
    

void Murty95Frame::potential()
{
    murty();
}

#ifdef _TEST

/* Main Program Begins */
class Murty95Frame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

