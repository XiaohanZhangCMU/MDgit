/*
  ewald.cpp
  by Hark Lee and Wei Cai, hblee@stanford.edu, caiwei@stanford.edu
  Last Modified : Thu Oct 18 18:20:05 2007

  FUNCTION  :  Ewald Summation (classical Ewald and Particle Mesh Ewald)

  To Do:
      1.  CE Virial stress still incorrect (perhaps from Self energy contribution)
      2.  CG relax dies early when both Y and Oxygen vacancy exist
      3.  PME still needs to be included

  Note:
      This file is included by md.cpp.  This file cannot be compiled by itself.
*/

/********************************************************************/
/* Ewald Summation for Coulomb Interaction */
/********************************************************************/

//###################################################################
//  CE Functions
//###################################################################

//-------------------------------------------------------------------
// Warning
// Note that if you want to switch conj_fixbox (from 1 to 0) or (from 
// 0 to 1) in the middle of simulation, CE should be reinitialized.
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// Classical Ewald
//-------------------------------------------------------------------
void MDFrame::CE()
{
    int i;
    CE_Real();
    CE_Rec();

    _EPOT_Ewald = _EPOT_Ewald_Real + _EPOT_Ewald_Rec + _EPOT_Ewald_Self;
    _VIRIAL_Ewald = _VIRIAL_Ewald_Real + _VIRIAL_Ewald_Rec;

    for(i=0;i<_NP;i++)
      {
        _F_Ewald[i] = _F_Ewald_Real[i] + _F_Ewald_Rec[i];
        _EPOT_IND_Ewald[i] = _EPOT_IND_Ewald_Real[i] + _EPOT_IND_Ewald_Rec[i];
      }
}

//-------------------------------------------------------------------
// Compute Real Part of Ewald Summation
//-------------------------------------------------------------------
void MDFrame::CE_Real()
{
    int i,j,ipt,jpt,isp,jsp; // integers for atom counting, atom indices, atom species
    double UC,r,r2,QQ,EXP;   // U = potential | r = radius | r2 = r^2 | ri6 = r^-6
    Vector3 sij, rij, fijC; // sij : scaled coordinates | rij : coordinates | fij : forces

    DUMP(HIG"Coulomb Potential (real)"NOR);

    // start time of real part
    CE_sREAL=clock();
    
    // refresh neighbor list
    refreshneighborlist();
    
    // initialize force, potential and virials
    _EPOT_Ewald_Real = 0.0;
    _VIRIAL_Ewald_Real.clear();

    for(i=0;i<_NP;i++)
    {
        _F_Ewald_Real[i].clear();
        _EPOT_IND_Ewald_Real[i]=0.0;
    }
    
    //-------------------------------------------------------------------
    //       Calculating Short Range Potential and Ewald Real Part
    //-------------------------------------------------------------------
    for(ipt=0;ipt<_NP;ipt++) // for all particles
    {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[ipt]==-1) continue;
        for(j=0;j<nn[ipt];j++) // for all neighbors
	{
            jpt=nindex[ipt][j]; // get the atom index of j-th neighbor
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            if(ipt>=jpt) continue; // only if i<j -> avoid redundancy
            sij=_SR[jpt]-_SR[ipt]; // scaled interatomics distance
            sij.subint(); // apply periodic boundary condition (PBC)
            rij=_H*sij; // get the interatomic distance vector (with PBC)
            r2=rij.norm2(); // get the square of interatomic distance
            r=sqrt(r2); // interatomic distance
            
            if(r<=Ewald_Rcut) // if the interatomic distance is shorter than the cut-off radius
	    {
                // get the species indices
                isp=species[ipt]; // species of ipt-atom
                jsp=species[jpt]; // species of jpt-atom
                QQ=P_CLMB * _ATOMCHARGE[isp]*_ATOMCHARGE[jsp]; // get the charge term
                
                // calculate potential
                UC=(erfc(Ewald_Alpha*r))/r*QQ;
		
                // calculate force
                EXP=(2.0/M_SQRT_PI)*Ewald_Alpha*exp(-SQR(Ewald_Alpha*r))*QQ;
                fijC=rij*((UC+EXP)/r2);
                
                // add pair terms to the accumulative sum
                _F_Ewald_Real[ipt]-=fijC; _F_Ewald_Real[jpt]+=fijC;

                _EPOT_Ewald_Real+=UC;

                _EPOT_IND_Ewald_Real[ipt]+=UC*0.5;
                _EPOT_IND_Ewald_Real[jpt]+=UC*0.5;
                _VIRIAL_Ewald_Real.addnvv(1.,fijC,rij);
	    }
	}
    }
    
    // end time of real part
    CE_eREAL=clock();
}

//-------------------------------------------------------------------
// Compute Reciprocal Part of Ewald Summation
//-------------------------------------------------------------------
void MDFrame::CE_Rec()
{    
    int i, ik, kal, kbl;
    int ka, kb, kc, np_real;
    double k2, kx, ky, kz;
    double pfac, tmp_H_L1, tmp_H_L2, tmp_H_L3;
    double UE, EXP;
    EComplex sincosk;

     // start time of reciprocal part
    CE_sREC=clock();
    
    // initialize force, potential and virials
   _EPOT_Ewald_Rec = 0.0;
    _VIRIAL_Ewald_Rec.clear();

    for(i=0;i<_NP;i++)
    {
        _F_Ewald_Rec[i].clear();
        _EPOT_IND_Ewald_Rec[i]=0.0;
    }
    
    // recalculate the system volume (in case the box is changing)
    if(conj_fixbox==0)
      {
        tmp_H_L1=sqrt(SQR(_H[0][0])+SQR(_H[1][0])+SQR(_H[2][0]));
        tmp_H_L2=sqrt(SQR(_H[0][1])+SQR(_H[1][1])+SQR(_H[2][1]));
        tmp_H_L3=sqrt(SQR(_H[0][2])+SQR(_H[1][2])+SQR(_H[2][2]));
        Ewald_cell_V=_H.det();        
      }
    else
      {
	// just to avoid warning
        tmp_H_L1=Ewald_H_L1;
        tmp_H_L2=Ewald_H_L2;
        tmp_H_L3=Ewald_H_L3;
      }

    // precalculate exp(iKx*rx), exp(iKy*ry), exp(iKz*rz)
    CE_fillSinCosTables();
    
    kal=-1; // an impossible value
    kbl=0;  // Just to avoid a warning
    for(ik=0;ik<CE_nKV;ik++) // for all K-vectors
      {
        // kx, ky, kz calculated in enumKV()
        ka=CE_KV[ik].a;
        kb=CE_KV[ik].b;
        kc=CE_KV[ik].c;
	
        if(conj_fixbox==0)
	  {
            kx=CE_KV[ik].kx*Ewald_H_L1/tmp_H_L1;
            ky=CE_KV[ik].ky*Ewald_H_L2/tmp_H_L2;
            kz=CE_KV[ik].kz*Ewald_H_L3/tmp_H_L3;
            k2=kx*kx+ky*ky+kz*kz;
	  }
        else
	  {
            kx=CE_KV[ik].kx;
            ky=CE_KV[ik].ky;
            kz=CE_KV[ik].kz;
            k2=CE_KV[ik].k2;
	  }
        
        //quick sckk [ Qj*exp(iKr) ] calculation
        if(kal!=ka || kbl!=kb)  // Re-calculate sc_axby
	  {
            kal=ka;
            kbl=kb;
            if(kb>=0) // if ky>=0
	      for(i=0;i<_NP;i++) // for all atoms
		EComplex::Mul(CE_scx[ka][i], CE_scy[kb][i], CE_sc_axby[i]);
            else      // if ky<0
	      for(i=0;i<_NP;i++) // for all atoms
		EComplex::MulConjugate(CE_scx[ka][i], CE_scy[-kb][i], CE_sc_axby[i]);
	  }
        // initialize structure factor 
        sincosk.Re=0.0;
        sincosk.Im=0.0;
        
        if(kc>=0) // if kz>=0
	  for(i=0;i<_NP;i++) 
            {
	      /* if fixed == -1, simply ignore this atom */
	      if(fixed[i]==-1) continue;
	      
	      EComplex::Mul(CE_sc_axby[i], CE_scz[kc][i], CE_sckk[i]);
	      CE_sckk[i]*=P_SQRT_CLMB*_ATOMCHARGE[species[i]];
	      sincosk+=CE_sckk[i];
            }
        else      // if kz<0
	  for(i=0;i<_NP;i++)
            {
	      /* if fixed == -1, simply ignore this atom */
	      if(fixed[i]==-1) continue;
	      
	      EComplex::MulConjugate(CE_sc_axby[i], CE_scz[-kc][i], CE_sckk[i]);
	      CE_sckk[i]*=P_SQRT_CLMB*_ATOMCHARGE[species[i]];
	      sincosk+=CE_sckk[i];
            }
	
        // calculate the reciprocal energy pre-factor
        EXP=((8.0*M_PI)/Ewald_cell_V)*exp(-k2/(4.0*SQR(Ewald_Alpha)))/k2;
	
        // calculate the reciprocal force
        for(i=0;i<_NP;i++)
	  {
            /* if fixed == -1, simply ignore this atom */
            if(fixed[i]==-1) continue;
            
            double ad=EXP*(CE_sckk[i].Im*sincosk.Re-CE_sckk[i].Re*sincosk.Im);
            _F_Ewald_Rec[i].x+=ad*((double)kx);
            _F_Ewald_Rec[i].y+=ad*((double)ky);
            _F_Ewald_Rec[i].z+=ad*((double)kz);
	  }
        pfac=1.0/k2+1/(4.0*SQR(Ewald_Alpha));
        UE=EXP/2.0*sincosk.Norm2();
	
        // Potential Energy Calculation        
        _EPOT_Ewald_Rec+=UE;
	
        /* Temporary: distribute rec energy equally to all atoms */
        np_real = 0;
        for(i=0;i<_NP;i++)
	  if(fixed[i]!=-1) np_real++;
	
        for(i=0;i<_NP;i++)
	  {        
            /* if fixed == -1, simply ignore this atom */
            if(fixed[i]==-1) continue;
            
            _EPOT_IND_Ewald_Rec[i]=UE/np_real;
	  }
        
        // Virial Stress Calculation
        _VIRIAL_Ewald_Rec[0][0]+=(UE*(1.0-2.0*(kx*kx)*pfac)); //xx
        _VIRIAL_Ewald_Rec[1][1]+=(UE*(1.0-2.0*(ky*ky)*pfac)); //yy
        _VIRIAL_Ewald_Rec[2][2]+=(UE*(1.0-2.0*(kz*kz)*pfac)); //zz
        _VIRIAL_Ewald_Rec[0][1]+=(UE*(-2.0*(kx*ky)*pfac));  //xy
        _VIRIAL_Ewald_Rec[1][0]+=(UE*(-2.0*(kx*ky)*pfac));  //yx
        _VIRIAL_Ewald_Rec[0][2]+=(UE*(-2.0*(kx*kz)*pfac));  //xz
        _VIRIAL_Ewald_Rec[2][0]+=(UE*(-2.0*(kx*kz)*pfac));  //zx
        _VIRIAL_Ewald_Rec[1][2]+=(UE*(-2.0*(ky*kz)*pfac));  //yz
        _VIRIAL_Ewald_Rec[2][1]+=(UE*(-2.0*(ky*kz)*pfac));  //zy
      }
    //INFO("[CE] Reciprocal Energy Calculated by CE = "<<_EPOT_Ewald_Rec);
    
    // start time of reciprocal part
    CE_eREC=clock();
}

//-------------------------------------------------------------------
// Initialize Ewald  (which will call CE_init() or PME_init
//-------------------------------------------------------------------
void MDFrame::Ewald_init()
{
    Matrix33 hinv;
    double tmp1, tmp2, tmp3, tmp4, tmp_RLIST;
    int i;
    
    // get the system volume
    Ewald_cell_V=_H.det();

    // Calculate inverse matrix and length of lattice vectors
    hinv = _H.inv();
    Ewald_H_L1=sqrt(SQR(_H[0][0])+SQR(_H[1][0])+SQR(_H[2][0]));
    Ewald_H_L2=sqrt(SQR(_H[0][1])+SQR(_H[1][1])+SQR(_H[2][1]));
    Ewald_H_L3=sqrt(SQR(_H[0][2])+SQR(_H[1][2])+SQR(_H[2][2]));

    // Allocate real space force
    Realloc(_F_Ewald,     Vector3,_NP);
    Realloc(_F_Ewald_Real,Vector3,_NP);
    Realloc(_F_Ewald_Rec ,Vector3,_NP);

    Realloc(_EPOT_IND_Ewald,     double,_NP);
    Realloc(_EPOT_IND_Ewald_Real,double,_NP);
    Realloc(_EPOT_IND_Ewald_Rec ,double,_NP);
    
    bindvar("F_Ewald",_F_Ewald,DOUBLE);
    bindvar("F_Ewald_Real",_F_Ewald_Real,DOUBLE);
    bindvar("F_Ewald_Rec", _F_Ewald_Rec, DOUBLE);
    
    bindvar("EPOT_IND_Ewald",_EPOT_IND_Ewald,DOUBLE);
    bindvar("EPOT_IND_Ewald_Real",_EPOT_IND_Ewald_Real,DOUBLE);
    bindvar("EPOT_IND_Ewald_Rec", _EPOT_IND_Ewald_Rec, DOUBLE);
    
    //-------------------------------------------------------------------    
    // Detemine the Ewald accuracy
    //-------------------------------------------------------------------
    // Ewaldprecision=sqrt(-ln(error*Ewaldprecision^2))
    //-------------------------------------------------------------------
    // sqrt(-ln(10^-4)) = 3.03485425877029
    // sqrt(-ln(10^-5)) = 3.39307021220756
    // sqrt(-ln(10^-6)) = 3.71692218884984
    // sqrt(-ln(10^-7)) = 4.01473481701573
    // sqrt(-ln(10^-8)) = 4.29193205257869
    // sqrt(-ln(10^-9)) = 4.55228138815544
    // sqrt(-ln(10^-10))= 4.79852591218808
    //-------------------------------------------------------------------

    if(Ewald_CE_or_PME==0)
    { /* CE: classical Ewald */
        if( Ewald_option_Alpha==0 )
        {
            // determine Alpha (for CE only)
            if( (Ewald_time_ratio==0) || (Ewald_precision==0))
            {
                ERROR("Ewald_time_ratio or Ewald_precision is Not Specified -> Check Your Script File");
            }
            // get the optimal alpha
            Ewald_Alpha = M_SQRT_PI*cbrt(sqrt(Ewald_time_ratio*_NP)/Ewald_cell_V);
            Ewald_Rcut = Ewald_precision / Ewald_Alpha;
            
            if (_SKIN <= 0) _SKIN = 1.0;
            tmp_RLIST = Ewald_Rcut +_SKIN;
            tmp1=Ewald_H_L1/3.0;
            tmp2=Ewald_H_L2/3.0;
            tmp3=Ewald_H_L3/3.0;
            tmp4=tmp_RLIST;
            
            if(tmp4>tmp1){tmp4=tmp1;}
            if(tmp4>tmp2){tmp4=tmp2;}
            if(tmp4>tmp3){tmp4=tmp3;}
            
            if(tmp_RLIST>tmp4)
            {
                WARNING("PP_RC Computed from Optimum Alpha is Too Big");
                WARNING("Automatically Setting _RLIST Equal to 1/3 of Minimum Cell Length");
                Ewald_Rcut = (tmp4*0.95) - _SKIN;
            }
        }
    
        if(Ewald_option_Alpha==1)
        {
            if( (Ewald_Rcut==0) || (Ewald_precision==0))
            {
                ERROR("Ewald_Rcut or Ewald_precision is Not Specified -> Check Your Script File");
            }
            // get the controlled alpha
            Ewald_Alpha = Ewald_precision / Ewald_Rcut;
        }
    }
    else
    { /* PME: particle mesh Ewald */        
        // check if Alpha has non-zero value
        if((Ewald_option_Alpha==2)&&(Ewald_Alpha==0))
        {
            WARNING("Ewald Parameter Alpha is Not Specified -> Check Your Script File");
        }
    }

    
    // Update Verlet list cut-off
    _RLIST = Ewald_Rcut + _SKIN;
    
    // calculate sum of square charge
    Ewald_qSQRsum=0.0;
    for(i=0;i<_NP;i++)
    {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[i]==-1) continue;
        Ewald_qSQRsum+=SQR(P_SQRT_CLMB*_ATOMCHARGE[species[i]]);
    }

    // compute the self-interaction energy
    _EPOT_Ewald_Self = - Ewald_Alpha / M_SQRT_PI * Ewald_qSQRsum;

    if(Ewald_CE_or_PME == 0)
    {
        CE_init();
    }
    else
    {
        PME_init();
    }
}

//-------------------------------------------------------------------
// Initialize Classical Ewald 
//-------------------------------------------------------------------
void MDFrame::CE_init()
{
    // Allocating sc_axby and sckk
    Realloc(CE_sckk,EComplex,(_NP*2));
    CE_sc_axby=CE_sckk+_NP;

    // determine the cutoff radius
    CE_Kc = 2*Ewald_Alpha*Ewald_precision; // reciprocal-space
      
    // Calculate the K-vectors
    CE_enumKV();

    // Print Out Ewald Information
    INFO("------------------------------------------------------------------------------");
    INFO("------------------------------------------------------------------------------");
    INFO("              Classical Ewald Parameters [Long Range Interaction] ");
    INFO("------------------------------------------------------------------------------");
    INFO("    Ewald Fourier Part Cutoff Radius (Kc)         = "<<CE_Kc);
    INFO("    Ewald Gaussian Factor (Alpha)                 = "<<Ewald_Alpha);
    INFO("    Total Number of K-Vectors                     = "<<CE_nKV);
    INFO("    Max Number of K-Vectors                       = "<<CE_MKV);
    INFO("    Max K-Points in A-direction                   = "<<CE_Mka);
    INFO("    Max K-Points in B-direction                   = "<<CE_Mkb);
    INFO("    Max K-Points in C-direction                   = "<<CE_Mkc);
    INFO("------------------------------------------------------------------------------");
    INFO("------------------------------------------------------------------------------");
}

//-------------------------------------------------------------------
//  Clear Classical Ewald
//-------------------------------------------------------------------
void MDFrame::CE_clear()
{
    // free variables used in the Classical Ewald summation    
    if(CE_scx) { free(CE_scx[0]); free(CE_scx); }
    if(CE_sckk){ free(CE_sckk); free(CE_KV); }
}

//-------------------------------------------------------------------
// Filling Sin & Cos Table
//-------------------------------------------------------------------
// calculate exp(-iKx*rx), exp(-iKy*ry), exp(-iKz*rz) 
// for all atoms and all K-vectors, should be updated every time step
//-------------------------------------------------------------------
void MDFrame::CE_fillSinCosTables()
{
    int i,ki;
    EComplex *sc1;
    for(i=0;i<_NP;i++)
    {
        // when K-vector is zero vector (kx=0, ky=0, kz=0)
        CE_scx[0][i].Re=1.0; CE_scx[0][i].Im=0.0;
        CE_scy[0][i].Re=1.0; CE_scy[0][i].Im=0.0;
        CE_scz[0][i].Re=1.0; CE_scz[0][i].Im=0.0;
    }
    if(CE_Mka>0) 
    {
        // calculate the exp(-ikx*rx) when kx=1 (for all atoms)
        sc1=CE_scx[1];
        for(i=0;i<_NP;i++)
            sc1[i].ExpI(2.0*M_PI*_SR[i].x);
        // calculate all of the term for kx>=2 (for all atoms)
        for(ki=2;ki<=CE_Mka;ki++) for(i=0;i<_NP;i++)
            EComplex::Mul(CE_scx[ki-1][i], sc1[i], CE_scx[ki][i]);
    }
    if(CE_Mkb>0) 
    {
        // calculate the exp(-iky*ry) when ky=1 (for all atoms)
        sc1=CE_scy[1];
        for(i=0;i<_NP;i++)
	    sc1[i].ExpI(2.0*M_PI*_SR[i].y);
        // calculate all of the term for ky>=2 (for all atoms)
        for(ki=2;ki<=CE_Mkb;ki++) for(i=0;i<_NP;i++)
            EComplex::Mul(CE_scy[ki-1][i], sc1[i], CE_scy[ki][i]);
    }
    if(CE_Mkc>0) 
    {
        // calculate the exp(-ikz*rz) when kz=1 (for all atoms)
        sc1=CE_scz[1];
        for(i=0;i<_NP;i++)
	    sc1[i].ExpI(2.0*M_PI*_SR[i].z);
        // calculate all of the terms when kz=>2 (for all atoms)
        for(ki=2;ki<=CE_Mkc;ki++) for(i=0;i<_NP;i++)
            EComplex::Mul(CE_scz[ki-1][i], sc1[i], CE_scz[ki][i]);
    }
}

//-------------------------------------------------------------------
// Enumerate K-Vectors
//-------------------------------------------------------------------
void MDFrame::CE_enumKV()
{
    //------------------------------------------------------------------- 
    //  In this function, following variables are prepared
    //  Mka, Mkb, Mkc; MKV, nKV, KV
    //-------------------------------------------------------------------
    const double INHOM=1.1; //The expected error in k-point estimation
    
    EComplex *p;
    int i;
    double kx, ky, kz, k2;
    int nka, nkb, nkc, ika, ikb, ikc;
        
    // Kc = (2*pi/Li)*ni and ni = Kc*Li/(2*pi)
    // So we can multiply Kc and 1/Gi to get the integer ni
    nka=(int)(CE_Kc*Ewald_H_L1/(2.0*M_PI));
    nkb=(int)(CE_Kc*Ewald_H_L2/(2.0*M_PI));
    nkc=(int)(CE_Kc*Ewald_H_L3/(2.0*M_PI));
    
    // Get the max-number of K-vectors (half sphere)
    ika=(int)(INHOM/3.0*CE_Kc*CE_Kc*CE_Kc*Ewald_cell_V/(2.0*M_PI*2.0*M_PI));

    // Allocate K-vectors
    if(CE_MKV < ika) {CE_MKV=ika; Realloc(CE_KV,K_rec,CE_MKV);}
    
    // find K-vectors
    CE_nKV=0;
    CE_Mka=CE_Mkb=CE_Mkc=0;
    for(ika=0;ika<nka+1;ika++)
        for(ikb= (ika==0? 0:-nkb);ikb<=nkb;ikb++)
            for(ikc=((ika==0 && ikb==0)?1:-nkc);ikc<=nkc;ikc++)
            {
                //Notice that this (kx,ky,kz) is row vector, mul in left (??)
                kx=(2.0*M_PI)*(ika/Ewald_H_L1);
                ky=(2.0*M_PI)*(ikb/Ewald_H_L2);
                kz=(2.0*M_PI)*(ikc/Ewald_H_L3);
                k2=kx*kx+ky*ky+kz*kz;
                if(k2 > CE_Kc*CE_Kc) continue;
                CE_KV[CE_nKV].kx=kx;
                CE_KV[CE_nKV].ky=ky;
                CE_KV[CE_nKV].kz=kz;
                CE_KV[CE_nKV].k2=k2;
                CE_KV[CE_nKV].a=ika;
                CE_KV[CE_nKV].b=ikb;
                CE_KV[CE_nKV].c=ikc;
                CE_nKV++;
		/*
		if(CE_nKV==10000)
		  {
		    INFO("_H_L1 = "<<Ewald_H_L1);
		    INFO("_H_L2 = "<<Ewald_H_L2);
		    INFO("_H_L3 = "<<Ewald_H_L3);
		    INFO("kx = "<<kx);
		    INFO("ky = "<<ky);
		    INFO("kz = "<<kz);
		    INFO("k2 = "<<k2);
		  }
		*/
	    	
                if(CE_nKV >= CE_MKV) INFO("Error : enumKV() -> Out of KV space");
                if(CE_Mka < ika) CE_Mka=ika;
                if(CE_Mkb < abs(ikb)) CE_Mkb=abs(ikb); // Mkb can be negative
                if(CE_Mkc < abs(ikc)) CE_Mkc=abs(ikc); // Mkc can be negative
            }
    if(CE_scx==NULL) p=0;
    else p=CE_scx[0]; // assign p to scx[0]
    Realloc(p,EComplex,(_NP*(CE_Mka+CE_Mkb+CE_Mkc+3)));
    INFO("CE - CE_enumKV() : Number of Sincos Table Entries = "<<(_NP*(CE_Mka+CE_Mkb+CE_Mkc+3)));
    if(p==NULL) INFO("Error : CE_enumKV() -> Out of memory allocating sin/cos tables");
    Realloc(CE_scx,EComplex*,(CE_Mka+CE_Mkb+CE_Mkc+3)); // need to study more this part
    if(CE_scx==NULL) INFO("Error : CE_enumKV() -> Out of memory allocating sin/cos tables");
    CE_scy=CE_scx+(CE_Mka+1);
    CE_scz=CE_scy+(CE_Mkb+1);
    
    CE_scx[0]=p; 
    for(i=1;i<=CE_Mka;i++) CE_scx[i]=CE_scx[i-1]+_NP;
    CE_scy[0]=CE_scx[0]+(CE_Mka+1)*_NP;
    for(i=1;i<=CE_Mkb;i++) CE_scy[i]=CE_scy[i-1]+_NP;
    CE_scz[0]=CE_scy[0]+(CE_Mkb+1)*_NP;
    for(i=1;i<=CE_Mkc;i++) CE_scz[i]=CE_scz[i-1]+_NP;
}


//###################################################################
//  PME Functions
//###################################################################

//-------------------------------------------------------------------
// Warning
// Note that if you want to switch conj_fixbox (from 1 to 0) or (from 
// 0 to 1) in the middle of simulation, PME should be reinitialized.
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// Particle Mesh Ewald
//-------------------------------------------------------------------
void MDFrame::PME()
{
    int i;
    CE_Real();
    PME_Rec();

    _EPOT_Ewald = _EPOT_Ewald_Real + _EPOT_Ewald_Rec + _EPOT_Ewald_Self;
    _VIRIAL_Ewald = _VIRIAL_Ewald_Real + _VIRIAL_Ewald_Rec;

    for(i=0;i<_NP;i++)
    {
        _F_Ewald[i] = _F_Ewald_Real[i] + _F_Ewald_Rec[i];
        _EPOT_IND_Ewald[i] = _EPOT_IND_Ewald_Real[i] + _EPOT_IND_Ewald_Rec[i];
    }
}

//-------------------------------------------------------------------
//  Initialize Particle Mesh Ewald
//-------------------------------------------------------------------
void MDFrame::PME_init()
{    
    int m1, m2, m3, tmp, tmp_index; // k-space indices
    ComplexType PME_BC_b1, PME_BC_b2, PME_BC_b3; // for the B matrix
    double tmp1, tmp2; // temporary dummy variable

    // k-space mesh setting
    PME_K23=PME_K2*PME_K3;
    PME_K123=PME_K1*PME_K23; // equal to the number of data points

    PME_K1d=((double)PME_K1);
    PME_K2d=((double)PME_K2);
    PME_K3d=((double)PME_K3);
    PME_K123d=((double)PME_K123);
    
    PME_K3S=(PME_K3/2)+1;
    PME_K23S=PME_K2*PME_K3S;
    PME_K123S=PME_K1*PME_K2*PME_K3S;
        
    // allocate array BC,Q,IQ,CONV, etc.
    Realloc(PME_B,double,PME_K123);
    Realloc(PME_C,double,PME_K123);
    Realloc(PME_BC,double,PME_K123);
    Realloc(PME_Q,double,(2*PME_K123S));
    Realloc(PME_IQ,ComplexType,PME_K123S);
    Realloc(PME_CONV,double,(2*PME_K123S));
    
    // scaled coefficients _UR
    Realloc(_UR,Vector3,_NP);
    
    // Mn(k+1) and factorial arrays
    Realloc(PME_bsp_Mk,double,(PME_bsp_n-1)); // caution : type double total (n-1) entries
    Realloc(PME_bsp_fac,double,(PME_bsp_n+1)); // caution : type integer total (n+1) entries
    
    // start and end indices of k-space mesh point
    Realloc(PME_m1s,int,_NP);
    Realloc(PME_m2s,int,_NP);
    Realloc(PME_m3s,int,_NP);
    
    // effective k-space indices
    tmp=PME_bsp_n*_NP;
    Realloc(PME_m1,int,PME_bsp_n);
    Realloc(PME_m2,int,PME_bsp_n);
    Realloc(PME_m3,int,PME_bsp_n);
    Realloc(PME_m1mod,int,tmp);
    Realloc(PME_m2mod,int,tmp);
    Realloc(PME_m3mod,int,tmp);
    
    // (u-k),Mn(u-k),Mn-1(u-k),Mn-1(u-k-1) arrays
    Realloc(PME_x1,double,tmp);
    Realloc(PME_x2,double,tmp);
    Realloc(PME_x3,double,tmp);
    Realloc(PME_MU1,double,tmp);
    Realloc(PME_MU2,double,tmp);
    Realloc(PME_MU3,double,tmp);
    Realloc(PME_d1MU1,double,PME_bsp_n);
    Realloc(PME_d1MU2,double,PME_bsp_n);
    Realloc(PME_d1MU3,double,PME_bsp_n);
    Realloc(PME_d2MU1,double,PME_bsp_n);
    Realloc(PME_d2MU2,double,PME_bsp_n);
    Realloc(PME_d2MU3,double,PME_bsp_n);

    // prepare data which will be used in the calculation
    PME_cal_bsp_fac();
    PME_cal_bsp_Mk();
    
    // calculate matrix B and C
    for (m1=0;m1<PME_K1;m1++) // from 0 to (K1-1)
    {
        for (m2=0;m2<PME_K2;m2++) // from 0 to (K2-1)
        {
            for (m3=0;m3<PME_K3;m3++) // from 0 to (K3-1)
            {
                tmp_index=m3+m2*PME_K3+m1*PME_K23;

                // calculate each B and C component
                PME_cal_bsp_Bm(m1,PME_K1); PME_BC_b1[0]=PME_R_bsp_Bm[0];  PME_BC_b1[1]=PME_R_bsp_Bm[1];
                PME_cal_bsp_Bm(m2,PME_K2); PME_BC_b2[0]=PME_R_bsp_Bm[0];  PME_BC_b2[1]=PME_R_bsp_Bm[1];
                PME_cal_bsp_Bm(m3,PME_K3); PME_BC_b3[0]=PME_R_bsp_Bm[0];  PME_BC_b3[1]=PME_R_bsp_Bm[1];
                tmp1=PME_cal_bsp_Cm(m1,m2,m3);
                tmp2=(SQR(PME_BC_b1[0])+SQR(PME_BC_b1[1]))*(SQR(PME_BC_b2[0])+SQR(PME_BC_b2[1]))*(SQR(PME_BC_b3[0])+SQR(PME_BC_b3[1]));
                
                // compute BC matrix entries
		PME_B[tmp_index]=tmp2;
		PME_C[tmp_index]=tmp1;
                PME_BC[tmp_index]=tmp1*tmp2;
            }
        }
    }

    // prepare fft plans
#if _USEFFTW
    PME_fft_plan1=fftw_plan_dft_r2c_3d(PME_K1,PME_K2,PME_K3,PME_Q,PME_IQ,FFTW_ESTIMATE);
    PME_fft_plan2=fftw_plan_dft_c2r_3d(PME_K1,PME_K2,PME_K3,PME_IQ,PME_CONV,FFTW_ESTIMATE);
#else
    PME_fft_plan1 = 0;
    PME_fft_plan2 = 0;
#endif
    

    // Print Out PME Information
    INFO("------------------------------------------------------------------------------");
    INFO("------------------------------------------------------------------------------");
    INFO("            Particle Mesh Ewald Parameters [Long Range Interaction] ");
    INFO("------------------------------------------------------------------------------");
    INFO("    Real Space Cutoff Radius (Rc)                 = "<<Ewald_Rcut);
    INFO("    Ewald Gaussian Factor (Alpha)                 = "<<Ewald_Alpha);
    INFO("    Mesh Points in X-direction (PME_K1)           = "<<PME_K1);
    INFO("    Mesh Points in Y-direction (PME_K2)           = "<<PME_K2);
    INFO("    Mesh Points in Z-direction (PME_K3)           = "<<PME_K3);
    INFO("    Order of B-spline Interpolation (PME_bsp_n)   = "<<PME_bsp_n);
    INFO("------------------------------------------------------------------------------");
    INFO("------------------------------------------------------------------------------");
}

//-------------------------------------------------------------------
//  Calculate Factorials
//-------------------------------------------------------------------
void MDFrame::PME_cal_bsp_fac()
{
    int i, j, mult;
    
    // calculate factorials (0! to PME_bsp_n!)
    PME_bsp_fac[0]=1.0; // 0!=1
    for (i=1;i<=PME_bsp_n;i++) // from 1 to n
    {
        mult=1;
        for (j=1;j<=i;j++)
        {
            mult=mult*j;
        }
        PME_bsp_fac[i]=((double)mult);
        //INFO("[PME] Calculating Factorial => "<<i<<"! = "<<mult);
    }
}

//-------------------------------------------------------------------
//  Calculate Mn(n,k+1)
//-------------------------------------------------------------------
void MDFrame::PME_cal_bsp_Mk()
{
    int i;
    
    // calculate Mn(k+1)
    for (i=0;i<(PME_bsp_n-1);i++) // from 0 to n-2
    {
        PME_bsp_Mk[i]=PME_cal_bsp_Mu(PME_bsp_n,((double)(i+1))); // Mn(n,i+1) => from 1 to n-1
        //INFO("[PME] bsp_Mk["<<i<<"] = "<<PME_bsp_Mk[i]);
    }
}

//-------------------------------------------------------------------
//  Calculate Mn(n,u)
//-------------------------------------------------------------------
double MDFrame::PME_cal_bsp_Mu(int n, double u)
{
    int i,k;
    double bsp_sign, bsp_Mu, bsp_Muk, tmp;

    // calculate Mn(u-k)
    bsp_sign=1.0; // start from +1
    bsp_Mu=0.0; // initialize sum quantity
    for (k=0;k<=n;k++) // from 0 to n
    {
        tmp=u-((double)k); // convert k to double and calculate (u-k)
        if(tmp>0.0) // only when (u-k)>0
        {
            // calculate (u-k)^(n-1)
            bsp_Muk=1.0;
            for (i=1;i<n;i++) // from 1 to n-1
            {
                bsp_Muk*=tmp;
            }
            // calculate n!/[k!(n-k)!]*(u-k)^(n-1) and sum it up
            bsp_Mu+=(bsp_sign*PME_bsp_fac[n]/(PME_bsp_fac[k]*PME_bsp_fac[n-k])*bsp_Muk);
            bsp_sign*=-1.0;
        }
    }
    bsp_Mu/=PME_bsp_fac[n-1]; // divide by (n-1)!
    return bsp_Mu;
}

//-------------------------------------------------------------------
//  Calculate bi(mi)
//-------------------------------------------------------------------
void MDFrame::PME_cal_bsp_Bm(int mi, int Ki)
{
    int k;
    double tmp;
    ComplexType bsp_bm, bsp_bms, bsp_bmn;
    
    // initialize the summation dummy complex number
    bsp_bms[0]=0.0;
    bsp_bms[1]=0.0;
    
    // calculate the B matrix coefficient Bi(mi)
    // compute denominator c+di
    for (k=0;k<(PME_bsp_n-1);k++) // from 0 to (n-2)
    {
        tmp=2.0*M_PI*(((double)(mi*k))/((double)Ki));
        bsp_bms[0]+=PME_bsp_Mk[k]*cos(tmp);
        bsp_bms[1]+=PME_bsp_Mk[k]*sin(tmp);
    }
    // compute numerator a+bi
    tmp=2.0*M_PI*(((double)((PME_bsp_n-1)*mi))/((double)Ki));
    bsp_bmn[0]=cos(tmp);
    bsp_bmn[1]=sin(tmp);
    
    // compute c^2+d^2
    tmp=SQR(bsp_bms[0])+SQR(bsp_bms[1]);
    
    // compute (a+bi)/(c+di) = [(ac+bd)+(bc-ad)i]/(c^2+d^2)
    bsp_bm[0]=(bsp_bmn[0]*bsp_bms[0]+bsp_bmn[1]*bsp_bms[1])/tmp; //(ac+bd)
    bsp_bm[1]=(bsp_bmn[1]*bsp_bms[0]-bsp_bmn[0]*bsp_bms[1])/tmp; //(bc-ad)
    PME_R_bsp_Bm[0]=bsp_bm[0];
    PME_R_bsp_Bm[1]=bsp_bm[1];
}

//-------------------------------------------------------------------
//  Calculate C(m1,m2,m3)
//-------------------------------------------------------------------
double MDFrame::PME_cal_bsp_Cm(int m1, int m2, int m3)
{
    double msq, m1s, m2s, m3s, tmp1, tmp2, bsp_Cm;
    
    if((m1==0)&&(m2==0)&&(m3==0))
    {
        bsp_Cm=0.0;
    }
    else
    {
        if(m1>(PME_K1/2)){m1=m1-PME_K1;}
        if(m2>(PME_K2/2)){m2=m2-PME_K2;}
        if(m3>(PME_K3/2)){m3=m3-PME_K3;}
        m1s=((double)m1)/Ewald_H_L1; // m1/L1 => double
        m2s=((double)m2)/Ewald_H_L2; // m2/L2 => double
        m3s=((double)m3)/Ewald_H_L3; // m3/L3 => double
        msq=SQR(m1s)+SQR(m2s)+SQR(m3s);
        tmp1=-SQR(M_PI)*msq/(SQR(Ewald_Alpha));
        tmp2=1.0/M_PI/Ewald_cell_V/msq;
        bsp_Cm=tmp2*exp(tmp1);
    }
    return bsp_Cm;
}

//-------------------------------------------------------------------
//  Particle Mesh Ewald Summation
//-------------------------------------------------------------------
void MDFrame::PME_Rec()
{
    int i, j, i1, i2, i3, m1, m2, m3, pm1, pm2, pm3;
    int tmp_ind_i, tmp_ind_1, tmp_ind_2, tmp_ind_3, tmp_index;
    double tmp_double, tmp_charge, m1s, m2s, m3s, msq, pfac, UE;

     // start time of reciprocal part
    PME_sREC=clock();

    // recalculate system size and volume (in case the box is changing)
    if(conj_fixbox==0)
      {
        Ewald_H_L1=sqrt(SQR(_H[0][0])+SQR(_H[1][0])+SQR(_H[2][0]));
        Ewald_H_L2=sqrt(SQR(_H[0][1])+SQR(_H[1][1])+SQR(_H[2][1]));
        Ewald_H_L3=sqrt(SQR(_H[0][2])+SQR(_H[1][2])+SQR(_H[2][2]));
        Ewald_cell_V=_H.det();
      }

    // initialize Q and recalculate C and BC if neccessary
    for (i1=0;i1<PME_K1;i1++) // from 0 to (K1-1)
      {
	for (i2=0;i2<PME_K2;i2++) // from 0 to (K2-1)
	  {
	    for (i3=0;i3<PME_K3;i3++) // from 0 to (K3-1)
	      {
		tmp_index=i3+i2*PME_K3+i1*PME_K23;
	
		if(conj_fixbox==0)
		  {
		    // calculate each B and C component
		    tmp_double=PME_cal_bsp_Cm(i1,i2,i3);
		
		    // compute BC matrix entries
		    PME_C[tmp_index]=tmp_double;
		    PME_BC[tmp_index]=PME_B[tmp_index]*tmp_double;
		  }

		// initialize Q matrix
		PME_Q[tmp_index]=0.0;
	      }
	  }
      }

    // initialize force, potential and virials
    _EPOT_Ewald_Rec = 0.0;
    _VIRIAL_Ewald_Rec.clear();
    for(i=0;i<_NP;i++)
      {
        _F_Ewald_Rec[i].clear();
        _EPOT_IND_Ewald_Rec[i]=0.0;
      }

    // get the scaled and manipulated coordinates UR
    for (i=0;i<_NP;i++)
      {
        _UR[i].x=(_SR[i].x+0.5)*PME_K1d;
        _UR[i].y=(_SR[i].y+0.5)*PME_K2d;
        _UR[i].z=(_SR[i].z+0.5)*PME_K3d;
      }
    
    // calculate Q matrix
    for (i=0;i<_NP;i++) // for all atoms
      {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[i]==-1) continue;
        
        tmp_ind_i=PME_bsp_n*i;
        tmp_charge=P_SQRT_CLMB*_ATOMCHARGE[species[i]];
        
        // since function bsp_Mu has bounded support, we can find the initial and
        // final k-space indices of each atom where Q matrix entries are non-zero
        // [ start index : floor(u)-n+1    end index : floor(u) ]
        
        PME_m1s[i]=((int)floor(_UR[i].x))-PME_bsp_n+1;
        PME_m2s[i]=((int)floor(_UR[i].y))-PME_bsp_n+1;
        PME_m3s[i]=((int)floor(_UR[i].z))-PME_bsp_n+1;
	
        for (j=0;j<PME_bsp_n;j++) // from 0 to (n-1)
	  {
            // calculate temporary index
            tmp_index =j+tmp_ind_i;
            
            // goes through effective range of k-space for atoms i
            PME_m1[j]=PME_m1s[i]+j; // in x-direction
            PME_m2[j]=PME_m2s[i]+j; // in y-direction
            PME_m3[j]=PME_m3s[i]+j; // in z-direction
            
            // apply periodic boundary condition
            PME_m1mod[tmp_index]=PME_int_mod(PME_m1[j],PME_K1);
            PME_m2mod[tmp_index]=PME_int_mod(PME_m2[j],PME_K2);
            PME_m3mod[tmp_index]=PME_int_mod(PME_m3[j],PME_K3);
            
            // calculate (u-k)
            PME_x1[tmp_index]=_UR[i].x-((double)PME_m1[j]);
            PME_x2[tmp_index]=_UR[i].y-((double)PME_m2[j]);
            PME_x3[tmp_index]=_UR[i].z-((double)PME_m3[j]);
            
            // calculate Mn(u) coefficients
            PME_MU1[tmp_index]=PME_cal_bsp_Mu(PME_bsp_n,PME_x1[tmp_index]);
            PME_MU2[tmp_index]=PME_cal_bsp_Mu(PME_bsp_n,PME_x2[tmp_index]);
            PME_MU3[tmp_index]=PME_cal_bsp_Mu(PME_bsp_n,PME_x3[tmp_index]);
	  }
        
        // scan the atomwise small cubic which has nonzero entries
        for (i1=0;i1<PME_bsp_n;i1++) // from 0 to (n-1)
	  {
            for (i2=0;i2<PME_bsp_n;i2++) // from 0 to (n-1)
	      {
                for (i3=0;i3<PME_bsp_n;i3++) // from 0 to (n-1)
		  {
                    tmp_ind_1=i1+tmp_ind_i;
                    tmp_ind_2=i2+tmp_ind_i;
                    tmp_ind_3=i3+tmp_ind_i;
                    tmp_index=PME_m3mod[tmp_ind_3]+PME_m2mod[tmp_ind_2]*PME_K3+PME_m1mod[tmp_ind_1]*PME_K23;
                    PME_Q[tmp_index]+=tmp_charge*PME_MU1[tmp_ind_1]*PME_MU2[tmp_ind_2]*PME_MU3[tmp_ind_3];
		  }
	      }
	  }
      }
    
    // out-of-place backward fft of matrix PME_Q (backward fft in the paper => forward fft for the fftw)
#ifdef _USEFFTW
    fftw_execute(PME_fft_plan1); // real to complex
#else
    INFO("FFTW not installed! Cannot use PME!");
#endif

    // calculate energy and stress and
    // matrix multiplication of PME_BC and PME_IQ -> (IQ*BC)
    for (m1=0;m1<PME_K1;m1++) // from 0 to (PME_K1-1)
      {
        for (m2=0;m2<PME_K2;m2++) // from 0 to (PME_K2-1)
	  {
            for (m3=0;m3<PME_K3S;m3++) // from 0 to (PME_K3S-1)
	      {
                tmp_ind_1=m3+m2*PME_K3S+m1*PME_K23S;
                tmp_index=m3+m2*PME_K3+m1*PME_K23;
                tmp_double=PME_BC[tmp_index]/PME_K123d;
		UE=0.0;
		
		// calculate the reciprocal energy
		if((m3==(PME_K3S-1))||(m3==0))
		  {UE=0.5*PME_BC[tmp_index]*(SQR(PME_IQ[tmp_ind_1][0])+SQR(PME_IQ[tmp_ind_1][1]));}
		else
		  {UE=PME_BC[tmp_index]*(SQR(PME_IQ[tmp_ind_1][0])+SQR(PME_IQ[tmp_ind_1][1]));}
		_EPOT_Ewald_Rec+=UE;
		
		// calculate BC*IQ
                PME_IQ[tmp_ind_1][0]*=tmp_double;
                PME_IQ[tmp_ind_1][1]*=tmp_double;
		
		// apply periodic condition
		if(m1>(PME_K1/2)){pm1=m1-PME_K1;}
		else{pm1=m1;}
		if(m2>(PME_K2/2)){pm2=m2-PME_K2;}
		else{pm2=m2;}
		if(m3>(PME_K3/2)){pm3=m3-PME_K3;}
		else{pm3=m3;}

		//calculate the pre-factor
		m1s=((double)pm1)/Ewald_H_L1; // m1/L1 => double
		m2s=((double)pm2)/Ewald_H_L2; // m2/L2 => double
		m3s=((double)pm3)/Ewald_H_L3; // m3/L3 => double
		msq=SQR(m1s)+SQR(m2s)+SQR(m3s);
		pfac=(1.0/msq+(SQR(M_PI))/(SQR(Ewald_Alpha)));
		
		// Virial Stress Calculation
		if(!((m1==0)&&(m2==0)&&(m3==0)))
		  {
		    _VIRIAL_Ewald_Rec[0][0]+=(UE*(1.0-2.0*((m1s*m1s))*pfac)); //xx
		    _VIRIAL_Ewald_Rec[1][1]+=(UE*(1.0-2.0*((m2s*m2s))*pfac)); //yy
		    _VIRIAL_Ewald_Rec[2][2]+=(UE*(1.0-2.0*((m3s*m3s))*pfac)); //zz
		    _VIRIAL_Ewald_Rec[0][1]+=(UE*(-2.0*((m1s*m2s))*pfac));  //xy
		    _VIRIAL_Ewald_Rec[1][0]+=(UE*(-2.0*((m1s*m2s))*pfac));  //yx
		    _VIRIAL_Ewald_Rec[0][2]+=(UE*(-2.0*((m1s*m3s))*pfac));  //xz
		    _VIRIAL_Ewald_Rec[2][0]+=(UE*(-2.0*((m1s*m3s))*pfac));  //zx
		    _VIRIAL_Ewald_Rec[1][2]+=(UE*(-2.0*((m2s*m3s))*pfac));  //yz
		    _VIRIAL_Ewald_Rec[2][1]+=(UE*(-2.0*((m2s*m3s))*pfac));  //zy
		  }
	      }
	  }
      }
    
    // in-place forward fft of matrix PME_CONV (forward fft in the paper => backward fft for the fftw)
#ifdef _USEFFTW
    fftw_execute(PME_fft_plan2); // complex to real
#else
    INFO("FFTW not installed! Cannot use PME!");
#endif

    /*
    // matrix multiplication of PME_Q and PME_CONV and calculate the reciprocal energy simultaneously
    for (m1=0;m1<PME_K1;m1++) // from 0 to (PME_K1-1)
    {
        for (m2=0;m2<PME_K2;m2++)  // from 0 to (PME_K2-1)
        {
            for (m3=0;m3<PME_K3;m3++) // from 0 to (PME_K3-1)
            {
                tmp_index=m3+m2*PME_K3+m1*PME_K23;
                _EPOT_Ewald_Rec+=(PME_Q[tmp_index]*PME_CONV[tmp_index]);
            }
        }
    }
    _EPOT_Ewald_Rec*=(PME_K123d/(2.0));
    //INFO("[PME] Reciprocal Energy Calculated by PME = "<<_EPOT_Ewald_Rec);
    */
    
    // force calculation - Calculate dQ/dr matrix
    for (i=0;i<_NP;i++) // for all atoms
      {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[i]==-1) continue;
        
        tmp_ind_i=PME_bsp_n*i;
        for (j=0;j<PME_bsp_n;j++) // from 0 to (n-1)
	  {
            // calculate temporary index
            tmp_index =j+tmp_ind_i;
            
            PME_d1MU1[j]=PME_cal_bsp_Mu((PME_bsp_n-1),PME_x1[tmp_index]); // Mn-1(u1-k1)
            PME_d1MU2[j]=PME_cal_bsp_Mu((PME_bsp_n-1),PME_x2[tmp_index]); // Mn-1(u2-k2)
            PME_d1MU3[j]=PME_cal_bsp_Mu((PME_bsp_n-1),PME_x3[tmp_index]); // Mn-1(u3-k3)
            
            PME_d2MU1[j]=PME_cal_bsp_Mu((PME_bsp_n-1),(PME_x1[tmp_index]-1.0)); // Mn-1(u1-k1-1)
            PME_d2MU2[j]=PME_cal_bsp_Mu((PME_bsp_n-1),(PME_x2[tmp_index]-1.0)); // Mn-1(u2-k2-1)
            PME_d2MU3[j]=PME_cal_bsp_Mu((PME_bsp_n-1),(PME_x3[tmp_index]-1.0)); // Mn-1(u3-k3-1)
	  }
        
        // scan the atomwise small cubic which has nonzero entries
        for (i1=0;i1<PME_bsp_n;i1++) // from 0 to (n-1)
	  {
            for (i2=0;i2<PME_bsp_n;i2++) // from 0 to (n-1)
	      {
                for (i3=0;i3<PME_bsp_n;i3++) // from 0 to (n-1)
		  {
                    tmp_ind_1=i1+tmp_ind_i;
                    tmp_ind_2=i2+tmp_ind_i;
                    tmp_ind_3=i3+tmp_ind_i;
                    tmp_index=PME_m3mod[tmp_ind_3]+PME_m2mod[tmp_ind_2]*PME_K3+PME_m1mod[tmp_ind_1]*PME_K23;
                    _F_Ewald_Rec[i].x+=((PME_d2MU1[i1]-PME_d1MU1[i1])*PME_MU2[tmp_ind_2]*PME_MU3[tmp_ind_3]*PME_CONV[tmp_index]);
                    _F_Ewald_Rec[i].y+=((PME_d2MU2[i2]-PME_d1MU2[i2])*PME_MU1[tmp_ind_1]*PME_MU3[tmp_ind_3]*PME_CONV[tmp_index]);
                    _F_Ewald_Rec[i].z+=((PME_d2MU3[i3]-PME_d1MU3[i3])*PME_MU1[tmp_ind_1]*PME_MU2[tmp_ind_2]*PME_CONV[tmp_index]);
		  }
	      }
	  }
        _F_Ewald_Rec[i].x*=(P_SQRT_CLMB*_ATOMCHARGE[species[i]]*PME_K123d*PME_K1d/Ewald_H_L1);
        _F_Ewald_Rec[i].y*=(P_SQRT_CLMB*_ATOMCHARGE[species[i]]*PME_K123d*PME_K2d/Ewald_H_L2);
        _F_Ewald_Rec[i].z*=(P_SQRT_CLMB*_ATOMCHARGE[species[i]]*PME_K123d*PME_K3d/Ewald_H_L3);
	//_F_Ewald_Rec[i].x*=PME_C1[species[i]];
	//_F_Ewald_Rec[i].y*=PME_C2[species[i]];
	//_F_Ewald_Rec[i].z*=PME_C3[species[i]];
      }
    
    // end time of reciprocal part
    PME_eREC=clock();
}

//-------------------------------------------------------------------
//  Clear Particle Mesh Ewald
//-------------------------------------------------------------------
void MDFrame::PME_clear()
{
    // free variables used in the Particle Mesh Ewald summation
    free(_UR); free(PME_bsp_Mk); free(PME_bsp_fac);
    
    free(PME_m1s); free(PME_m2s); free(PME_m3s);
    free(PME_m1); free(PME_m2); free(PME_m3);
    free(PME_m1mod); free(PME_m2mod); free(PME_m3mod);
    free(PME_x1); free(PME_x2); free(PME_x3);
    
    free(PME_MU1); free(PME_MU2); free(PME_MU3);
    free(PME_d1MU1); free(PME_d1MU2); free(PME_d1MU3);
    free(PME_d2MU1); free(PME_d2MU2); free(PME_d2MU3);
    
    free(PME_Q); free(PME_IQ);
    free(PME_BC); free(PME_CONV);

#ifdef _USEFFTW    
    fftw_destroy_plan(PME_fft_plan1);
    fftw_destroy_plan(PME_fft_plan2);
#endif
}

