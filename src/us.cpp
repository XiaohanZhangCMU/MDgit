/*
  us.cpp 
  taken from ffs.cpp by Seunghwa Ryu, shryu@stanford.edu
  Last Modified : Sun Aug  5 17:16:38 2012

  FUNCTION  :  Umbrella Sampling for free energy barrier calculations

  To Do:
      1.  calcrystalorder() needs to be improved

  Note:
      This file is included by md.cpp.  This file cannot be compiled by itself.
*/

/********************************************************************/
/* Umbrella Sampling (US) function (Seunghwa Ryu)              */
/********************************************************************/

void MDFrame::initUMB()
{
    allocUMB();
    if ((react_coord_type == 0) || (react_coord_type == 2))
    { 
       setconfig4();
       copyStoS4();
    }

    refreshneighborlist(); 
    cal_react_coord();

    kinetic_flag = 0;
}

void MDFrame::allocUMB()
{
    int i;
    Realloc(Narray,int,2*delta_n+1);
    Realloc(Occurence,int,2*delta_n+1);
    Realloc(Probability,double,2*delta_n+1);
    for (i=0;i<2*delta_n+1;i++)
    {
	Narray[i]=n_center-delta_n+i;
	Occurence[i]=0;
	Probability[i]=0;
    }
}

/* write histogram data to file */
void MDFrame::printhist()
{
    FILE *fp, *fp1, *fp2, *fp3;
    int i;
    fp  = fopen("Narray.txt","w");  
    fp1 = fopen("Freq.txt","w");
    fp2 = fopen("Prob.txt","w");
    fp3 = fopen("Nsample.txt","w");
    for (i=0;i<2*delta_n+1;i++)
    {
	fprintf(fp,"%d\n",Narray[i]);
	if (UMB_K_in_Kelvin >= 0 ) {
	fprintf(fp1,"%d\n",Occurence[i]);
	fprintf(fp2,"%f\n",Probability[i]);
	} else {
	fprintf(fp1,"%f\n",Probability[i]);
	fprintf(fp2,"%f\n",Probability[i]);
	}
    }
    fprintf(fp3,"%d\n%d\n",MC_UMB_num_of_trials,MC_UMB_accept_tot);
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
}

/* read histogram data from file */
void MDFrame::readhist()
{
    FILE *fp, *fp1, *fp2, *fp3;
    int i,tempint0,tempint1;
    double temp;
    fp  = fopen("Narray.txt","r");
    fp1 = fopen("Freq.txt","r");
    fp2 = fopen("Prob.txt","r");
    fp3 = fopen("Nsample.txt","w");
    for (i=0;i<2*delta_n+1;i++)
    {
	fscanf(fp,"%d",&tempint0);
        fscanf(fp1,"%d",&tempint1);
        fscanf(fp2,"%lf",&temp);
        Narray[i]=tempint0;
        Occurence[i]=tempint1;
        Probability[i]=temp;
    }
    fscanf(fp3,"%d",&MC_UMB_num_of_trials);
    fscanf(fp3,"%d",&MC_UMB_accept_tot);
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
}

void MDFrame::cal_react_coord() {

   if (react_coord_type == 0)
   { /* number of atoms inside dislocation loop */
     if ( L_for_QLM != 0 ) 
     {
        INFO("L_for_QLM = "<<L_for_QLM);
        FATAL("We have agreed to set L_for_QLM = 0 for caldislocation order");
     }
     caldislocationorder();
     react_coord = N_lgst_cluster;
   }
   else if (react_coord_type == 1)
   { /* normalized distance from config1 to state config2 times (chainlength-1) */
      calconstrainS_inst();
      react_coord = constrainS_inst * (react_coord_ngrid-1);
   }
   else if (react_coord_type == 2)
   { /* number of atoms inside largest crystal cluster */
      if ( L_for_QLM == 0 )
      {
        FATAL("Attempt to compute crystal order parameter while L_for_QLM = "<<L_for_QLM);
      }
      calcrystalorder();
     
      react_coord = N_lgst_cluster;
   }
   else if (react_coord_type == 3)
   { /* compute reaction coordinate based on relaxed MEP */
      cal_react_coord_from_MEP();
   }
   else
   {
      ERROR("cal_react_coord() unknown react_coord_type ("<<react_coord_type<<")");
   }

   //react_coord_old = react_coord;
}

void MDFrame::cal_react_coord_from_MEP() {
   int n, i, j, size;
   double *TangMag2, r2, fr, fr_old;
   void *c; Vector3 **dTang, dr;

   if(_Rc==NULL)
   {
      FATAL("_Rc is NULL");
   }
   /* compute tangent vectors along MEP */ 
   n=constrainatoms[0]; 
   size=sizeof(Vector3 *)*(_CHAINLENGTH+1)+sizeof(Vector3)*(_CHAINLENGTH+1)*n + 10;
   //INFO_Printf("n = %d size = %d\n",n,size);
   c=malloc(size); memset(c,0,size);
   dTang = (Vector3 **)c; dTang[0]=(Vector3 *)(dTang+(_CHAINLENGTH+1));
   for(i=1;i<_CHAINLENGTH+1;i++) dTang[i]=dTang[i-1]+n;
   TangMag2=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
   for(j=1;j<_CHAINLENGTH;j++) // j; index of chains
   {
       /* A more sophisticated approach uses the energy long the chain, see nebrelax() */
       r2 = 0;
       for(i=0;i<n;i++) // n: number of constrained atoms
       {
           //INFO_Printf("Rc[%d][%d], Rc[%d][%d]\n",j+1,i,j-1,i);
           dTang[j][i] = _Rc[j+1][i]-_Rc[j-1][i];
           r2+=dTang[j][i].norm2();
       }
       TangMag2[j]=r2;
       //INFO_Printf("TangMag2[%d] = %g\n",j,TangMag2[j]);
   }

   /* compute progress along MEP */
   SHtoR(); fr=0; fr_old=0;
   for(j=1;j<_CHAINLENGTH;j++)
   {
       fr=0;            
       for(i=0;i<n;i++)
       {
           dr = _R[constrainatoms[i+1]]-_Rc[j][i];
           fr+=dot(dr,dTang[j][i]);
       }
       fr/=sqrt(TangMag2[j]);
       //INFO_Printf("fr = %g\n",fr);
       if (fabs(fr)<1e-10)
       {
           react_coord = round( 1.0*j/(_CHAINLENGTH-1) * (react_coord_ngrid-1) );
           INFO_Printf("found exact match: react_coord = %g\n",react_coord);
           return;
       }
       else if (fr<0)
       {
           break;
       }
       else
       { 
           fr_old = fr;
       }
   }
   if (fr<0)
   {
       if (j<=1)
       { 
          fr_old=0;            
          for(i=0;i<n;i++)
          {
             dr = _Rc[j-1][i]-_Rc[j][i];
             fr_old+=dot(dr,dTang[j][i]);
          }
          fr_old/=sqrt(TangMag2[j]);
          react_coord = round((1.0 - fr/fr_old)/(_CHAINLENGTH-1) * (react_coord_ngrid-1));
       }
       else
       {
          react_coord = round((fr_old*j - fr*(j-1))/(fr_old-fr)/(_CHAINLENGTH-1) * (react_coord_ngrid-1));
       }
  }
  else
  {
      fr_old=0;            
      for(i=0;i<n;i++)
      {
         dr = _Rc[_CHAINLENGTH][i]-_Rc[_CHAINLENGTH][i];
         fr_old+=dot(dr,dTang[_CHAINLENGTH-1][i]);
      }
      fr_old/=sqrt(TangMag2[j]);
      react_coord = round(1.0 - (1.0 - fr/fr_old)/(_CHAINLENGTH-1) * (react_coord_ngrid-1));
  }

  //INFO_Printf("react_coord = %g\n",react_coord);
  free(dTang); free(TangMag2);
}

/* To do: move to eam.cpp, for not being general enough */
/* Wei Cai, 8/5/2012 */
/* UMB used by MD simulation run() */
void MDFrame::UMBprocess()
{
    double umb_temp,R;
    int umb_deltan;

    if ( UMB_K_in_Kelvin==0 ) {
	if ( N_lgst_cluster < Narray[0] )
	{
	    copyS4toS();
	    Occurence[0]++;
	    Probability[0]++;
	}
	else if ( N_lgst_cluster > Narray[2*delta_n])
	{
	    copyS4toS();
	    Occurence[2*delta_n]++;
	    Probability[2*delta_n]++;
	}
	else
	{
	    umb_deltan=N_lgst_cluster-Narray[0];
	    Occurence[umb_deltan]++;
	    Probability[umb_deltan]++;
	    copyStoS4();
	    N_lgst_temp=N_lgst_cluster;
	}
    }

    if ( UMB_K_in_Kelvin > 0 && kinetic_flag == 0 )
    {
	umb_temp=(N_lgst_cluster-N_lgst_temp)*(N_lgst_cluster+N_lgst_temp-2*n_center);
	R=exp(-0.5*UMB_K_in_Kelvin*(umb_temp)/_TDES);
//        INFO("check");
//        INFO_Printf("R=%f, umb_temp=%f, n_center=%d N_lgst=%d, N_lgst_temp=%d\n",R,umb_temp,n_center,N_lgst_cluster,N_lgst_temp);
	acc_UMB = (drand48()<=R);
	if (acc_UMB) {
	    umb_deltan=(N_lgst_cluster-Narray[0]);
	    umb_temp=exp(0.5*UMB_K_in_Kelvin*pow((N_lgst_cluster-n_center),2)/_TDES);
	    if (YES_HMC!=1) copyStoS4();
	    else 
	    {
		if (acc_HMC==1) copyStoS4();
	    }		
	    if (umb_deltan>=0 && umb_deltan <2*delta_n && Kinetic < 1 && curstep > UMB_equilstep) {
	        Occurence[umb_deltan]++;
	        Probability[umb_deltan]+=umb_temp;
	    }
//            INFO_Printf("UMB accepted\n");
//            INFO_Printf("R=%f, umb_temp=%f, _TDES=%f \n",R,umb_temp,_TDES);
//            INFO_Printf("N_lgst=%d, N_lgst_temp=%d \n",N_lgst_cluster,N_lgst_temp);
	    N_lgst_temp=N_lgst_cluster;
//            perturbevelocity();
//            initvelocity();
    	}
	else
	{
	    umb_deltan=(N_lgst_temp-Narray[0]);
	    umb_temp=exp(0.5*UMB_K_in_Kelvin*pow((N_lgst_temp-n_center),2)/_TDES);
	    copyS4toS();
//            initvelocity();
//            input=30.0;
//            perturbevelocity();
	    if (umb_deltan >=0 && umb_deltan <= 2*delta_n && Kinetic < 1 && curstep > UMB_equilstep ) {
	    	Occurence[umb_deltan]++;
		Probability[umb_deltan]+=umb_temp;
	    }
//            INFO_Printf("UMB rejected\n");
   	}
    } 
/*
    if ( Kinetic == 1 && kinetic_flag == 1 ) {
	kinetic0++;
	n_array[kinetic0]=n_array[kinetic0]+(N_lgst_cluster-n_center)*(N_lgst_cluster-n_center);
	if ( kinetic0 > (int) 1.0*Kinetic_Time/savepropfreq ) kinetic_flag = 0;
    }

    if ( Kinetic == 1 && N_lgst_cluster == n_center && kinetic_flag == 0 && curstep < totalsteps - Kinetic_Time - 1)
    {
	kinetic_flag=1;
	Kinetic_Swip++;
	kinetic0=0;
    }
*/
}

int MDFrame::ADVsample()
{
    int key;
    key = 0;
	    /* Forward Flux Sampling (FFS) Seunghwa Ryu */
            if (saveFFScn == 1 && curstep > step0 && (curstep%plotfreq)==0)
            {
 	        FFScurstep=curstep;
                key=FFSprocess();
                if ( key==1 ) {
		    if ( FFShist == 1 ) printFFShist();
		}
            }
  
	    /* Umbrella Sampling */
	    if (YES_UMB==1 && (curstep%savepropfreq)==0 && YES_HMC!=1 && curstep>0 ) {
		UMB_curstep=curstep;
                if (drand48()<UMB_tryrate) UMBprocess(); 
                if ( curstep%(10*savepropfreq)==0 )
			printhist();
            }
  
            /* Hybrid Monte Carlo */
            if ( YES_HMC==1 && (curstep%Ns_HMC)==0 && curstep>0 )
            {
		if (UMB_K_in_Kelvin >=0) HMCprocess();
		if (UMB_K_in_Kelvin ==-1) UMBprocess();
	    }	
   return key;
}



void MDFrame::allocUMBorder()
{
    if ( UMB_order_param == 103 ) {
	L_for_QLM = 3;allocQLM();plot_color_axis=8;
    }
    else if ( UMB_order_param == 104 ) {
        L_for_QLM = 4;allocQLM();plot_color_axis=8;
    }
    else if ( UMB_order_param == 106 ) {
        L_for_QLM = 6;allocQLM();plot_color_axis=8;
    }
    else if ( UMB_order_param < 100 ) {
        if ( UMB_order_param/10 == 4 ) assignUMBplane();
        L_for_QLM = 0;allocQLM();plot_color_axis=9;
    }
#ifdef _GSL
    else if ( UMB_order_param/10 == 11 ) {
    // setting for 111:HCPonly,112:FCConly,113:HCP+FCC,114:LIQonly
        L_for_QLM = 6;l_for_qlm = 4;allocQLM();allocqlm();
        plot_color_axis=7;
    }
#endif
    else {
    INFO_Printf("UMB_order_param Lists : \n\n- - - - SLIP and DISLOCATION PARAMETER - - - - \n10: slip analysis for dislocation nucleation \n      (caldislocationorder(), plot_color_axis =9 ) \n      defined by Ngan maximum relative displacement order parameter\n\n20: kink pair analysis for kink pair nucleation \n      (caldislocationorder(), plot_color_axis =9 )\n      defined by Ryu-Cai slip order parameter\n\n30: slip analysis for cross slip nucleation ( plot_color_axis =9 )\n      (caldislocationorder(), plot_color_axis =9 )\n      defined by Ryu-Cai slip order parameter\n\n40: general slip analysis for kink pair nucleation \n      (caldislocationorder(), plot_color_axis =9 )\n      defined by Ryu-Cai slip order parameter\n\n        41: homogeneous dislocation,\n        42: nanorod dislocation,\n        43: cross slip\n\n- - - - CRYSTAL ORDER PARAMETER - - - - \n\n10X : crystal order parameter by \n         Tianshu Lis definition : \n         (calcrystalorder() with L=6, plot_color_axis =8 ) \n\n103: crystal order with L = 3, find DC exclusively\n104: crystal order with L = 4, find DC+FCC+HCP+BCC order together\n106: crystal order with L = 6, find DC+FCC+HCP+BCC order together\n\n11X: crystal order parameter by\n        HCP - FCC, recognition :\n        (calqlwl() function with L=6 & l=4, plot_color_axis=7)\n        TOPOL is assigned by   HCP : - 0.5  /  FCC : +0.5  /  LIQ : 0.0\n\n111: HCP only\n112: FCC only\n113: HCP+FCC only\n114: LIQ only\n");
    ERROR("UMB_order_param ("<< UMB_order_param <<") is not a option\n");return;
    }
}

/*
  Alloc memory for md.cpp::
 **UMB_nindex_ref, *UMB_nn_ref, *UMB_nindex_ref_mem 
*/
void MDFrame::assignUMBindexref()
{
  int i, j,jpt, n;
  double rrij, temp1;
  Vector3 sij, rij;
  int shft1,shft2, mx, mz;
  mx = _NP;
  mz = _NNM;
  shft1=mx*mz*sizeof(int);
  shft2=mx*sizeof(int *);
  if(shft1+shft2==0) return;
    
  Realloc(UMB_nindex_ref_mem,char,(shft1+shft2));

  UMB_nindex_ref=(int **)(UMB_nindex_ref_mem+shft1);

  memset(UMB_nindex_ref_mem,0,shft1+shft2);
  for(i=0;i<mx;i++)
    {
      UMB_nindex_ref[i]=(int *)(UMB_nindex_ref_mem+i*mz*sizeof(int));
    }

  Realloc(UMB_nn_ref,int,mx);

  //printf("I am here 0\n");
  //Now we assign values to the arrs.
  
  refreshneighborlist(); 
  
  for (i=0;i< _NP;i++)
    {
      n = 0;
      //printf("I am here 0.0, nn[i] = %d\n", nn[i]);
      for(j=0;j<nn[i];j++)
        {
	  //printf("I am here 0.1\n");
	  jpt=nindex[i][j];
	  if(jpt==i) continue;
	  sij=_SR2[jpt]-_SR2[i]; sij.subint();
	  rij=_H*sij; rrij=rij.norm();

	  /* select nearest neighbor atoms */
	  if (rrij > Rc_for_QLM) continue;
	  //printf("I am here 0.2\n");
	  /* select atoms on the adjacent plane */
	  temp1 = dot(rij,UMB_nvec)/( UMB_thickness );
	  if (temp1 < UMB_HMIN || temp1 > UMB_HMAX) continue; //glide partial
	  //printf("I am here 0.3\n");
	  UMB_nindex_ref[i][n] = jpt;
	  n++;
	}
      UMB_nn_ref[i] = n;
     // printf("i = %d, n = %d\n", i,n);
    }
  //printf("I am here 1\n");
}

void MDFrame::assignUMBslipvec()
{
    UMB_slipvec.set(input[0],input[1],input[2]);
    UMB_slipvec=UMB_slipvec/UMB_slipvec.norm();
}

void MDFrame::assignUMBplane()
{
    for(int i=0;i<_NP;i++) _R[i]=_H*_SR2[i];

    Vector3 RU0, RU1, RU2, RD0, NV, NN, dR1, dR2, dR;
    int u0ind, u1ind, u2ind, d0ind;
    double distance;

   if ( UMB_order_param >= 41 && UMB_order_param <= 44 ) {
    if ( UMB_order_param < 41 || UMB_order_param > 43 ) {
        INFO("read in atom IDs from input");
        u0ind=(int)round(input[0]);
        u1ind=(int)round(input[1]);
        u2ind=(int)round(input[2]);
        d0ind=(int)round(input[3]);
    }
    else if ( UMB_order_param == 41 ) {
        WARNING("Beware of hard coded atom IDs");
	u0ind=10473;u1ind=6728;u2ind=8547;d0ind=8548;
        UMB_order_group=102506;
    }
    else if ( UMB_order_param == 42 ) {
        WARNING("Beware of hard coded atom IDs");
        u0ind=16843;u1ind=9667;u2ind=17299;d0ind=12365;
        UMB_order_group=102506;
    }
    else if ( UMB_order_param == 43 ) {
        WARNING("Beware of hard coded atom IDs");
        u0ind=190823;u1ind=144738;u2ind=132216;d0ind=166316;
        UMB_order_group=102506;
    }
    else {
        u0ind=-1; u1ind=-1; u2ind=-1; d0ind=-1;
        FATAL("uninitliazed u0ind,u1ind,u2ind,d0ind in assignUMBplane!");
    }
    RU0.set(_R[u0ind].x,_R[u0ind].y,_R[u0ind].z);
    RU1.set(_R[u1ind].x,_R[u1ind].y,_R[u1ind].z);    
    RU2.set(_R[u2ind].x,_R[u2ind].y,_R[u2ind].z);    
    RD0.set(_R[d0ind].x,_R[d0ind].y,_R[d0ind].z);

    dR1=RU1-RU0;dR2=RU2-RU0;
    NV=cross(dR1,dR2);
    NN=NV/NV.norm();
  
    dR=RD0-RU0;
    distance=dot(NN,dR);
    if (distance < 0 ) NN=NN*(-1.0);
    UMB_R0 = RU0;
    UMB_nvec = NN;
    UMB_thickness = fabs(distance);

    if ( UMB_order_param == 41 )
	UMB_slipvec.set(-1.0,0.0,0.0);
    else if ( UMB_order_param == 42 )
    {
        dR.set(-1.0/sqrt(2.0),1.0/sqrt(2.0),0);
        dR1=cross(UMB_nvec,dR);
	UMB_slipvec=dR1;
    }
    else if ( UMB_order_param == 43 )
    {
	UMB_slipvec.set(4.0,-2.0,-1.0);
        UMB_slipvec=UMB_slipvec/UMB_slipvec.norm();
    }
    UMB_HMIN = 0.5;
    UMB_HMAX = 1.5;
    UMB_NNBR = 3;
   }
   else if ( UMB_order_param == 45 ) {
        INFO("set R0, nvec, thickness, slipvec from input");       
        UMB_R0.set  ( input[0], input[1], input[2] );
        UMB_nvec.set( input[3], input[4], input[5] );
        UMB_thickness = input[6];
        UMB_slipvec.set( input[7], input[8], input[9] );
        UMB_order_group=102506;
        UMB_HMIN = input[10];
        UMB_HMAX = input[11];
        UMB_NNBR = (int)input[12];
   } 
   else {
        FATAL("unknown UMB_order_param! (try 45)");
   }

   INFO_Printf("UMB_R0        = (%g, %g, %g)\n", UMB_R0.x, UMB_R0.y, UMB_R0.z);
   INFO_Printf("UMB_nvec      = (%g, %g, %g)\n", UMB_nvec.x, UMB_nvec.y, UMB_nvec.z);
   INFO_Printf("UMB_thickness =  %g \n", UMB_thickness);
   INFO_Printf("UMB_slipvec   = (%g, %g, %g)\n", UMB_slipvec.x, UMB_slipvec.y, UMB_slipvec.z);

    for(int i=0;i<_NP;i++) 
    {
        dR = _R[i]-UMB_R0; /* note this _R is computed from _SR2 */
	distance=dot(UMB_nvec,dR);
        if ( fabs(distance/UMB_thickness) < 0.5 ){ 
            group[i]=UMB_order_group; 
        //    INFO_Printf("%d goes into group %d\n", i, UMB_order_group);
        }
    }

    SHtoR(); /* now return R to H*_SR */

}


void MDFrame::HMCprocess()
{
    double tmp,R, Hamil_now, Hamil_pre;
    Hamil_pre=EPOT_temp+KATOM_temp;
    Hamil_now=_EPOT+_KATOM;
    tmp=( (Hamil_pre) - (Hamil_now) )/(KB*_TDES);
    R=exp(tmp);
INFO_Printf("EPOT_pre=%e, EPOT_now=%e\n",EPOT_temp,_EPOT);
INFO_Printf("KATOM_pre=%e, KATOM_now=%e\n",KATOM_temp,_KATOM);
INFO_Printf("T_HMC    =%f, T_inst   =%f\n",T_HMC,_T);
INFO_Printf("Hamil_pre=%e, Hamil_now=%e, R=%e\n",Hamil_pre,_HELM,R);
    if (tmp<-40) acc_HMC=0;
    else
    {
	if (drand48() < R)
	    acc_HMC=1;
	else
	    acc_HMC=0;
    }
    if (acc_HMC)
    {
	INFO("HMC accepted");
        if (YES_UMB==0)
	{
	    copyStoS4();
            EPOT_temp=_EPOT;
            T_HMC=_TDES;
//            perturbevelocity();
            KATOM_temp=_KATOM;
	}
        else if (YES_UMB==1)
	{
	    UMBprocess();
	    if (acc_UMB || UMB_K_in_Kelvin==-1 )
	    {
		EPOT_temp=_EPOT;
                INFO("HMC&UMB accepted");
                T_HMC=_TDES;
	    }
//            perturbevelocity();
	    KATOM_temp=_KATOM;
	}
    }
    else
    {
   	INFO("HMC rejected");
	copyS4toS();
//        perturbevelocity();
//        initvelocity();
        if (YES_UMB==1)
	{
	N_lgst_cluster=N_lgst_temp;
	UMBprocess();
	}
        KATOM_temp=_KATOM;
    }
}

void MDFrame::initHMC()
{
    setconfig4();
    copyStoS4();
    call_potential();
    calcprop();
    EPOT_temp=_EPOT;
    KATOM_temp=_KATOM;
    T_HMC=_TDES;
}

void MDFrame::setconfig4()
{
    Realloc(_SR4,Vector3,_NP);
    Realloc(_R4,Vector3,_NP);
    Realloc(_VSR4,Vector3,_NP);
    Realloc(_VR4,Vector3,_NP);
}

void MDFrame::copyStoS4()
{
    memcpy(_SR4,_SR,sizeof(Vector3)*_NP);
    memcpy(_R4,_R,sizeof(Vector3)*_NP);
    memcpy(_VSR4,_VSR,sizeof(Vector3)*_NP);
    memcpy(_VR4,_VR,sizeof(Vector3)*_NP);
    _H4=_H;
}

void MDFrame::copyS4toS()
{
    memcpy(_SR,_SR4,sizeof(Vector3)*_NP);
    memcpy(_R,_R4,sizeof(Vector3)*_NP);
    memcpy(_VSR,_VSR4,sizeof(Vector3)*_NP);
    memcpy(_VR,_VR4,sizeof(Vector3)*_NP);
    _H=_H4;
}

