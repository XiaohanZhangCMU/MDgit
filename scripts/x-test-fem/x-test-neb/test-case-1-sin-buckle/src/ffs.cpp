/*
  ffs.cpp
  by Seunghwa Ryu, shryu@stanford.edu
  Last Modified : Mon Sep  8 17:16:38 2008

  FUNCTION  :  Foward Flux Sampling for rare transition events

  To Do:
      1.  calcrystalorder() needs to be improved

  Note:
      This file is included by md.cpp.  This file cannot be compiled by itself.
*/

/********************************************************************/
/* Forward Flux Sampling (FFS) function (Seunghwa Ryu)              */
/********************************************************************/

void MDFrame::assign_Lam()
{
    int i, n;
    n = (int)input[0];
    Realloc(_Lam_array,int,n);
    Realloc(_Lam_check,int,n);
    for (i=1;i<=n;i++)
    {
	_Lam_array[i-1]=(int)input[i];
	_Lam_check[i-1]=0;
    }
    FFS0_check=1;
}

int MDFrame::FFSprocess()
{
    int i,Ntemp;
    double reaction;
    if (FFShist==1) logFFShist();
    if (FFSautoend==1) {
        if (FFSbackward==0){ 
	Ntemp = (int) _Lam_array[FFSoption*2];
	_Lam_array[FFSoption]=Ntemp;
        } else if ( FFSbackward==1 ) {
        Ntemp = (int) _Lam_array[0];
        _Lam_array[FFSoption]=Ntemp;
        }
    } 
    if (YES_SEP<2)
        NORM_REACTION=(1.0)*N_lgst_cluster/_Lam_array[FFSoption];
    else
	NORM_REACTION=(1.0-SEPA_RATIO)*N_lgst_cluster/_Lam_array[FFSoption]+
		SEPA_RATIO*SEPA_ORDER;
    reaction=NORM_REACTION;

 if (FFScommitor==0) {

  if (FFSbackward==0) {
    if (FFSoption==0) 
    {
	if (reaction >= 1.0 && FFS0_check == 0 )
	{
	    FFScn.write(this,zipfiles,true);
	    FFS0_check=1;
	}
        if (YES_SEP!=2)
        {
	if (FFS0_check==1 && N_lgst_cluster <= lambda_A )
	    FFS0_check=0;
        }
        else
        {
        if (FFS0_check==1 && reaction <= 0.4 )
            FFS0_check=0;
        }
        return 0;
    }
    else
    {
	if (reaction >= 1.0)
	{
	    FFS0_check = 1;
	    FFScn.write(this,zipfiles,true);
	    return 1;
	}
	
	if ( FFSpruning > 0 )
	{
	    for (i=FFSoption;i>1;i--)
	    {
		if ( N_lgst_cluster<=_Lam_array[i-2] && _Lam_check[i-2]==0 )
		{
		    if (drand48()>FFS_Pp)
		    {
			FFScn_weight=FFScn_weight/(1-FFS_Pp);
			_Lam_check[i-2]=1;
		    }
		    else
		    {
		 	FFS0_check = -1;
		    }
		}
		if ( N_lgst_cluster <= lambda_A ) FFS0_check=-1;
	    }
	    if ( N_lgst_cluster <= lambda_A ) FFS0_check=-1;
	}
	else
	{
	    if ( FFScurstep > step0 + totalsteps || N_lgst_cluster <= lambda_A )
		FFS0_check = -1;
	}
	if (FFS0_check<0) return 1;
	else return 0;
    }
} 
else if (FFSbackward==1) {
    if (FFSoption==0) 
    {
	if (N_lgst_cluster < B_lambda_0 && FFS0_check == 0 )
        {
	   FFScn.write(this,zipfiles,true);
	   FFS0_check=1;
 	}
        if (FFS0_check==1 && N_lgst_cluster >= B_lambda_B )
        {
	   FFS0_check=0;
	}
        if (N_lgst_cluster <= B_lambda_cut ) 
            return 0;
        else 
            return 1; 
    }
    else {
	if (N_lgst_cluster < _Lam_array[FFSoption])
        {
	    FFS0_check = 1;
            FFScn.write(this,zipfiles,true);
            return 1;
        }
        if ( FFScurstep > step0 + totalsteps || N_lgst_cluster >= B_lambda_B )
	    FFS0_check = -1;
        if (FFS0_check<0) return 1;
        else return 0;
    }
  }
 } else {
  if ( FFScurstep > step0 + totalsteps ) return 1;

  if (N_lgst_cluster <= lambda_A)
  {
        FFS0_check = -1;
	return 1;
  } else if (N_lgst_cluster >= B_lambda_B) {
        FFS0_check = 1;
	return 1;
  } else {
	FFS0_check = 0;
	return 0;
  }
 }

 WARNING("FFSProcess: should have returned before this line!");
 return -1;
}

void MDFrame::initFFShist()
{
    FILE *fp,*fp0;
    int i, no, temp;

    no = (int) 1.0*(MAX_LAM-MIN_LAM)/Delta_LAM+2;
    LAM_MAX_INDEX=no-1;
    Realloc(FFSfullarray,int,no);
    Realloc(FFSfullhistogram,int,no);

    fp0 = fopen("FFSfullhistogram.txt","r+");
    for (i=0;i<no;i++)
    {
	FFSfullarray[i]=MIN_LAM+Delta_LAM*i;
        FFSfullhistogram[i]=0;
    }
    if (fp0 == NULL)
    { 
    fp = fopen("FFSfullarray.txt","w");
    for (i=0;i<=LAM_MAX_INDEX;i++)
        fprintf(fp,"%d\n",FFSfullarray[i]);
    fclose(fp);
    }
    else
    {
        fclose(fp0);
        fp = fopen("FFSfullhistogram.txt","r");
        for (i=0;i<=LAM_MAX_INDEX;i++)
        {
	    fscanf(fp,"%d",&temp);
            INFO_Printf("%d",temp);
	    FFSfullhistogram[i]=FFSfullhistogram[i]+temp;
        }
        fclose(fp);
    }
}

void MDFrame::logFFShist()
{
    int i;
    if ( N_lgst_cluster<=FFSfullarray[0] )
	FFSfullhistogram[0]=FFSfullhistogram[0]+1;
    else if (N_lgst_cluster>FFSfullarray[LAM_MAX_INDEX-1]+Delta_LAM)
	FFSfullhistogram[LAM_MAX_INDEX]=FFSfullhistogram[LAM_MAX_INDEX]+1;
    else
    {
	i=(int) floor( 1.0*(N_lgst_cluster-MIN_LAM)/Delta_LAM);
        FFSfullhistogram[i]=FFSfullhistogram[i]+1;	
    }
}

void MDFrame::printFFShist()
{
    FILE *fp;
    int i;
    fp=fopen("FFSfullhistogram.txt","w");
    for (i=0;i<=LAM_MAX_INDEX;i++)
	fprintf(fp,"%d\n",FFSfullhistogram[i]);
    fclose(fp);
}

void MDFrame::initKinetic()
{
    allocKinetic();
    kinetic_flag = 0;
}


void MDFrame::allocKinetic()
{
    int i;
    kineticNmax = (int) (1.0*Kinetic_Time/savepropfreq+1.0);
    Realloc(time_array,double,kineticNmax+1);
    Realloc(n_array,int,kineticNmax+1);
    for (i=0;i<kineticNmax+1;i++)
    {
        time_array[i]=1.0*i*savepropfreq*_TIMESTEP;
        n_array[i]=0;
    }
}


int MDFrame::Kineticprocess()
{
    if ( kinetic0 > kineticNmax )
  {
        kinetic0 = 0;
        return 1;
   }
  else
  {
    if ( Kinetic == 1 && kinetic_flag == 0 ) {
        KN_Center = N_lgst_cluster;
        kinetic_flag = 1;
        kinetic0=0;
    } else if ( Kinetic == 1 && kinetic_flag == 1 && curstep%(savepropfreq)==0 ) {
        kinetic0++;
        INFO_Printf("kinetic0=%d Nlgst=%d KNCenter=%d\n",kinetic0,N_lgst_cluster,KN_Center);
        n_array[kinetic0]=n_array[kinetic0]+(N_lgst_cluster-KN_Center);
    }
  return 0;
  }
}


void MDFrame::kineticprint()
{
    FILE *fp, *fp1, *fp2;
    int i;
    fp = fopen("time_array.txt","w");
    fp1= fopen("n_array.txt","a+");
    fp2= fopen("KN_LOG.txt","a+");
    for (i=0;i< kineticNmax + 1;i++) {
        fprintf(fp,"%10.4e ",time_array[i]);
        fprintf(fp1,"%5d ",n_array[i]);
    }
        fprintf(fp,"\n");
        fprintf(fp1,"\n");
        fprintf(fp2,"%d\n",KN_Center);
    fclose(fp);fclose(fp1);fclose(fp2);
}

