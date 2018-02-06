/*
  crystal_order.cpp 
  taken from ffs.cpp by Seunghwa Ryu, shryu@stanford.edu
  Last Modified : Sun Aug  5 17:16:38 2012

  FUNCTION  :  compute crystalline order parameters for atoms

  To Do:
      1.  calcrystalorder() needs to be improved

  Note:
      This file is included by md.cpp.  This file cannot be compiled by itself.
*/

void MDFrame::calcrystalorder()
{
    int i, j, jpt, m, M, Nb,Nc_i;
    Vector3 sij, rij;
    double rrij, xij, yij, zij, r_l, rc,ro;
    double metric[30];

    int N_cluster, Size_cluster[_NP];    
    int k0, N_IMP;

    M = 2*L_for_QLM + 2;
    N_solid_P = 0;
    N_skin_P = 0;
    N_lgst_cluster = 0;
    N_lgst_skin = 0;
    N_cluster = 0;
    N_IMP = 0;
    Nc_i = 0;

    if (L_for_QLM==3)
    {
        metric[0]=5.0;
	metric[1]=metric[0];
	metric[2]=30.0;
	metric[3]=metric[2];
	metric[4]=3.0;
	metric[5]=metric[4];
	metric[6]=2.0;
    }
    else if (L_for_QLM==4)
    {
        metric[0] = 35.0;
        metric[1] = metric[0];
        metric[2] = 280.0;
        metric[3] = metric[2];
        metric[4] = 20.0;
        metric[5] = metric[4];
        metric[6] = 40.0;
        metric[7] = metric[6];
        metric[8] = 2.0;
    }
     else if (L_for_QLM==6)
    {
	metric[0] = 231.0;
	metric[1] = metric[0];
	metric[2] = 2772.0;
	metric[3] = metric[2];
	metric[4] = 126.0;
	metric[5] = metric[4];
	metric[6] = 420.0;
	metric[7] = metric[6];
	metric[8] = 105.0;
	metric[9] = metric[8];
	metric[10]= 168.0;
	metric[11]= metric[10];
	metric[12]= 2.0;
    }
    else
    {
        ERROR("calcrystalorder(): L = "<<L_for_QLM<<" is not implemented");
        return;
    }

    /* initialize QLM and QL to zero */
    for(i=0 ; i < _NP ; i++)
    {
        if (species[i]==IMP_INDEX) N_IMP++;
	QLM_mark[i]=0;
        QLM_solid[i]=0;
	Size_cluster[i]=0;
        _QL[i]=0;
        for(m=0;m<M;m++)
            _QLM[i*M+m]=0;
    }

    for(i=0;i<_NP;i++)
    {
        Nb=0;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(jpt==i) continue;
            sij=_SR[jpt]-_SR[i]; sij.subint();
            rij=_H*sij; rrij=rij.norm();

            if(rrij > Rc_for_QLM) continue;
	    xij=rij.x;
	    yij=rij.y;
	    zij=rij.z;
	    r_l=pow(rrij,(double)L_for_QLM);

            /* rij.x, rij.y, rij.z
             * rij[0], rij[1], rij[2]
             */
	    Nb++;
	    if (L_for_QLM == 3)
            {	
                 _QLM[i*M+ 0 ] += ( pow(xij,3.0) - 3.0*xij*pow(yij,2.0) )/r_l;
		 _QLM[i*M+ 1 ] += ( pow(yij,3.0) - 3.0*yij*pow(xij,2.0) )/r_l;
		 _QLM[i*M+ 2 ] += ( pow(xij,2.0) - pow(yij,2.0) )*zij/r_l;
		 _QLM[i*M+ 3 ] += -2.0*xij*yij*zij/r_l;
		 _QLM[i*M+ 4 ] += (4.0*zij*zij-xij*xij-yij*yij)*xij/r_l;
		 _QLM[i*M+ 5 ] += -(4.0*zij*zij-xij*xij-yij*yij)*yij/r_l;
		 _QLM[i*M+ 6 ] += zij*(2.0*zij*zij-3.0*xij*xij-3.0*yij*yij)/r_l;            
	    }
            else if (L_for_QLM ==4)
            {
                _QLM[i*M+ 0 ] += ( pow(xij,4.0) - 6.0*pow(xij,2.0)*pow(yij,2.0) + pow(yij,4.0) )/r_l;
                _QLM[i*M+ 1 ] += ( -4.0*xij*yij*( xij*xij-yij*yij ) )/r_l;
                _QLM[i*M+ 2 ] += ( xij*zij*( xij*xij - 3.0*yij*yij) )/r_l;
                _QLM[i*M+ 3 ] += ( -1.0*yij*zij*(3.0*xij*xij - yij*yij) )/r_l;
                _QLM[i*M+ 4 ] += ( (yij*yij-xij*xij)*(xij*xij+yij*yij-6.0*zij*zij) )/r_l;
                _QLM[i*M+ 5 ] += ( 2.0*xij*yij*( xij*xij + yij*yij - 6.0*zij*zij) )/r_l;
                _QLM[i*M+ 6 ] += ( -1.0*xij*zij*(3.0*xij*xij+3.0*yij*yij-4.0*zij*zij) )/r_l;
                _QLM[i*M+ 7 ] += ( yij*zij*( 3.0*xij*xij + 3.0*yij*yij - 4.0*zij*zij) )/r_l;
                _QLM[i*M+ 8 ] += ( 35.0*pow(zij,4.0)-30.0*zij*zij*(xij*xij+yij*yij+zij*zij)+3.0*r_l )/r_l;
            }
            else if (L_for_QLM ==6)
	    {
		_QLM[i*M+ 0 ] += (pow(xij,6.0)-15*pow(xij,4.0)*pow(yij,2.0)+15*pow(xij,2.0)*pow(yij,4.0)-pow(yij,6.0))/r_l;
		_QLM[i*M+ 1 ] += 2*xij*yij*(3*pow(xij,4.0)-10*pow(xij,2.0)*pow(yij,2.0)+3*pow(yij,4.0))/r_l;
		_QLM[i*M+ 2 ] += xij*zij*(pow(xij,4.0)-10*pow(xij,2.0)*pow(yij,2.0)+5*pow(yij,4.0))/r_l;
		_QLM[i*M+ 3 ] += yij*zij*(5*pow(xij,4.0)-10*pow(xij,2.0)*pow(yij,2.0)+pow(yij,4.0))/r_l;
		_QLM[i*M+ 4 ] += -(-10*pow(xij,4.0)*pow(zij,2.0)+pow(xij,6.0)-5*pow(xij,4.0)*pow(yij,2.0)+60*pow(xij*yij*zij,2.0)-5*pow(xij,2.0)*pow(yij,4.0)-10*pow(yij,4.0)*pow(zij,2.0)+pow(yij,6.0))/r_l;
		_QLM[i*M+ 5 ] += -4*(-10*pow(zij,2.0)+xij*xij+yij*yij)*xij*yij*(xij*xij-yij*yij)/r_l;
		_QLM[i*M+ 6 ] += -xij*zij*(-8*pow(xij,2.0)*pow(zij,2.0)+3*pow(xij,4.0)-6*pow(xij,2.0)*pow(yij,2.0)+24*pow(yij,2.0)*pow(zij,2.0)-9*pow(yij,4.0))/r_l;
		_QLM[i*M+ 7 ] += -zij*(-8*pow(zij,2.0)+3*xij*xij+3*yij*yij)*yij*(3*xij*xij-yij*yij)/r_l;
		_QLM[i*M+ 8 ] += (16*pow(xij,2.0)*pow(zij,4.0)-16*pow(xij,4.0)*pow(zij,2.0)+pow(xij,6.0)+pow(xij,4.0)*pow(yij,2.0)-pow(xij,2.0)*pow(yij,4.0)-16*pow(yij,2.0)*pow(zij,4.0)+16*pow(yij,4.0)*pow(zij,2.0)-pow(yij,6.0))/r_l;
		_QLM[i*M+ 9 ] += 2*(16*pow(zij,4.0)-16*pow(xij,2.0)*pow(zij,2.0)-16*pow(yij,2.0)*pow(zij,2.0)+pow(xij,4.0)+2*pow(xij*yij,2.0)+pow(yij,4.0))*xij*yij/r_l;
		_QLM[i*M+ 10] += xij*zij*(8*pow(zij,4.0)-20*pow(xij*zij,2.0)-20*pow(yij*zij,2.0)+5*pow(xij,4.0)+10*pow(xij*yij,2.0)+5*pow(yij,4.0))/r_l;
		_QLM[i*M+ 11] += yij*zij*(8*pow(zij,4.0)-20*pow(xij*zij,2.0)-20*pow(yij*zij,2.0)+5*pow(xij,4.0)+10*pow(xij*yij,2.0)+5*pow(yij,4.0))/r_l;
		_QLM[i*M+ 12] += -(-16*pow(zij,6.0)+120*pow(xij,2.0)*pow(zij,4.0)+120*pow(yij,2.0)*pow(zij,4.0)-90*pow(xij,4.0)*pow(zij,2.0)-180*pow(xij*yij*zij,2.0)-90*pow(yij,4.0)*pow(zij,2.0)+5*pow(xij,6.0)+15*pow(xij,4.0)*pow(yij,2.0)+15*pow(xij,2.0)*pow(yij,4.0)+5*pow(yij,6.0))/r_l;
	    }
        }
        /* the last entry is the magnitude of the vector */
	for (m=0;m<M-1;m++)
	    _QLM[i*M+M-1] += metric[m]*_QLM[i*M+m]*_QLM[i*M+m];
	_QLM[i*M+M-1]=sqrt(_QLM[i*M+M-1]);
    }

    for(i=0;i<_NP;i++)
    {
        Nb=0;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(jpt==i) continue;
            sij=_SR[jpt]-_SR[i]; sij.subint();
            rij=_H*sij; rrij=rij.norm();  
            if(rrij > Rc_for_QLM) continue;
	    Nb++; 
            for(m=0;m<M-1;m++)
            {
                if ( fabs(_QLM[jpt*M+M-1]) > 0.0000001 )
                    _QL[i] += _QLM[i*M+m] * _QLM[jpt*M+m] * metric[m] / _QLM[jpt*M+M-1];
            }
        }
        if (Nb>0)
	{
	_QL[i]/= fabs((_QLM[i*M+M-1]*Nb));
		if (L_for_QLM==6)
			_QL[i]=-_QL[i];
	}
        else
        _QL[i] = 2; 

        /* put results in TOPOL array for visualization */
	_TOPOL[i] = _QL[i];
	if (_QL[i] <= QLM_cutoff )
	{
		N_solid_P++;
		QLM_solid[i]=1;
	}
    }

    k0=0;
    for (i=0;i<_NP;i++)
    {
	if ( QLM_mark[i]==0 && QLM_solid[i]==1 )
	{
	    N_cluster++;
	    N_cluster_temp = 0;
	    if (N_cluster == 1)
		Nc_i=i;
 	    if (wSKIN>0) DFS1(i,0,0,0);
	    else DFS(i,0,0);
	    Size_cluster[N_cluster]=N_cluster_temp;
            if (N_cluster_temp > N_lgst_cluster)
            {
		N_lgst_cluster = N_cluster_temp;
 	        Nc_i=i;
            }
//	    printf("Size_cluster = %d\n", N_cluster_temp);
//	    printf("N_largest_cluster = %d\n", N_lgst_cluster);	    
/*********************************************/
	    if (YES_UMB==1 && UMB_K_in_Kelvin < 0 && N_cluster_temp > 0 )
	    {
		for (j=0;j<=2*delta_n;j++) if (Narray[j]==N_cluster_temp) {
		if ( UMB_curstep > UMB_equilstep ) {
		    Occurence[j]++;
		    Probability[j]++;
		}
		}
	    }
/*********************************************/
	}
}
/*********************************************/    
    if ( YES_UMB == 1 && UMB_K_in_Kelvin < 0 )
    {   
        k0=_NP-N_solid_P;
	for (j=0;j<=2*delta_n;j++) if (Narray[j]==0) {
	if (UMB_curstep > UMB_equilstep ) {
	    Occurence[j]=Occurence[j]+k0;
	    Probability[j]=Probability[j]+k0;
	}
	} 
    }
/*********************************************/    
    _CLUSTER_CM.clear();
    if ( N_lgst_index != 0 )
    {
	for (i=0;i<_NP;i++)
	    QLM_mark[i]=0;
        _SRforCM=_SR[Nc_i];
	if (wSKIN >0 )
	    DFS1(Nc_i,1,N_lgst_index,N_lgst_skin_index);
	else
	    DFS(Nc_i,1,N_lgst_index);
    }
    if (YES_SEP==2)
    {
	ro=pow(3.0*_Lam_array[FFSoption]*_OMEGA/4.0/M_PI/_NP,1.0/3.0)+0.5*Rc_for_QLM;
        rc=pow(3.0*_OMEGA/4.0/M_PI,1.0/3.0);
        SEPA_TARGET=4.0*M_PI*N_IMP/_OMEGA/(1.0-pow(ro/rc,3.0))*
                    (pow(rc,3.0-IMP_R_EXP)-pow(ro,3.0-IMP_R_EXP))/(3.0-IMP_R_EXP);
    }
    if (YES_SEP>0) 
    {
       _CLUSTER_CM=_CLUSTER_CM/N_IMP;
       _CLUSTER_CM=_CLUSTER_CM+_SRforCM;
       SEPA_ORDER=0.0;
       for (i=0;i<_NP;i++)
       {
   	   if (species[i]!=IMP_INDEX) continue;
           if (IMP_TOPOL_INDEX!=0)
           {
              _QL[i]=IMP_TOPOL_INDEX;
              _TOPOL[i]=IMP_TOPOL_INDEX;
           }
           if (species[Nc_i]==IMP_INDEX) continue;
           sij=_SR[i]-_CLUSTER_CM;
           sij.subint();
           rij=_H*sij;rrij=rij.norm();
           rrij=pow(rrij,IMP_R_EXP);
           SEPA_ORDER=SEPA_ORDER+1.0/rrij;
       }
       if (SEPA_ORDER>0.00001)
       {
	   SEPA_ORDER=1.0/(1.0*SEPA_ORDER/SEPA_TARGET);
       }
       NC_TIMES_SEPA = 1.0*N_lgst_cluster*SEPA_ORDER;
       if (YES_SEP==1) N_lgst_cluster = (int) NC_TIMES_SEPA;
    }
}



#ifdef _GSL
void MDFrame::calqlwl()
{
    int i, j, jpt, m, M, Nb,Nc_i, m1, m2, m3, o,p,q;
    Vector3 sij, rij;
    double rrij, xij, yij, zij, r_l, rc,ro;
    double metric[30];

    int N_cluster, Size_cluster[_NP];    
    int k0, N_IMP;
//    INFO("\n\n\n calqlwl!!!  \n\n\n");

    if ( UMB_order_param >= 110 && UMB_order_param < 120 ) 
    { 
	calcrystalorder();
    }

    M = 2*l_for_qlm + 2;
    N_solid_P = 0;
    N_skin_P = 0;
    N_lgst_cluster = 0;
    N_lgst_skin = 0;
    N_cluster = 0;
    N_IMP = 0;
    Nc_i = 0;
//	printf ("debug\n");
    if (l_for_qlm==3)
    {
        metric[0]=sqrt(35.0/M_PI)/8.0;
	metric[1]=metric[0];
	metric[2]=sqrt(105.0/2.0/M_PI)/4.0;
	metric[3]=metric[2];
	metric[4]=sqrt(21.0/M_PI)/8.0;
	metric[5]=metric[4];
	metric[6]=sqrt(7.0/M_PI)/4.0;
    }
    else if (l_for_qlm==4)
    {
        metric[0] = sqrt(35.0/M_PI/2.0)/16.0*3.0;
        metric[1] = metric[0];
        metric[2] = sqrt(35.0/M_PI)/8.0*3.0;
        metric[3] = metric[2];
        metric[4] = sqrt(5.0/M_PI/2.0)/8.0*3.0;
        metric[5] = metric[4];
        metric[6] = sqrt(5.0/M_PI)/8.0*3.0;
        metric[7] = metric[6];
        metric[8] = sqrt(1.0/M_PI)/16.0*3.0;
    }
     else if (l_for_qlm==6)
    {
	metric[0] = sqrt(3003.0/M_PI)/64.0;
	metric[1] = metric[0];
	metric[2] = sqrt(1001.0/M_PI)/32.0*3.0;
	metric[3] = metric[2];
	metric[4] = sqrt(91.0/M_PI/2.0)/32.0*3.0;
	metric[5] = metric[4];
	metric[6] = sqrt(1365.0/M_PI)/32.0;
	metric[7] = metric[6];
	metric[8] = sqrt(1365.0/M_PI)/64.0;
	metric[9] = metric[8];
	metric[10]= sqrt(273.0/M_PI/2.0)/16.0;
	metric[11]= metric[10];
	metric[12]= sqrt(13.0/M_PI)/32.0;
    }
    else
    {
        ERROR("calcrystalorder(): L = "<<l_for_qlm<<" is not implemented");
        return;
    }

    /* initialize QLM and QL to zero */
    for(i=0 ; i < _NP ; i++)
    {
        if (species[i]==IMP_INDEX) N_IMP++;
	QLM_mark[i]=0;
        QLM_solid[i]=0;
	Size_cluster[i]=0;
        _ql[i]=0;
        _wl[i]=0;
        _ql_[i]=0;
        _wl_[i]=0;
//        _QL[i]=0;
        _Ql[i]=0;
        _TOPOL[i]=0;
        for(m=0;m<M;m++)
        { 
	    _qlmr[i*M+m]=0;
            _qlmi[i*M+m]=0;
            _qlm_r[i*M+m]=0;
            _qlm_i[i*M+m]=0;
            _qlm[i*M+m]=0;
        }
    }
//	INFO("calcry1.5\n");
    for(i=0;i<_NP;i++)
    {
        Nb=0;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(jpt==i) continue;
            sij=_SR[jpt]-_SR[i]; sij.subint();
            rij=_H*sij; rrij=rij.norm();

            if(rrij > Rc_for_QLM) continue;
	    xij=rij.x;
	    yij=rij.y;
	    zij=rij.z;
	    r_l=pow(rrij,(double)l_for_qlm);

            /* rij.x, rij.y, rij.z
             * rij[0], rij[1], rij[2]
             */
	    Nb++;
	    if (l_for_qlm == 3)
            {	
                 _qlm[i*M+ 0 ] += ( pow(xij,3.0) - 3.0*xij*pow(yij,2.0) )/r_l;
		 _qlm[i*M+ 1 ] += ( pow(yij,3.0) - 3.0*yij*pow(xij,2.0) )/r_l;
		 _qlm[i*M+ 2 ] += ( pow(xij,2.0) - pow(yij,2.0) )*zij/r_l;
		 _qlm[i*M+ 3 ] += -2.0*xij*yij*zij/r_l;
		 _qlm[i*M+ 4 ] += (4.0*zij*zij-xij*xij-yij*yij)*xij/r_l;
		 _qlm[i*M+ 5 ] += -(4.0*zij*zij-xij*xij-yij*yij)*yij/r_l;
		 _qlm[i*M+ 6 ] += zij*(2.0*zij*zij-3.0*xij*xij-3.0*yij*yij)/r_l;       
	    }
            else if (l_for_qlm ==4)
            {
                _qlm[i*M+ 0 ] += ( pow(xij,4.0) - 6.0*pow(xij,2.0)*pow(yij,2.0) + pow(yij,4.0) )/r_l;
                _qlm[i*M+ 1 ] += ( -4.0*xij*yij*( xij*xij-yij*yij ) )/r_l;
                _qlm[i*M+ 2 ] += ( xij*zij*( xij*xij - 3.0*yij*yij) )/r_l;
                _qlm[i*M+ 3 ] += ( -1.0*yij*zij*(3.0*xij*xij - yij*yij) )/r_l;
                _qlm[i*M+ 4 ] += ( (yij*yij-xij*xij)*(xij*xij+yij*yij-6.0*zij*zij) )/r_l;
                _qlm[i*M+ 5 ] += ( 2.0*xij*yij*( xij*xij + yij*yij - 6.0*zij*zij) )/r_l;
                _qlm[i*M+ 6 ] += ( -1.0*xij*zij*(3.0*xij*xij+3.0*yij*yij-4.0*zij*zij) )/r_l;
                _qlm[i*M+ 7 ] += ( yij*zij*( 3.0*xij*xij + 3.0*yij*yij - 4.0*zij*zij) )/r_l;
                _qlm[i*M+ 8 ] += ( 35.0*pow(zij,4.0)-30.0*zij*zij*(xij*xij+yij*yij+zij*zij)+3.0*r_l )/r_l;
            }
            else if (l_for_qlm ==6)
	    {
		_qlm[i*M+ 0 ] += (pow(xij,6.0)-15*pow(xij,4.0)*pow(yij,2.0)+15*pow(xij,2.0)*pow(yij,4.0)-pow(yij,6.0))/r_l;
		_qlm[i*M+ 1 ] += 2*xij*yij*(3*pow(xij,4.0)-10*pow(xij,2.0)*pow(yij,2.0)+3*pow(yij,4.0))/r_l;
		_qlm[i*M+ 2 ] += xij*zij*(pow(xij,4.0)-10*pow(xij,2.0)*pow(yij,2.0)+5*pow(yij,4.0))/r_l;
		_qlm[i*M+ 3 ] += yij*zij*(5*pow(xij,4.0)-10*pow(xij,2.0)*pow(yij,2.0)+pow(yij,4.0))/r_l;
		_qlm[i*M+ 4 ] += -(-10*pow(xij,4.0)*pow(zij,2.0)+pow(xij,6.0)-5*pow(xij,4.0)*pow(yij,2.0)+60*pow(xij*yij*zij,2.0)-5*pow(xij,2.0)*pow(yij,4.0)-10*pow(yij,4.0)*pow(zij,2.0)+pow(yij,6.0))/r_l;
		_qlm[i*M+ 5 ] += -4*(-10*pow(zij,2.0)+xij*xij+yij*yij)*xij*yij*(xij*xij-yij*yij)/r_l;
		_qlm[i*M+ 6 ] += -xij*zij*(-8*pow(xij,2.0)*pow(zij,2.0)+3*pow(xij,4.0)-6*pow(xij,2.0)*pow(yij,2.0)+24*pow(yij,2.0)*pow(zij,2.0)-9*pow(yij,4.0))/r_l;
		_qlm[i*M+ 7 ] += -zij*(-8*pow(zij,2.0)+3*xij*xij+3*yij*yij)*yij*(3*xij*xij-yij*yij)/r_l;
		_qlm[i*M+ 8 ] += (16*pow(xij,2.0)*pow(zij,4.0)-16*pow(xij,4.0)*pow(zij,2.0)+pow(xij,6.0)+pow(xij,4.0)*pow(yij,2.0)-pow(xij,2.0)*pow(yij,4.0)-16*pow(yij,2.0)*pow(zij,4.0)+16*pow(yij,4.0)*pow(zij,2.0)-pow(yij,6.0))/r_l;
		_qlm[i*M+ 9 ] += 2*(16*pow(zij,4.0)-16*pow(xij,2.0)*pow(zij,2.0)-16*pow(yij,2.0)*pow(zij,2.0)+pow(xij,4.0)+2*pow(xij*yij,2.0)+pow(yij,4.0))*xij*yij/r_l;
		_qlm[i*M+ 10] += xij*zij*(8*pow(zij,4.0)-20*pow(xij*zij,2.0)-20*pow(yij*zij,2.0)+5*pow(xij,4.0)+10*pow(xij*yij,2.0)+5*pow(yij,4.0))/r_l;
		_qlm[i*M+ 11] += yij*zij*(8*pow(zij,4.0)-20*pow(xij*zij,2.0)-20*pow(yij*zij,2.0)+5*pow(xij,4.0)+10*pow(xij*yij,2.0)+5*pow(yij,4.0))/r_l;
		_qlm[i*M+ 12] += -(-16*pow(zij,6.0)+120*pow(xij,2.0)*pow(zij,4.0)+120*pow(yij,2.0)*pow(zij,4.0)-90*pow(xij,4.0)*pow(zij,2.0)-180*pow(xij*yij*zij,2.0)-90*pow(yij,4.0)*pow(zij,2.0)+5*pow(xij,6.0)+15*pow(xij,4.0)*pow(yij,2.0)+15*pow(xij,2.0)*pow(yij,4.0)+5*pow(yij,6.0))/r_l;
	    }
        }

/* assign qlmr (qlm real) and qlmi (qlm imaginary) */
        for (m=0;m<M/2;m++)
        {
	    _qlmr[i*M+m]=_qlm[i*M+2*m]*metric[2*m];
            _qlmi[i*M+m]=_qlm[i*M+2*m+1]*metric[2*m];
            _qlmr[i*M+M-2-m]=pow(-1.0,m-l_for_qlm)*_qlmr[i*M+m];
            _qlmi[i*M+M-2-m]=-1.0*pow(-1.0,m-l_for_qlm)*_qlmi[i*M+m];
        }

//if (i==0) INFO_Printf("Nb=%d, qlmr[i*M+M-1]=%f qlm_r[i*M+M-1]=%f",Nb,_qlmr[i*M+M-1],_qlm_r[i*M+M-1]); 

        if (Nb==0) continue;
/* assign sum |qlm|^2 to qlmr[2l+1] */
        for (m=0;m<M-1;m++) 
        {
	    _qlmr[i*M+m] /=(double)Nb;
            _qlmi[i*M+m] /=(double)Nb;
	    _qlmr[i*M+M-1] += _qlmr[i*M+m]*_qlmr[i*M+m] + _qlmi[i*M+m]*_qlmi[i*M+m];
        }

/* compute ql and assign to _ql */
        _ql[i]=sqrt( 4.0*M_PI/(2.0*l_for_qlm+1.0) * _qlmr[i*M+M-1] );

/* o=m1+l p=m2+l q=m3+l. compute wl assign to wl */
        for (o=0;o<=2*l_for_qlm;o++)
	{
	    for (p=0;p<=2*l_for_qlm;p++)
	    {
	        m1=o-l_for_qlm;
	        m2=p-l_for_qlm;
	        m3=-m1-m2;
	        if ( abs(m3) > l_for_qlm ) continue;
                q=m3+l_for_qlm;
               
                _wl[i] += gsl_sf_coupling_3j(2*l_for_qlm,2*l_for_qlm,2*l_for_qlm,2*m1,2*m2,2*m3)
                       *( _qlmr[i*M+o]*_qlmr[i*M+p]*_qlmr[i*M+q]-_qlmr[i*M+o]*_qlmi[i*M+p]*_qlmi[i*M+q]
                         -_qlmi[i*M+o]*_qlmr[i*M+p]*_qlmi[i*M+q]-_qlmi[i*M+o]*_qlmi[i*M+p]*_qlmr[i*M+q]);
	    }
	}
        _wl[i]/=pow(_qlmr[i*M+M-1],1.5);
    }


/* Now, compute average properties */
    for(i=0;i<_NP;i++)
    {
        Nb=0;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(jpt==i) continue;
            sij=_SR[jpt]-_SR[i]; sij.subint();
            rij=_H*sij; rrij=rij.norm();  
            if(rrij > Rc_for_QLM) continue;
	    Nb++; 

            for (m=0;m<M-1;m++) 
            {
	        _qlm_r[i*M+m] += _qlmr[jpt*M+m];
                _qlm_i[i*M+m] += _qlmi[jpt*M+m];
	    }

        }
        if (Nb==0) continue;

/* assign sum over |qlm_|^2 */
        for (m=0;m<M-1;m++) 
        {
	        _qlm_r[i*M+m] /= (double)Nb;
                _qlm_i[i*M+m] /= (double)Nb;
		_qlm_r[i*M+M-1] += _qlm_r[i*M+m]*_qlm_r[i*M+m] + _qlm_i[i*M+m]*_qlm_i[i*M+m];
	}

/* assign ql_ */
        _ql_[i]=sqrt( 4.0*M_PI/(2.0*l_for_qlm+1.0) * _qlm_r[i*M+M-1] );

/* assign wl_ */
        for (o=0;o<=2*l_for_qlm;o++)
	{
	    for (p=0;p<=2*l_for_qlm;p++)
	    {
	        m1=o-l_for_qlm;
	        m2=p-l_for_qlm;
	        m3=-m1-m2;
	        if ( abs(m3) > l_for_qlm ) continue;
                q=m3+l_for_qlm;
               
                _wl_[i] += gsl_sf_coupling_3j(2*l_for_qlm,2*l_for_qlm,2*l_for_qlm,2*m1,2*m2,2*m3)
                       *( _qlm_r[i*M+o]*_qlm_r[i*M+p]*_qlm_r[i*M+q]-_qlm_r[i*M+o]*_qlm_i[i*M+p]*_qlm_i[i*M+q]
                         -_qlm_i[i*M+o]*_qlm_r[i*M+p]*_qlm_i[i*M+q]-_qlm_i[i*M+o]*_qlm_i[i*M+p]*_qlm_r[i*M+q]);
	    }
	}
        _wl_[i]/=pow(_qlm_r[i*M+M-1],1.5);

        if (qlm_id==1) _TOPOL[i]=_ql[i];
        else if (qlm_id==2) _TOPOL[i]=_wl[i];
        else if (qlm_id==3) _TOPOL[i]=_ql_[i];
        else if (qlm_id==4) _TOPOL[i]=_wl_[i];
        else _TOPOL[i]=_ql[i];

 
        if (_Ql[i]<1.9)    
        {
           if ( _ql_[i]>0.05 && _wl_[i]> 0.00 && _QL[i] < QLM_cutoff ) 
           {
               _Ql[i]=-0.5; // HCP
//               INFO_Printf("%d HCP\n",i);
           }
           else if ( _ql_[i]>0.10 && _wl_[i]<-0.10 && _QL[i] < QLM_cutoff )
           {
               _Ql[i]=0.5; // FCC
//               INFO_Printf("%d FCC\n",i);
           }
           else
               _Ql[i]=0.0; // Liquid
        }


        if ( UMB_order_param >= 110 && UMB_order_param < 120 )
        {
	    _TOPOL[i]=_Ql[i];
        }

        if ( UMB_order_param == 111  )
	    if (_Ql[i]>-1.0 && _Ql[i]<-0.25)
            {
                N_solid_P++;
                QLM_solid[i]=1;
            }
        if ( UMB_order_param == 112  )
  	   if (_Ql[i]>0.25 && _Ql[i]<1.0)
           {
                N_solid_P++;
                QLM_solid[i]=1;
           }
        if ( UMB_order_param == 113  )
           if (fabs(_Ql[i])>0.25 && fabs(_Ql[i])<1.0)
           {
                N_solid_P++;
                QLM_solid[i]=1;
           }
        if ( UMB_order_param == 114  )
	   if (_Ql[i]>-0.25 && _Ql[i]<0.25)
           {
                N_solid_P++;
                QLM_solid[i]=1;
           }
//if (i==_NP-1) INFO_Printf("N_solid_P=%d\n",N_solid_P);

    }
 

//    INFO("calcry3");
    k0=0;
    for (i=0;i<_NP;i++)
    {
	if ( QLM_mark[i]==0 && QLM_solid[i]==1 )
	{
	    N_cluster++;
	    N_cluster_temp = 0;
	    if (N_cluster == 1)
		Nc_i=i;
 	    if (wSKIN>0) DFS1(i,0,0,0);
	    else DFS(i,0,0);
	    Size_cluster[N_cluster]=N_cluster_temp;
            if (N_cluster_temp > N_lgst_cluster)
            {
		N_lgst_cluster = N_cluster_temp;
 	        Nc_i=i;
            }
//	    printf("Size_cluster = %d\n", N_cluster_temp);
//	    printf("N_largest_cluster = %d\n", N_lgst_cluster);	    
/*********************************************/
	    if (YES_UMB==1 && UMB_K_in_Kelvin < 0 && N_cluster_temp > 0 )
	    {
		for (j=0;j<=2*delta_n;j++) if (Narray[j]==N_cluster_temp) {
		if ( UMB_curstep > UMB_equilstep ) {
		    Occurence[j]++;
		    Probability[j]++;
		}
		}
	    }
/*********************************************/
	}
    }
//INFO_Printf("N_lgst_cluster=%d\n",N_lgst_cluster);
/*********************************************/    
    if ( YES_UMB == 1 && UMB_K_in_Kelvin < 0 )
    {   
        k0=_NP-N_solid_P;
	for (j=0;j<=2*delta_n;j++) if (Narray[j]==0) {
	if (UMB_curstep > UMB_equilstep ) {
	    Occurence[j]=Occurence[j]+k0;
	    Probability[j]=Probability[j]+k0;
	}
	} 
    }
/*********************************************/    
    _CLUSTER_CM.clear();
    if ( N_lgst_index != 0 )
    {
	for (i=0;i<_NP;i++)
	    QLM_mark[i]=0;
        _SRforCM=_SR[Nc_i];
	if (wSKIN >0 )
	    DFS1(Nc_i,1,N_lgst_index,N_lgst_skin_index);
	else
	    DFS(Nc_i,1,N_lgst_index);
    }

    
   
}

void MDFrame::allocqlm()
{
    int size;

    if (l_for_qlm < 0)
    {
        ERROR("l_for_qlm = "<<l_for_qlm<<" skip allocqlm");
        return;
    }
    if (_NP <= 0)
    {
        ERROR("NP = "<<_NP<<" skip allocqlm");
        return;
    }


    size = (2*l_for_qlm + 2)*_NP;

    Realloc(_qlmr,double,size);
    Realloc(_qlmi,double,size);
    Realloc(_qlm_r,double,size);
    Realloc(_qlm_i,double,size);
    Realloc(_qlm,double,size);

    Realloc(_ql,double,_NP);
    Realloc(_ql_,double,_NP);
    Realloc(_wl,double,_NP);
    Realloc(_wl_,double,_NP);
    Realloc(_Ql,double,_NP);

    if (strcmp(qlm_type,"ql")==0)
        qlm_id = 1;
    else if (strcmp(qlm_type,"wl")==0)
        qlm_id = 2;
    else if (strcmp(qlm_type,"ql_")==0)
        qlm_id = 3;
    else if (strcmp(qlm_type,"wl_")==0)
        qlm_id = 4;
    else 
	qlm_id = 1;

}

#endif



void MDFrame::caldislocationorder()
{
    int i, Nc_i, j, jpt, k0, n;
    int N_cluster, Size_cluster[_NP],nn_rc, N_bulk;
    Vector3 sij0, sij, rij0, rij, ri00;
    Vector3 drtemp, nvec, tempvec, tempslipvec;
    double rrij0, rrij, temp, temp1, temp2;
    double a[5];
    int CSgroup = 950129;
    if ( UMB_order_param < 10)
    { /* 10: slip analysis for dislocation nucleation
         20: kink pair analysis for kink pair nucleation
         30: slip analysis for cross slip nucleation
         40: general slip analysis with general slip plane and vector
       */
       ERROR("caldislocationorder: UMB_order_param ("
             <<UMB_order_param<<") needs to be >= 10");
    }

    N_cluster = 0;
    N_lgst_cluster = 0;
    N_solid_P = 0;
    k0=0;
    Nc_i = 0;
    for (i=0;i< _NP;i++)
    {
        QLM_mark[i]=0;
        QLM_solid[i]=0;
        Size_cluster[i]=0;
	_QL[i]=0;
    }

    if (_SR2==NULL)
    {
 	FATAL("cal_dislocation_order: _SR2 == NULL");
	return;
    }
    N_bulk = _NP;
    for (i=0;i< _NP;i++)
    {
        _TOPOL[i] = _QL[i] = 0;
        if ( UMB_order_param/10 == 2 || UMB_order_param/10 == 4 )
        {  /* kink pair analysis */  /* general slip analysis */
 	    if (group[i] != UMB_order_group ) continue;
            if (HETERO_DN != 0 && YES_UMB == 1 && n_center < 20  )
            {
                ri00 = _H*_SR[i];
                if ( ri00.x < 0.0 || ri00.y > 0.0 ) continue;
            }
            if ( UMB_order_param == 43 && fabs(_SR2[i][1]) > 0.10 && YES_UMB == 1 ) continue;
        }
        if ( UMB_order_param == 30)
        {  /* cross slip analysis */
	   if (group[i] != CSgroup ) continue;
        }

       /*********************************************/
       if ( UMB_order_param == 10 )
       {/* slip analysis for dislocation nucleation */
        if (group[i] != UMB_order_group ) continue;
            if (HETERO_DN != 0 && YES_UMB == 1 && n_center < 20  )
            {
                ri00 = _H*_SR[i];
                if ( ri00.x < 0.0 || ri00.y > 0.0 ) continue;
            }
        nn_rc = 0;
        for (j=0;j<nn[i];j++)
        {
	    jpt = nindex[i][j];
	    if (jpt == i ) continue;

	    sij0 = _SR2[jpt]-_SR2[i]; sij0.subint();
            rij0 = _H*sij0; rrij0=rij0.norm();

            sij  = _SR[jpt] - _SR[i]; sij.subint();
	    rij = _H*sij; rrij=rij.norm();

            drtemp = rij-rij0;
	    temp = fabs(rrij0-rrij);

	    if (rrij > Rc_for_QLM) continue;
            nn_rc ++;
	    if (temp > _QL[i]) _QL[i]=temp;
	}
            if ( HETERO_DN != 0 ) {
                if (nn_rc < 10 )
	        {
                    _QL[i]=0;
                    N_bulk--;
            	}
            }
       } 
       /*********************************************/
       else if ( UMB_order_param == 20 )
       {/* kink pair analysis */
#if 0 /* two choices of coefficients */
           a[0]=1.0;
           a[1]=1.0;
           a[2]=-1.0;
           a[3]=1.0;
           a[4]=1.0;
#else
           a[0]=1.0/2.5;
           a[1]=3.0/2.5;
           a[2]=-2.0/2.5;
           a[3]=3.0/2.5;
           a[4]=1.0/2.5;
#endif
	    if ( ( i % 2) == 0) 
            {
                temp = a[0]*(_SR[i-4955][2]-_SR2[i-4955][2])
                      +a[1]*(_SR[i-4952][2]-_SR2[i-4952][2]) 
                      +a[2]*(_SR[i][2]-_SR2[i][2])
                      +a[3]*(_SR[i-6][2]   -_SR2[i-6][2])
                      +a[4]*(_SR[i+2][2]-_SR2[i+2][2]);
//                INFO_Printf("                %f                \n",temp);
            }
            else
            {
       	        temp = a[0]*(_SR[i-4946][2]-_SR2[i-4946][2])
                      +a[1]*(_SR[i-4952][2]-_SR2[i-4952][2])
                      +a[2]*(_SR[i][2]-_SR2[i][2])
                      +a[3]*(_SR[i+3][2]-_SR2[i+3][2])
                      +a[4]*(_SR[i+2][2]-_SR2[i+2][2]);
//                INFO_Printf("                %f                \n",temp);
            }
	    _QL[i]=temp;
       }
       else if ( UMB_order_param == 21 )
       {/* kink pair analysis for twice longer dislocation */
#if 0 /* two choices of coefficients */
           a[0]=1.0;
           a[1]=1.0;
           a[2]=-1.0;
           a[3]=1.0;
           a[4]=1.0;
#else
           a[0]=1.0/2.5;
           a[1]=3.0/2.5;
           a[2]=-2.0/2.5;
           a[3]=3.0/2.5;
           a[4]=1.0/2.5;
#endif
            if ( ( i % 2) == 0)
            {
                temp = a[0]*(_SR[i-9905][2]-_SR2[i-9905][2])
                      +a[1]*(_SR[i-9902][2]-_SR2[i-9902][2])
                      +a[2]*(_SR[i][2]-_SR2[i][2])
                      +a[3]*(_SR[i-6][2]   -_SR2[i-6][2])
                      +a[4]*(_SR[i+2][2]-_SR2[i+2][2]);
		if ( printUMBorder == 1 )
	            INFO_Printf("                %f                \n",temp);
            }
            else
            {
                temp = a[0]*(_SR[i-9896][2]-_SR2[i-9896][2])
                      +a[1]*(_SR[i-9902][2]-_SR2[i-9902][2])
                      +a[2]*(_SR[i][2]-_SR2[i][2])
                      +a[3]*(_SR[i+3][2]-_SR2[i+3][2])
                      +a[4]*(_SR[i+2][2]-_SR2[i+2][2]);
		if ( printUMBorder == 1 )
                    INFO_Printf("                %f                \n",temp);
            }
            _QL[i]=temp;
            temp2=temp/KinkDmax*10.0;
            _QLweighted[i]=temp2;
                if ( printUMBorder == 1 )
		    INFO_Printf("qlweighted = %f\n",_QLweighted[i]);
       }
       /*********************************************/
       else if ( UMB_order_param == 30 )
       {/* cross slip analysis */
            /* find the neighboring atoms on the adjacent plane */
            n = 0;
            nvec.set(0.0,1.0/3.0,-sqrt(8.0)/3.0);
            temp = 0;
            for(j=0;j<nn[i];j++)
            {
                jpt=nindex[i][j];
                if(jpt==i) continue;
                sij=_SR2[jpt]-_SR2[i]; sij.subint();
                rij=_H*sij; rrij=rij.norm();

                /* select nearest neighbor atoms */
                if (rrij > Rc_for_QLM) continue;

                /* select atoms on the adjacent plane */
                temp1 = dot(rij,nvec)/(latticeconst[0]*sqrt(3.0)/4.0);
                if (temp1 < 0.5 || temp1 > 1.5) continue;

                /* compute the averaged displacement of these neighbor atoms */
                temp += _SR[jpt][0]-_SR2[jpt][0];
                n++;
   
                /* graphically verify we got the correct atoms */
                //fixed[jpt] = 1;
            }
            _QL[i] = (_SR[i][0]-_SR2[i][0] - temp/n)*_H[0][0];
//            INFO_Printf("_QL[%d] = %g  n = %d\n",i,_QL[i],n);
       }
       /*********************************************/
       else if ( UMB_order_param/10 == 4 )
       {/* cross slip analysis */
            /* find the neighboring atoms on the adjacent plane */
	   if (UMB_nn_ref[i] == UMB_NNBR ) {  // glide set in Diamond
	      //printf("UMB_nn_ref[%d] = %d, UMB_NNBR = %d\n", i, UMB_nn_ref[i], UMB_NNBR);
              tempvec.clear();
              tempslipvec.clear();

              for(j=0;j<UMB_nn_ref[i];j++)
                {
                 jpt = UMB_nindex_ref[i][j];
                 /* compute the averaged displacement of these neighbor atoms */
                 tempvec += _SR[jpt]-_SR2[jpt];
                }

	        tempslipvec = _SR[i]-_SR2[i]; tempvec=tempvec/UMB_nn_ref[i];
	        tempslipvec = _H*(tempslipvec - tempvec );
                if ( UMB_order_param != 43 ) {
                    _QL[i] = dot(UMB_slipvec,tempslipvec);
		}
                else 
                { 
		    nvec.x=-1.0*UMB_slipvec.x; nvec.y=UMB_slipvec.y; nvec.z=UMB_slipvec.z;
                    if ( _SR2[i][1] >= 0.0 )
			_QL[i] = dot(UMB_slipvec,tempslipvec);
                    else
			_QL[i] = dot(nvec,tempslipvec);
                }

            } else {
            _QL[i]=0;
            }
/* 
           if ( group[i] == UMB_order_group && _QL[i] > QLM_cutoff )
           {
            INFO_Printf("\nindex= %d, group = %d\n",i,group[i]);
            INFO_Printf("slipvec.x = %f slipvec.y = %f slipvec.z = %f\n",tempslipvec.x,tempslipvec.y,tempslipvec.z);
//            INFO_Printf("nvec.x = %f nvec.y = %f nvec.z = %f\n",UMB_nvec.x,UMB_nvec.y,UMB_nvec.z);
            INFO_Printf("_QL[%d] = %g  n = %d UMB_thickness =%f\n\n",i,_QL[i],n,UMB_thickness);
           }
*/
       }
       /*********************************************/
       else {
          ERROR("unknown UMB_order_param ("<<UMB_order_param<<")!");
       }

        _TOPOL[i]=_QL[i];
	
        if ( UMB_noslipdirection == 0 ) 
        {
            if ( _QL[i] > QLM_cutoff )
            {
                    N_solid_P++;
                    QLM_solid[i]=1;
            }
        } 
        else
        {
            if ( fabs(_QL[i]) > QLM_cutoff )
            {
                    N_solid_P++;
                    QLM_solid[i]=1;
            }
	}
    }
      
    for (i=0;i<_NP;i++)
    {      
        if ( QLM_mark[i]==0 && QLM_solid[i]==1 )
        {
	    N_cluster++;
	    N_cluster_temp = 0;
            N_lgst_inDOUBLE = 0;
	    if (N_cluster == 1)
	 	Nc_i = i;
            if ( UMB_order_group == 0 )
	        DFS(i,0,0);
            else
	        DFS2(i,UMB_order_group,0,0);
            if ( UMB_order_param/10 == 2 )
		N_cluster_temp=int(N_lgst_inDOUBLE);	    

	    Size_cluster[N_cluster] = N_cluster_temp;
	    if (N_cluster_temp > N_lgst_cluster)
	    {
		N_lgst_cluster = N_cluster_temp;
		Nc_i=i;
	    }
	    if (YES_UMB == 1 && UMB_K_in_Kelvin < 0 && N_cluster_temp > 0 && curstep > UMB_equilstep)
	    {
	        for (j=0;j<=2*delta_n;j++)
	            if (Narray[j]==N_cluster_temp) {
		        Occurence[j]++;
		        Probability[j]++;
	            }
	    }
	}
	else 
	    k0++;
    }
    if (YES_UMB == 1 && UMB_K_in_Kelvin < 0 && curstep > UMB_equilstep )
	for (j=0;j<=2*delta_n;j++) 
		if (Narray[j]==0) {
		    Occurence[j]=Occurence[j]+k0;
		    Probability[j]=Probability[j]+k0;
		}

    if (N_lgst_index !=0)
    {
	for(i=0;i<_NP;i++)
	    QLM_mark[i]=0;
	if (UMB_order_group == 0 )
	    DFS(Nc_i,1,N_lgst_index);
	else
	    DFS2(Nc_i,UMB_order_group,1,N_lgst_index);
    }
}

void MDFrame::DFS(int i, int On_Off, double index )
{
    int j, jpt;
    double rrij;
    Vector3 sij, rij;

    if (QLM_mark[i]==1 || QLM_solid[i]==0)
	return;
    else
    {
	QLM_mark[i]=1;
	if ( On_Off > 0 ) 
	{
	    _QL[i]=index;
	    _TOPOL[i]=_QL[i];
            if ( YES_SEP > 0 )
            {
               sij=_SR[i]-_SRforCM; sij.subint();
               _CLUSTER_CM+=sij;
            }
        }
	N_cluster_temp++;
	for (j=0;j<nn[i];j++)
	{
	    jpt = nindex[i][j];
	    if ( jpt == i || QLM_solid[jpt] == 0 || QLM_mark[jpt] == 1)
		continue;
	    sij=_SR[jpt]-_SR[i]; sij.subint();
	    rij=_H*sij; rrij=rij.norm();
	   
	    if (rrij > Rc_for_QLM) continue;
	    DFS(jpt, On_Off, index);
        }
    }
}

void MDFrame::DFS2(int i, int groupindex, int On_Off, double index )
{
    int j, jpt; //, iD, rD;
    double rrij;
    Vector3 sij, rij;

    if (QLM_mark[i]==1 || QLM_solid[i]==0)
	return;
    else
    {
        if ( group[i]!=groupindex ) return;

        QLM_mark[i]=1;
        if ( On_Off > 0 )
        {
	    _QL[i]=index;
	    _TOPOL[i]=_QL[i];
	}
	N_cluster_temp++;
        N_lgst_inDOUBLE=N_lgst_inDOUBLE+_QLweighted[i];
	for (j=0;j<nn[i];j++)
	{
	    jpt = nindex[i][j];
	    if ( jpt == i || QLM_solid[jpt] == 0 || QLM_mark[jpt] == 1)
		continue;
	    sij=_SR[jpt]-_SR[i]; sij.subint();
	    rij=_H*sij; rrij=rij.norm();
	   
	    if (rrij > Rc_for_QLM) continue;
	    DFS2(jpt, groupindex, On_Off, index);
        }
    }
}

void MDFrame::DFS1(int i, int On_Off, double index, double index1 )
{
    int j, jpt;
    double rrij;
    Vector3 sij, rij;

    if (QLM_mark[i]==1 && QLM_solid[i]==0)
	return;
    else if (QLM_mark[i]==0 && QLM_solid[i]==0)
    {
        QLM_mark[i]=1;
	if (On_Off > 0 )
	{
	    _QL[i]=index1;
	    _TOPOL[i]=_QL[i];
	    N_lgst_skin++;
            if ( YES_SEP > 0 )
            {
               sij=_SR[i]-_SRforCM; sij.subint();
               _CLUSTER_CM+=sij;
            }
	}
        N_skin_P++;
        N_cluster_temp++;
        return; 
    }
    else
    {
	QLM_mark[i]=1;
	if ( On_Off > 0 )
	{
	    _QL[i]=index;
	    _TOPOL[i]=_QL[i];
	}
	N_cluster_temp++;
	for (j=0;j<nn[i];j++)
	{
	    jpt = nindex[i][j];
	    if ( jpt == i || QLM_mark[jpt] == 1)
		continue;
	    sij=_SR[jpt]-_SR[i]; sij.subint();
	    rij=_H*sij; rrij=rij.norm();
	   
	    if (rrij > Rc_for_QLM) continue;
	    DFS1(jpt, On_Off, index, index1);
        }
    }
}

void MDFrame::allocQLM()
{
    int size;

    if (L_for_QLM < 0)
    {
        ERROR("L_for_QLM = "<<L_for_QLM<<" skip allocQLM");
        return;
    }
    if (_NP <= 0)
    {
        ERROR("NP = "<<_NP<<" skip allocQLM");
        return;
    }

    size = (2*L_for_QLM + 2)*_NP;

    Realloc(_QLM,double,size);
   
    Realloc(QLM_solid,int,_NP);

    Realloc(QLM_mark,int,_NP);

    Realloc(_QL,double,_NP);
 
    Realloc(_QLweighted,double,_NP);

}

