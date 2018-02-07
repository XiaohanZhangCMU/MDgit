/* OpenMP implementation */

double PhaseFieldFrame::FreeEnergy_single_omp() /* traditional (single) phase field model */
{
    int n, i, j, k, ip, im, jp, jm, kp, km, ind;
    double h2, f, g, phi, dphidx, dphidy, dphidz, d2phi, dphi_SQR;
    double absdphi, EPS_loc, EPS_prefactor, triplesum;
    double nx_lab, ny_lab, nz_lab, nx_cryst, ny_cryst, nz_cryst;
    double dnx_lab_dphidx, dny_lab_dphidx, dnz_lab_dphidx; 
    double dnx_lab_dphidy, dny_lab_dphidy, dnz_lab_dphidy; 
    double dnx_lab_dphidz, dny_lab_dphidz, dnz_lab_dphidz; 
    double dEPS_dnx_lab, dEPS_dny_lab, dEPS_dnz_lab;
    double dEPS_dnx_cryst, dEPS_dny_cryst, dEPS_dnz_cryst;
    double dEPS_dphidx, dEPS_dphidy, dEPS_dphidz;
    double avg_dPHIdt0;
    
#ifdef _OPENMP
    /* only local arrays can be shared by openmp threads */
    double *shm_PHI[MAXNUMFIELDS], *shm_d2PHI[MAXNUMFIELDS]; 
    double *shm_dPHIdx[MAXNUMFIELDS], *shm_dPHIdy[MAXNUMFIELDS], *shm_dPHIdz[MAXNUMFIELDS];
    double *shm_dFdPHI[MAXNUMFIELDS], *shm_EPS_loc[MAXNUMFIELDS];
    double *shm_tmp_x[MAXNUMFIELDS], *shm_tmp_y[MAXNUMFIELDS], *shm_tmp_z[MAXNUMFIELDS];
    double *shm_d2dFdPHI[MAXNUMFIELDS];
    double *shm_dPHIdt0[MAXNUMFIELDS], *shm_dPHIdt[MAXNUMFIELDS];
    for(n=0;n<num_fields;n++)
    {
        shm_PHI[n]   = _PHI[n]; shm_d2PHI[n] = _d2PHI[n];
        shm_dPHIdx[n] = _dPHIdx[n]; shm_dPHIdy[n] = _dPHIdy[n]; shm_dPHIdz[n] = _dPHIdz[n]; 
        shm_dFdPHI[n] = _dFdPHI[n]; shm_EPS_loc[n] = _EPS_loc[n];
        shm_tmp_x[n] = tmp_x[n]; shm_tmp_y[n] = tmp_y[n]; shm_tmp_z[n] = tmp_z[n];
        shm_d2dFdPHI[n] = _d2dFdPHI[n];
        shm_dPHIdt0[n] = _dPHIdt0[n];
        shm_dPHIdt[n]  = _dPHIdt[n];
    }
#else
    
    /* make the code compilable when not using -openmp */
#define shm_PHI     _PHI
#define shm_d2PHI   _d2PHI
#define shm_dPHIdx  _dPHIdx
#define shm_dPHIdy  _dPHIdy
#define shm_dPHIdz  _dPHIdz
#define shm_dFdPHI  _dFdPHI
#define shm_EPS_loc _EPS_loc
#define shm_tmp_x    tmp_x
#define shm_tmp_y    tmp_y
#define shm_tmp_z    tmp_z
#define shm_d2dFdPHI _d2dFdPHI
#define shm_dPHIdt   _dPHIdt
#define shm_dPHIdt0  _dPHIdt0
    
#endif
    
    _F = 0;
    _Fraw = 0;
    
    if (gridsize<=0)
        FATAL("FreeEnergy: gridesize = "<<gridsize);
    
    /* compute spatial gradients */
    h2 = gridsize*gridsize;
    for(n=0;n<num_fields;n++)
    {
        double loc_F, sum_F;
        sum_F = 0;
#pragma omp parallel for reduction(+:sum_F) \
private(n,i,j,k,ip,im,jp,jm,kp,km,ind,loc_F,phi,d2phi,dphidx,dphidy,dphidz,dphi_SQR) \
private(f,g,EPS_loc,EPS_prefactor) \
private(nx_lab,ny_lab,nz_lab,nx_cryst,ny_cryst,nz_cryst,triplesum,absdphi) \
private(dnx_lab_dphidx,dny_lab_dphidx,dnz_lab_dphidx) \
private(dnx_lab_dphidy,dny_lab_dphidy,dnz_lab_dphidy) \
private(dnx_lab_dphidz,dny_lab_dphidz,dnz_lab_dphidz) \
private(dEPS_dnx_cryst,dEPS_dny_cryst,dEPS_dnz_cryst) \
private(dEPS_dnx_lab,  dEPS_dny_lab,  dEPS_dnz_lab  ) \
private(dEPS_dphidx,   dEPS_dphidy,   dEPS_dphidz   ) \
shared(shm_PHI,shm_d2PHI,shm_dPHIdx,shm_dPHIdy,shm_dPHIdz) \
shared(shm_dFdPHI,shm_EPS_loc,shm_tmp_x,shm_tmp_y,shm_tmp_z)
        
        for(i=0;i<_NX;i++)
        {
            loc_F = 0;
            ip = (i+1)%_NX;
            im = (i-1+_NX)%_NX;
            for (j=0;j<_NY;j++)
            {
                jp = (j+1)%_NY;
                jm = (j-1+_NY)%_NY;
                for (k=0;k<_NZ;k++)
                {
                    kp = (k+1)%_NZ;
                    km = (k-1+_NZ)%_NZ;
                    
                    ind = i*_NY*_NZ+j*_NZ+k;
                    phi = _PHI[n][ind];
                    d2phi = (shm_PHI[n][ip*_NY*_NZ+j*_NZ+k] + shm_PHI[n][im*_NY*_NZ+j*_NZ+k]
                             +shm_PHI[n][i*_NY*_NZ+jp*_NZ+k] + shm_PHI[n][i*_NY*_NZ+jm*_NZ+k]
                             +shm_PHI[n][i*_NY*_NZ+j*_NZ+kp] + shm_PHI[n][i*_NY*_NZ+j*_NZ+km]
                             -shm_PHI[n][i*_NY*_NZ+j*_NZ+k] * 6.0 ) / h2;
                    
                    shm_d2PHI[n][ind] = d2phi;
                    
                    dphidx = (shm_PHI[n][ip*_NY*_NZ+j*_NZ+k] - shm_PHI[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize);
                    dphidy = (shm_PHI[n][i*_NY*_NZ+jp*_NZ+k] - shm_PHI[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize);
                    dphidz = (shm_PHI[n][i*_NY*_NZ+j*_NZ+kp] - shm_PHI[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize);
                    
                    shm_dPHIdx[n][ind] = dphidx;
                    shm_dPHIdy[n][ind] = dphidy;
                    shm_dPHIdz[n][ind] = dphidz;
                    
                    f = phi*phi*(1-phi)*(1-phi);
                    g = 2.0*phi*(1-phi)*(1-2.0*phi);
                    
                    //_F += f * _U[n][n]; 
                    loc_F += f * _U[n][n]; 
                    shm_dFdPHI[n][ind] = g * _U[n][n];
                    
                    if ((_EPS1[n][n] == 0)&&(_EPS2[n][n] == 0)&&(_EPS3[n][n] == 0))
                    { /* isotropic interface energy */
                        EPS_loc = _EPS[n][n];
                        EPS_prefactor = 0.5*SQR(EPS_loc);
                        
                        loc_F += (dphidx*dphidx + dphidy*dphidy + dphidz*dphidz) * EPS_prefactor;
                        shm_dFdPHI[n][ind] += - 2.0 * d2phi * EPS_prefactor;
                        
                    }
                    else
                    { /* cubic anisotropic interface energy */
                        dphi_SQR   = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);
                        absdphi = sqrt(dphi_SQR + 1e-24);
                        nx_lab = dphidx / absdphi;
                        ny_lab = dphidy / absdphi;
                        nz_lab = dphidz / absdphi;
                        
#define _R _ROT_MATRIX 
                        nx_cryst = _R[0][0]*nx_lab + _R[0][1]*ny_lab + _R[0][2]*nz_lab;
                        ny_cryst = _R[1][0]*nx_lab + _R[1][1]*ny_lab + _R[1][2]*nz_lab;
                        nz_cryst = _R[2][0]*nx_lab + _R[2][1]*ny_lab + _R[2][2]*nz_lab;
                        
                        dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
                        dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
                        dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);
                        
                        dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
                        dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
                        dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);
                        
                        dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
                        dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
                        dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);
                        
                        triplesum = SQR(nx_cryst)*SQR(ny_cryst) + SQR(ny_cryst)*SQR(nz_cryst) + SQR(nz_cryst)*SQR(nx_cryst);
                        EPS_loc   = _EPS[n][n] + _EPS1[n][n]*triplesum + _EPS2[n][n]*SQR(nx_cryst*ny_cryst*nz_cryst)
                        + _EPS3[n][n]*SQR(triplesum);
                        shm_EPS_loc[n][ind]   = EPS_loc;
                        
                        loc_F += 0.5 * SQR(EPS_loc) * dphi_SQR;
                        /* dFdPHI += -d/dx ( SQR(EPS_loc) * dphidx + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidx) )
                         -d/dy ( SQR(EPS_loc) * dphidy + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidy) )
                         -d/dz ( SQR(EPS_loc) * dphidz + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidz) ) ;
                         */
                        
                        dEPS_dnx_cryst = 2.*_EPS1[n][n]*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))
                        + 2.*_EPS2[n][n]*nx_cryst*SQR(ny_cryst)*SQR(nz_cryst)
                        + 4.*_EPS3[n][n]*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))*triplesum;
                        dEPS_dny_cryst = 2.*_EPS1[n][n]*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))
                        + 2.*_EPS2[n][n]*ny_cryst*SQR(nz_cryst)*SQR(nx_cryst)
                        + 4.*_EPS3[n][n]*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))*triplesum;
                        dEPS_dnz_cryst = 2.*_EPS1[n][n]*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))
                        + 2.*_EPS2[n][n]*nz_cryst*SQR(nx_cryst)*SQR(ny_cryst)
                        + 4.*_EPS3[n][n]*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))*triplesum;
                        
                        dEPS_dnx_lab = _R[0][0]*dEPS_dnx_cryst + _R[1][0]*dEPS_dny_cryst + _R[2][0]*dEPS_dnz_cryst;
                        dEPS_dny_lab = _R[0][1]*dEPS_dnx_cryst + _R[1][1]*dEPS_dny_cryst + _R[2][1]*dEPS_dnz_cryst;
                        dEPS_dnz_lab = _R[0][2]*dEPS_dnx_cryst + _R[1][2]*dEPS_dny_cryst + _R[2][2]*dEPS_dnz_cryst;
                        
                        dEPS_dphidx = dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx;
                        dEPS_dphidy = dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy;
                        dEPS_dphidz = dEPS_dnx_lab*dnx_lab_dphidz + dEPS_dny_lab*dny_lab_dphidz + dEPS_dnz_lab*dnz_lab_dphidz;
                        
                        //tmp_x[n][ind] = SQR(EPS_loc)*dphidx + dphi_SQR*EPS_loc*dEPS_dphidx;
                        //tmp_y[n][ind] = SQR(EPS_loc)*dphidy + dphi_SQR*EPS_loc*dEPS_dphidy;
                        //tmp_z[n][ind] = SQR(EPS_loc)*dphidz + dphi_SQR*EPS_loc*dEPS_dphidz;
                        
                        /* new scheme for computing variational derivative, 2012/04/27  */
                        shm_tmp_x[n][ind] = dphi_SQR*EPS_loc*dEPS_dphidx;
                        shm_tmp_y[n][ind] = dphi_SQR*EPS_loc*dEPS_dphidy;
                        shm_tmp_z[n][ind] = dphi_SQR*EPS_loc*dEPS_dphidz;
                        
                        shm_dFdPHI[n][ind] += - 1.0 * SQR(EPS_loc) * d2phi;
                        
                    }
                }
            }
            sum_F = sum_F + loc_F;
            
        }
        _F = sum_F * CUBE(gridsize);
        _Fraw = _F;
        
        if ((_EPS1[n][n] != 0)||(_EPS2[n][n] != 0)||(_EPS3[n][n] != 0))
        { /* cubic anisotropic interface energy */
#pragma omp parallel for default(none) \
private(n,i,j,k,ip,im,jp,jm,kp,km,ind) \
shared(shm_dFdPHI,shm_EPS_loc,shm_tmp_x,shm_tmp_y,shm_tmp_z) \
shared(shm_dPHIdx,shm_dPHIdy,shm_dPHIdz)
            for(i=0;i<_NX;i++)
            {
                ip = (i+1)%_NX;
                im = (i-1+_NX)%_NX;
                for (j=0;j<_NY;j++)
                {
                    jp = (j+1)%_NY;
                    jm = (j-1+_NY)%_NY;
                    for (k=0;k<_NZ;k++)
                    {
                        kp = (k+1)%_NZ;
                        km = (k-1+_NZ)%_NZ;
                        
                        ind = i*_NY*_NZ+j*_NZ+k;
                        shm_dFdPHI[n][ind] +=
                        - (shm_tmp_x[n][ip*_NY*_NZ+j*_NZ+k]-shm_tmp_x[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize)
                        - (shm_tmp_y[n][i*_NY*_NZ+jp*_NZ+k]-shm_tmp_y[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize)
                        - (shm_tmp_z[n][i*_NY*_NZ+j*_NZ+kp]-shm_tmp_z[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize) ;
                        
                        /* new scheme for computing variational derivative, 2012/04/27  */
                        shm_dFdPHI[n][ind] +=
                        - (SQR(shm_EPS_loc[n][ip*_NY*_NZ+j*_NZ+k])-SQR(shm_EPS_loc[n][im*_NY*_NZ+j*_NZ+k])) / (2.0*gridsize) * shm_dPHIdx[n][ind]
                        - (SQR(shm_EPS_loc[n][i*_NY*_NZ+jp*_NZ+k])-SQR(shm_EPS_loc[n][i*_NY*_NZ+jm*_NZ+k])) / (2.0*gridsize) * shm_dPHIdy[n][ind]
                        - (SQR(shm_EPS_loc[n][i*_NY*_NZ+j*_NZ+kp])-SQR(shm_EPS_loc[n][i*_NY*_NZ+j*_NZ+km])) / (2.0*gridsize) * shm_dPHIdz[n][ind];
                    }
                }
            }
        }
        
        
        if (dynamics_type == 1) /* Cahn-Hillard equation */
        {
#pragma omp parallel for default(none) \
private(n,i,j,k,ip,im,jp,jm,kp,km,ind,h2) \
shared(shm_d2dFdPHI,shm_dFdPHI)
            for(i=0;i<_NX;i++)
            {
                ip = (i+1)%_NX;
                im = (i-1+_NX)%_NX;
                for (j=0;j<_NY;j++)
                {
                    jp = (j+1)%_NY;
                    jm = (j-1+_NY)%_NY;
                    
                    for (k=0;k<_NZ;k++)
                    {
                        kp = (k+1)%_NZ;
                        km = (k-1+_NZ)%_NZ;
                        
                        ind = i*_NY*_NZ+j*_NZ+k;
                        shm_d2dFdPHI[n][ind] = 
                        (shm_dFdPHI[n][ip*_NY*_NZ+j*_NZ+k] + shm_dFdPHI[n][im*_NY*_NZ+j*_NZ+k]
                         +shm_dFdPHI[n][i*_NY*_NZ+jp*_NZ+k] + shm_dFdPHI[n][i*_NY*_NZ+jm*_NZ+k]
                         +shm_dFdPHI[n][i*_NY*_NZ+j*_NZ+kp] + shm_dFdPHI[n][i*_NY*_NZ+j*_NZ+km]
                         -shm_dFdPHI[n][ind] * 6.0 ) / h2;
                    }
                }
            }
        }
        
#pragma omp parallel for default(none) \
private(n,ind)  shared(shm_d2dFdPHI,shm_dFdPHI,shm_dPHIdt0)
        for(ind=0;ind<_NX*_NY*_NZ;ind++)
        {
            if ((dynamics_type == 0) || (dynamics_type == 2)) /* Cahn-Hillard equation */
            {
                shm_dPHIdt0[n][ind] = -1.0 * Mob_GL * shm_dFdPHI[n][ind];
            }
            else if (dynamics_type == 1) /* Ginzburg-Landau equation */
            {
                shm_dPHIdt0[n][ind] = -1.0 * Mob_D * shm_d2dFdPHI[n][ind];
            }
        }
        
        if ((dynamics_type == 0) || (dynamics_type == 1)) /* no additional constraints */
        {
            /* 0: Ginzburg-Landau equation: not conserved */
            /* 1: Cahn-Hilliard equation: conserved */
#pragma omp parallel for default(none) \
private(n,ind)  shared(shm_dPHIdt,shm_dPHIdt0)
            for(ind=0;ind<_NX*_NY*_NZ;ind++)
                shm_dPHIdt[n][ind] = shm_dPHIdt0[n][ind];
        }
        else if (dynamics_type == 2) 
        {
            /* Ginzburg-Landau equation but with constant volume constraint */ 
            
            avg_dPHIdt0 = 0;
#pragma omp parallel for reduction(+:avg_dPHIdt0) \
private(n,ind)  shared(shm_dPHIdt0)
            for(ind=0;ind<_NX*_NY*_NZ;ind++)
                avg_dPHIdt0 += shm_dPHIdt0[n][ind];
            
            avg_dPHIdt0 /= (_NX*_NY*_NZ);
            
#pragma omp parallel for default(none) \
private(n,ind) \
shared(avg_dPHIdt0,shm_dPHIdt,shm_dPHIdt0)
            for(ind=0;ind<_NX*_NY*_NZ;ind++)
                _dPHIdt[n][ind] = _dPHIdt0[n][ind] - avg_dPHIdt0;
        }
        else
        {
            ERROR("unknown dynamics_type = "<<dynamics_type);
        }
        
    }/* end of for(n=0;n<num_fields;n++) */
    
    
    /* Need to make this work for OpenMP */
    /* calculate the max value of dphidt */ 
    double Gmiddle = 0;
    
#pragma omp parallel for default(none) \
private(n,ind) \
shared(shm_dPHIdt,Gmiddle)  
    
    for(n=0;n<num_fields;n++)
    {
        for(ind=0;ind<_NX*_NY*_NZ;ind++)
        {
            if ( Gmiddle <= fabs(_dPHIdt[n][ind]) )  Gmiddle = fabs(_dPHIdt[n][ind]);
        }
    }
    _G = Gmiddle;
    
    return _F;
}




double PhaseFieldFrame::FreeEnergy_multi_omp() /* multi phase field model */
{
    int n, p, i, j, k, ip, ip2, ip3, im, im2, im3, jp, jp2, jp3, jm, jm2, jm3, kp, kp2, kp3, km, km2, km3, ind;
    double h2, f, g, phi, dphidx, dphidy, dphidz, d2phi, dphi_SQR, fp;
    double absdphi, EPS_loc, EPS_prefactor, triplesum;
    double nx_lab, ny_lab, nz_lab, nx_cryst, ny_cryst, nz_cryst;
    double dnx_lab_dphidx, dny_lab_dphidx, dnz_lab_dphidx; 
    double dnx_lab_dphidy, dny_lab_dphidy, dnz_lab_dphidy; 
    double dnx_lab_dphidz, dny_lab_dphidz, dnz_lab_dphidz; 
    double dEPS_dnx_lab, dEPS_dny_lab, dEPS_dnz_lab;
    double dEPS_dnx_cryst, dEPS_dny_cryst, dEPS_dnz_cryst;
    double dEPS_dphidx, dEPS_dphidy, dEPS_dphidz;
    double eps1, eps2, eps3, dEPS_dphidx_loc, dEPS_dphidy_loc, dEPS_dphidz_loc;
    double grad_loc_x, grad_loc_y, grad_loc_z, grad_loc_SQR;
    double grad_term, avg_dPHIdt0;
    
#ifdef _OPENMP
    /* only local arrays can be shared by openmp threads */
    double *shm_PHI[MAXNUMFIELDS], *shm_d2PHI[MAXNUMFIELDS]; 
    double *shm_dPHIdx[MAXNUMFIELDS], *shm_dPHIdy[MAXNUMFIELDS], *shm_dPHIdz[MAXNUMFIELDS];
    double *shm_dFdPHIi[MAXNUMFIELDS], *shm_EPS_loc[MAXNUMFIELDS];
    double *shm_tmp_x[MAXNUMFIELDS], *shm_tmp_y[MAXNUMFIELDS], *shm_tmp_z[MAXNUMFIELDS];
    double *shm_d2dFdPHI[MAXNUMFIELDS];
    double *shm_dPHIdt0[MAXNUMFIELDS], *shm_dPHIdt[MAXNUMFIELDS], *shm_dPHIdtCH[MAXNUMFIELDS], *shm_d2dPHIdt0[MAXNUMFIELDS];
    double *shm_K_S, *shm_K_L, *shm_K_LS, *shm_K_LV, *shm_K_SV, *shm_K_other, *shm_K_D;
    double *shm_matrix_x, *shm_matrix_y;
    double *shm_NW_orig;
    
    for(n=0;n<num_fields;n++)
    {
        shm_PHI[n]   = _PHI[n]; shm_d2PHI[n] = _d2PHI[n];
        shm_dPHIdx[n] = _dPHIdx[n]; shm_dPHIdy[n] = _dPHIdy[n]; shm_dPHIdz[n] = _dPHIdz[n]; 
        shm_dFdPHIi[n] = _dFdPHIi[n]; shm_EPS_loc[n] = _EPS_loc[n];
        shm_tmp_x[n] = tmp_x[n]; shm_tmp_y[n] = tmp_y[n]; shm_tmp_z[n] = tmp_z[n];
        shm_d2dFdPHI[n] = _d2dFdPHI[n];
        shm_dPHIdt0[n] = _dPHIdt0[n];
        shm_dPHIdt[n]  = _dPHIdt[n];
        shm_dPHIdtCH[n]  = _dPHIdtCH[n];
        shm_d2dPHIdt0[n] = _d2dPHIdt0[n];
    }
    
    shm_K_S = _K_S;
    shm_K_L = _K_L;
    shm_K_LS = _K_LS;
    shm_K_LV = _K_LV;
    shm_K_SV = _K_SV;
    shm_K_other = _K_other;
    shm_K_D = _K_D;
    shm_matrix_x = matrix_x;
    shm_matrix_y = matrix_y;
    shm_NW_orig = _NW_orig;

    
#else
    
    /* make the code compilable when not using -openmp */
#define shm_PHI     _PHI
#define shm_d2PHI   _d2PHI
#define shm_dPHIdx  _dPHIdx
#define shm_dPHIdy  _dPHIdy
#define shm_dPHIdz  _dPHIdz
#define shm_dFdPHIi  _dFdPHIi
#define shm_EPS_loc _EPS_loc
#define shm_tmp_x    tmp_x
#define shm_tmp_y    tmp_y
#define shm_tmp_z    tmp_z
#define shm_d2dFdPHI _d2dFdPHI
#define shm_dPHIdt   _dPHIdt
#define shm_dPHIdt0  _dPHIdt0
#define shm_d2dPHIdt0 _d2dPHIdt0
#define shm_dPHIdtCH _dPHIdtCH
#define shm_K_L _K_S
#define shm_K_LS _K_LS
#define shm_K_LV _K_LV
#define shm_K_SV _K_SV
#define shm_K_other _K_other
#define shm_K_D _K_D
#define shm_matrix_x matrix_x
#define shm_matrix_y matrix_y
#define shm_NW_orig _NW_orig    
#endif
    
    _F = 0;
    _Fraw = 0;
    
    if (gridsize<=0)
        FATAL("FreeEnergy: gridesize = "<<gridsize);
    
    /* Zero out time derivative of phase fields and assgin kinetic coefficient for C-H*/
    
    
#pragma omp parallel for default(none) \
private(n,ind) \
shared(shm_dPHIdt0,shm_dFdPHIi,shm_dPHIdtCH,shm_K_D) 
    
    for (ind=0; ind<_NX*_NY*_NZ; ind++)
    {
        shm_dPHIdt0[0][ind] = shm_dPHIdt0[1][ind] = shm_dPHIdt0[2][ind] = 0;
        shm_dFdPHIi[0][ind] = shm_dFdPHIi[1][ind] = shm_dFdPHIi[2][ind] = 0;
        shm_dPHIdtCH[0][ind] = shm_dPHIdtCH[1][ind] = shm_dPHIdtCH[2][ind] = 0;
        shm_K_D[ind] = Mob_D;
    }
    
    /* Reset liquid chemical potential, will be updated by constraint */
    if (dynamics_type == 2 || dynamics_type == 3 || dynamics_type == 8)
    {
        _M[0][1] = 0; _M[0][2] = _M[1][2];
        _M[1][0] = 0; _M[2][0] = _M[2][1];
        _MU[0] = _MU[1];
    }
    
    if (dynamics_type == 4 || dynamics_type == 5 || dynamics_type == 6 || dynamics_type == 7)
    {
        _MU[0] = _MU[1];
        _MU[2] = _MU[1];
    }
    
    
    /* compute spatial gradients */
    h2 = gridsize*gridsize;
    
#pragma omp parallel for default(none) \
private(n,i,j,k,ip,im,jp,jm,kp,km,ind,phi,d2phi,dphidx,dphidy,dphidz) \
shared(h2,shm_PHI,shm_d2PHI,shm_dPHIdx,shm_dPHIdy,shm_dPHIdz) 
    
    for(i=0;i<_NX;i++)
    {
        ip = (i+1)%_NX;
        im = (i-1+_NX)%_NX;
        for (j=0;j<_NY;j++)
        {
            jp = (j+1)%_NY;
            jm = (j-1+_NY)%_NY;
            for (k=0;k<_NZ;k++)
            {
                kp = (k+1)%_NZ;
                km = (k-1+_NZ)%_NZ;
                
                ind = i*_NY*_NZ+j*_NZ+k;
                for(n=0;n<num_fields;n++)
                { 
                    phi = shm_PHI[n][ind];
                    d2phi = ( shm_PHI[n][ip*_NY*_NZ+j*_NZ+k] + shm_PHI[n][im*_NY*_NZ+j*_NZ+k]
                             +shm_PHI[n][i*_NY*_NZ+jp*_NZ+k] + shm_PHI[n][i*_NY*_NZ+jm*_NZ+k]
                             +shm_PHI[n][i*_NY*_NZ+j*_NZ+kp] + shm_PHI[n][i*_NY*_NZ+j*_NZ+km]
                             -shm_PHI[n][ind] * 6.0 ) / h2;
                    
                    shm_d2PHI[n][ind] = d2phi;
                    
                    dphidx = (shm_PHI[n][ip*_NY*_NZ+j*_NZ+k] - shm_PHI[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize);
                    dphidy = (shm_PHI[n][i*_NY*_NZ+jp*_NZ+k] - shm_PHI[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize);
                    dphidz = (shm_PHI[n][i*_NY*_NZ+j*_NZ+kp] - shm_PHI[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize);
                    
                    shm_dPHIdx[n][ind] = dphidx;
                    shm_dPHIdy[n][ind] = dphidy;
                    shm_dPHIdz[n][ind] = dphidz;
                    
                }
            }
        }
    }
    
    
    
    /* Compute Free energy and time derivative of phase fields */
#pragma omp parallel for default(none) \
private(i,j,k,ip,ip2,ip3,im,im2,im3,jp,jp2,jp3,jm,jm2,jm3,kp,kp2,kp3,km,km2,km3,ind) \
shared(shm_PHI,shm_K_S,shm_K_L,shm_K_LS,shm_K_LV,shm_K_SV,shm_K_other)
    
    for (i=0;i<_NX;i++)
    {
        ip  = (i+1)%_NX;
        ip2 = (i+2)%_NX;
        ip3 = (i+3)%_NX;
        im  = (i-1+_NX)%_NX;
        im2 = (i-2+_NX)%_NX;
        im3 = (i-3+_NX)%_NX;
        for (j=0;j<_NY;j++)
        {
            jp  = (j+1)%_NY;
            jp2 = (j+2)%_NY;
            jp3 = (j+3)%_NY;
            jm  = (j-1+_NY)%_NY;
            jm2 = (j-2+_NY)%_NY;
            jm3 = (j-3+_NY)%_NY;
            for (k=0;k<_NZ;k++)
            {
                kp  = (k+1)%_NZ;
                kp2 = (k+2)%_NZ;
                kp3 = (k+3)%_NZ;
                km  = (k-1+_NZ)%_NZ;
                km2 = (k-2+_NZ)%_NZ; 
                km3 = (k-3+_NZ)%_NZ;
                ind = i*_NY*_NZ+j*_NZ+k;
                
                /* define the solid phase boundary */
                if (     (((shm_PHI[1][ind]-0.5)*(shm_PHI[1][ip*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][im*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][ip2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][im2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][ip3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][im3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+jp*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+jm*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+jp2*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+jm2*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+jp3*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+jm3*_NZ+k]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+j*_NZ+kp]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+j*_NZ+km]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+j*_NZ+kp2]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+j*_NZ+km2]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+j*_NZ+kp3]-0.5)<=0)
                          || ((shm_PHI[1][ind]-0.5)*(shm_PHI[1][i*_NY*_NZ+j*_NZ+km3]-0.5)<=0)) )
                { 
                    shm_K_S[ind] = 1;
                }
                else
                {
                    shm_K_S[ind] = 0;
                }
                
                /* define the liquid phase boundary */
                if (        (((shm_PHI[0][ind]-0.5)*(shm_PHI[0][ip*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][im*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][ip2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][im2*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][ip3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][im3*_NY*_NZ+j*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+jp*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+jm*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+jp2*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+jm2*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+jp3*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+jm3*_NZ+k]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+j*_NZ+kp]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+j*_NZ+km]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+j*_NZ+kp2]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+j*_NZ+km2]-0.5)<=0) 
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+j*_NZ+kp3]-0.5)<=0)
                             || ((shm_PHI[0][ind]-0.5)*(shm_PHI[0][i*_NY*_NZ+j*_NZ+km3]-0.5)<=0)))
                { 
                    shm_K_L[ind] = 1;
                }
                else
                {
                    shm_K_L[ind] = 0;
                }
                /* define the l-s interface */
                if (( shm_K_L[ind]==1 ) && ( shm_K_S[ind]==1 )) 
                {
                    shm_K_LS[ind] = 1;
                }
                else
                {
                    shm_K_LS[ind] = 0;
                }
                /* define the s-v interface */
                if (( shm_K_L[ind]==0 ) && ( shm_K_S[ind]==1 ))
                {
                    shm_K_SV[ind] = 1;
                }
                else
                {
                    shm_K_SV[ind] = 0;
                }
                /* define the l-v interface */
                if (( shm_K_L[ind]==1 ) && ( shm_K_S[ind]==0 ))
                {
                    shm_K_LV[ind] = 1;
                }
                else
                {
                    shm_K_LV[ind] = 0;
                }
                
                shm_K_other[ind] = 1 - shm_K_LV[ind] - shm_K_SV[ind] - shm_K_LS[ind];
                
            }
        }
    }
    
    
    
    
    double loc_F, loc_Fraw, sum_F, sum_Fraw;
    sum_F = 0; sum_Fraw = 0;    
#pragma omp parallel for reduction(+:sum_F,sum_Fraw) \
private(n,i,j,k,ip,im,jp,jm,kp,km,ind,loc_F,loc_Fraw,phi,grad_term,f,fp) \
private(EPS_loc,EPS_prefactor,dphidx,dphidy,dphidz,dphi_SQR,absdphi) \
private(nx_lab,ny_lab,nz_lab,nx_cryst,ny_cryst,nz_cryst) \
private(eps1,eps2,eps3,triplesum) \
private(dnx_lab_dphidx,dny_lab_dphidx,dnz_lab_dphidx) \
private(dnx_lab_dphidy,dny_lab_dphidy,dnz_lab_dphidy) \
private(dnx_lab_dphidz,dny_lab_dphidz,dnz_lab_dphidz) \
private(dEPS_dnx_cryst,dEPS_dny_cryst,dEPS_dnz_cryst) \
private(dEPS_dnx_lab,  dEPS_dny_lab,  dEPS_dnz_lab  ) \
private(dEPS_dphidx,   dEPS_dphidy,   dEPS_dphidz   ) \
private(dEPS_dphidx_loc,dEPS_dphidy_loc,dEPS_dphidz_loc ) \
shared(shm_dPHIdx,shm_dPHIdy,shm_dPHIdz,shm_PHI,shm_dFdPHIi,shm_d2PHI) \
shared(shm_tmp_x,shm_tmp_y,shm_tmp_z) \
shared(shm_matrix_x,shm_matrix_y,shm_EPS_loc,shm_NW_orig)   
    for (ind=0; ind<_NX*_NY*_NZ; ind++)
    {    
        loc_F = 0; loc_Fraw = 0;
        for(n=0;n<num_fields;n++)
        { 
            grad_term = SQR(shm_dPHIdx[n][ind]) + SQR(shm_dPHIdy[n][ind]) + SQR(shm_dPHIdz[n][ind]) ;
            
            f = SQR(shm_PHI[n][ind]) * SQR(1-shm_PHI[n][ind]);
            //fp = _PHI[n][ind];
            fp = tanh((shm_PHI[n][ind]-0.5)/0.10)/2 + 0.5;
            loc_F += f * _U[n][n] + fp * _MU[n];
            loc_Fraw += f * _U[n][n];
            /* the influence of the external force (only for liquid) */ 
            if ( n == 0 )
            { 
                loc_F += Fext_x * fp * shm_matrix_x[ind] + Fext_y * fp * shm_matrix_y[ind]; 
            }
            
            // Penalty term for three phase co-existing    
            loc_F += Penalty_3phase/num_fields*SQR(shm_PHI[0][ind])*SQR(shm_PHI[1][ind])*SQR(shm_PHI[2][ind]);
            
            // Penalty term for deviation from original NW shape
            loc_F += Penalty_NW/num_fields*SQR(shm_PHI[1][ind] - shm_NW_orig[ind]);
            
            if ( ((_EPS1[n][n] == 0)&&(_EPS2[n][n] == 0)&&(_EPS3[n][n] == 0)) )
            { /* isotropic interface energy */
                EPS_loc = 1.0;
                EPS_prefactor = SQR(EPS_loc);
                
                loc_F += _EPS[n][n] * EPS_prefactor * grad_term;
                loc_Fraw += _EPS[n][n] * EPS_prefactor * grad_term;
                /* the equation of motion should follow the Steinbach 1999 formulation */
                shm_dFdPHIi[n][ind] = _U[n][n]*(4.0*shm_PHI[n][ind]*shm_PHI[n][ind]*shm_PHI[n][ind]
                                                + 2.0*shm_PHI[n][ind] - 6.0*shm_PHI[n][ind]*shm_PHI[n][ind]) 
                + _MU[n]*(5.0-5.0*SQR(tanh(10*shm_PHI[n][ind]-5.0)))
                - 2.0*_EPS[n][n]*EPS_prefactor*shm_d2PHI[n][ind]; 
                
                /* with the choice of C(phi) = phi */  
                /* _dFdPHIi[n][ind] = _U[n][n]*(4.0*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
                 + 2.0*_PHI[n][ind] - 6.0*_PHI[n][ind]*_PHI[n][ind]) 
                 + _MU[n] - 2.0*_EPS[n][n]*EPS_prefactor*_d2PHI[n][ind]; */
                
                /* with the external force */
                if ( n==0 )
                { shm_dFdPHIi[0][ind] += Fext_x * (5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_matrix_x[ind]
                    + Fext_y * (5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_matrix_y[ind];
                }
                
            }
            else
            {   /* anisotropic interface energy */                          
                /* the EPS term will be computed and added below */
                /* surface normal defined by solid phase gradient */
                dphidx = shm_dPHIdx[n][ind];
                dphidy = shm_dPHIdy[n][ind];
                dphidz = shm_dPHIdz[n][ind];
                
                dphi_SQR   = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);
                absdphi = sqrt(dphi_SQR + 1e-24);
                
                nx_lab = dphidx / absdphi;
                ny_lab = dphidy / absdphi;
                nz_lab = dphidz / absdphi;
                
                nx_cryst = _R[0][0]*nx_lab + _R[0][1]*ny_lab + _R[0][2]*nz_lab;
                ny_cryst = _R[1][0]*nx_lab + _R[1][1]*ny_lab + _R[1][2]*nz_lab;
                nz_cryst = _R[2][0]*nx_lab + _R[2][1]*ny_lab + _R[2][2]*nz_lab;
                
                eps1 = _EPS1[n][n]; /* for solid(1) phase only */
                eps2 = _EPS2[n][n]; /* for solid(1) phase only */
                eps3 = _EPS3[n][n]; /* for solid(1) phase only */
                
                triplesum = SQR(nx_cryst)*SQR(ny_cryst) + SQR(ny_cryst)*SQR(nz_cryst) + SQR(nz_cryst)*SQR(nx_cryst);
                
                /* Notice the difference from the Matlab implementation */
                EPS_loc   = 1.0 + eps1*triplesum + eps2*SQR(nx_cryst*ny_cryst*nz_cryst) + eps3*SQR(triplesum);
                shm_EPS_loc[n][ind] = EPS_loc;
                
                /* Free energy for anisotropic surface energy */                
                loc_F += _EPS[n][n] * SQR(shm_EPS_loc[n][ind]) * grad_term;
                loc_Fraw += _EPS[n][n] * SQR(shm_EPS_loc[n][ind]) * grad_term;
                
                dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
                dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
                dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);
                
                dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
                dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
                dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);
                
                dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
                dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
                dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);
                
                dEPS_dnx_cryst = 2.*eps1*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))
                + 2.*eps2*nx_cryst*SQR(ny_cryst)*SQR(nz_cryst)
                + 4.*eps3*nx_cryst*(SQR(ny_cryst)+SQR(nz_cryst))*triplesum;
                dEPS_dny_cryst = 2.*eps1*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))
                + 2.*eps2*ny_cryst*SQR(nz_cryst)*SQR(nx_cryst)
                + 4.*eps3*ny_cryst*(SQR(nz_cryst)+SQR(nx_cryst))*triplesum;
                dEPS_dnz_cryst = 2.*eps1*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))
                + 2.*eps2*nz_cryst*SQR(nx_cryst)*SQR(ny_cryst)
                + 4.*eps3*nz_cryst*(SQR(nx_cryst)+SQR(ny_cryst))*triplesum;
                
                dEPS_dnx_lab = _R[0][0]*dEPS_dnx_cryst + _R[1][0]*dEPS_dny_cryst + _R[2][0]*dEPS_dnz_cryst;
                dEPS_dny_lab = _R[0][1]*dEPS_dnx_cryst + _R[1][1]*dEPS_dny_cryst + _R[2][1]*dEPS_dnz_cryst;
                dEPS_dnz_lab = _R[0][2]*dEPS_dnx_cryst + _R[1][2]*dEPS_dny_cryst + _R[2][2]*dEPS_dnz_cryst;
                
                dEPS_dphidx = dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx;
                dEPS_dphidy = dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy;
                dEPS_dphidz = dEPS_dnx_lab*dnx_lab_dphidz + dEPS_dny_lab*dny_lab_dphidz + dEPS_dnz_lab*dnz_lab_dphidz;
                

                dEPS_dphidx_loc = dEPS_dphidx;
                dEPS_dphidy_loc = dEPS_dphidy;
                dEPS_dphidz_loc = dEPS_dphidz;
 
                
                grad_loc_x = shm_dPHIdx[n][ind];
                grad_loc_y = shm_dPHIdy[n][ind];
                grad_loc_z = shm_dPHIdz[n][ind];
                grad_loc_SQR = SQR(grad_loc_x) + SQR(grad_loc_y) + SQR(grad_loc_z);
                
                /* new scheme for computing variational derivative, 2012/04/27  */
                shm_tmp_x[n][ind] = grad_loc_SQR*EPS_loc*dEPS_dphidx_loc;
                shm_tmp_y[n][ind] = grad_loc_SQR*EPS_loc*dEPS_dphidy_loc;
                shm_tmp_z[n][ind] = grad_loc_SQR*EPS_loc*dEPS_dphidz_loc;
                
                shm_dFdPHIi[n][ind] += _U[n][n]*(4*shm_PHI[n][ind]*shm_PHI[n][ind]*shm_PHI[n][ind]
                                                 + 2*shm_PHI[n][ind] - 6*shm_PHI[n][ind]*shm_PHI[n][ind]) 
                + _MU[n]*(5.0-5.0*SQR(tanh(10*shm_PHI[n][ind]-5.0))); 
                
                /* for the choice of C(phi) = phi */
                /* _dFdPHIi[n][ind] += _U[n][n]*(4*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
                 + 2*_PHI[n][ind] - 6*_PHI[n][ind]*_PHI[n][ind]) 
                 + _MU[n]; */
                
                shm_dFdPHIi[n][ind] += (-2.0) * _EPS[n][n]*SQR(EPS_loc) *(shm_d2PHI[n][ind]);
                
                // with the external force term
                if ( n==0 )
                { shm_dFdPHIi[0][ind] += Fext_x * (5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_matrix_x[ind]
                    + Fext_y * (5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_matrix_y[ind];
                }
            }
            
        }
        sum_F = sum_F + loc_F;
        sum_Fraw = sum_Fraw + loc_Fraw;
        
        shm_dFdPHIi[0][ind] += 2*Penalty_3phase*shm_PHI[0][ind]*SQR(shm_PHI[1][ind])*SQR(shm_PHI[2][ind]);
        shm_dFdPHIi[1][ind] += 2*Penalty_3phase*shm_PHI[1][ind]*SQR(shm_PHI[0][ind])*SQR(shm_PHI[2][ind]);
        shm_dFdPHIi[2][ind] += 2*Penalty_3phase*shm_PHI[2][ind]*SQR(shm_PHI[0][ind])*SQR(shm_PHI[1][ind]);
        
        /* variational derivative due to the NW shape penalty term */
        shm_dFdPHIi[1][ind] += 2*Penalty_NW*(shm_PHI[1][ind] - shm_NW_orig[ind]);
        
    }
    
    _F = sum_F*CUBE(gridsize);
    _Fraw = sum_Fraw*CUBE(gridsize);
       
    for(n=0;n<num_fields;n++)
    {    
	    if ((_EPS1[n][n] != 0)||(_EPS2[n][n] != 0)||(_EPS3[n][n] != 0))
	    {        
	#pragma omp parallel for default(none) \
	private(i,j,k,ip,im,jp,jm,kp,km,ind) \
	private(grad_loc_x,grad_loc_y,grad_loc_z) \
	shared(shm_dFdPHIi,shm_tmp_x,shm_tmp_y,shm_tmp_z,shm_dPHIdx,shm_dPHIdy,shm_dPHIdz,shm_EPS_loc,n) 
		
		/* cubic anisotropic interface energy only for solid phase */
		for(i=0;i<_NX;i++)
		{
		    ip = (i+1)%_NX;
		    im = (i-1+_NX)%_NX;
		    for (j=0;j<_NY;j++)
		    {
		        jp = (j+1)%_NY;
		        jm = (j-1+_NY)%_NY;
		        for (k=0;k<_NZ;k++)
		        {
		            kp = (k+1)%_NZ;
		            km = (k-1+_NZ)%_NZ;
		            
		            ind = i*_NY*_NZ+j*_NZ+k;
		            shm_dFdPHIi[n][ind] += _EPS[n][n]*(-2.0)*(shm_tmp_x[n][ip*_NY*_NZ+j*_NZ+k]-shm_tmp_x[n][im*_NY*_NZ+j*_NZ+k]) / (2.0*gridsize)
		            + _EPS[n][n]*(-2.0)*(shm_tmp_y[n][i*_NY*_NZ+jp*_NZ+k]-shm_tmp_y[n][i*_NY*_NZ+jm*_NZ+k]) / (2.0*gridsize)
		            + _EPS[n][n]*(-2.0)*(shm_tmp_z[n][i*_NY*_NZ+j*_NZ+kp]-shm_tmp_z[n][i*_NY*_NZ+j*_NZ+km]) / (2.0*gridsize);
		            
		            /* new scheme for computing variational derivative, 2012/04/27  */
		            grad_loc_x = shm_dPHIdx[n][ind];
		            grad_loc_y = shm_dPHIdy[n][ind];
		            grad_loc_z = shm_dPHIdz[n][ind];
		            
		            
		            shm_dFdPHIi[n][ind] += _EPS[n][n]*(-4.0)*shm_EPS_loc[n][ind]*(shm_EPS_loc[n][ip*_NY*_NZ+j*_NZ+k]
		                                                                          - shm_EPS_loc[n][im*_NY*_NZ+j*_NZ+k])/ (2.0*gridsize) * grad_loc_x
		            + _EPS[n][n]*(-4.0)*shm_EPS_loc[n][ind]*(shm_EPS_loc[n][i*_NY*_NZ+jp*_NZ+k]
		                                                     - shm_EPS_loc[n][i*_NY*_NZ+jm*_NZ+k])/ (2.0*gridsize) *grad_loc_y                                 
		            + _EPS[n][n]*(-4.0)*shm_EPS_loc[n][ind]*(shm_EPS_loc[n][i*_NY*_NZ+j*_NZ+kp]
		                                                     - shm_EPS_loc[n][i*_NY*_NZ+j*_NZ+km])/ (2.0*gridsize) * grad_loc_z;                
		            
		        }
		    }
		}
	    }
    }
    
    
    /* equation of motion based on Steinback 1999 formulation */  
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_dFdPHIi,shm_dPHIdt0,shm_dPHIdtCH) \
shared(shm_K_SV,shm_K_LV,shm_K_LS,shm_K_other,shm_K_D)  
    for (ind=0; ind<_NX*_NY*_NZ; ind++)
    {
        
        /* G-L equation of motion */
        shm_dPHIdt0[0][ind] = - shm_K_SV[ind]*(Mob_GL*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[1][ind]) 
                                               + Mob_GL*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[2][ind])) 
        - shm_K_LV[ind]*(Mob_GL*((shm_dFdPHIi[0][ind]+shm_dFdPHIi[2][ind])/2 - shm_dFdPHIi[1][ind]) 
                         + Mob_LV*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[2][ind])) 
        - shm_K_LS[ind]*(Mob_GL*((shm_dFdPHIi[0][ind]+shm_dFdPHIi[1][ind])/2 - shm_dFdPHIi[2][ind])
                         + Mob_LS*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[1][ind])) 
        - shm_K_other[ind]*(Mob_GL*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[1][ind])
                            + Mob_GL*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[2][ind]));
        
        shm_dPHIdt0[1][ind] = - shm_K_LV[ind]*(Mob_GL*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[2][ind])
                                               + Mob_GL*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[0][ind])) 
        - shm_K_SV[ind]*(Mob_GL*((shm_dFdPHIi[1][ind]+shm_dFdPHIi[2][ind])/2 - shm_dFdPHIi[0][ind])
                         + Mob_SV*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[2][ind]))
        - shm_K_LS[ind]*(Mob_GL*((shm_dFdPHIi[0][ind]+shm_dFdPHIi[1][ind])/2 - shm_dFdPHIi[2][ind])
                         + Mob_LS*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[0][ind])) 
        - shm_K_other[ind]*(Mob_GL*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[0][ind])
                            + Mob_GL*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[2][ind]));
        
        shm_dPHIdt0[2][ind] = - shm_K_LS[ind]*(Mob_GL*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[1][ind])
                                               + Mob_GL*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[0][ind])) 
        - shm_K_LV[ind]*(Mob_GL*((shm_dFdPHIi[0][ind]+shm_dFdPHIi[2][ind])/2 - shm_dFdPHIi[1][ind])
                         + Mob_LV*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[0][ind]))
        - shm_K_SV[ind]*(Mob_GL*((shm_dFdPHIi[1][ind]+shm_dFdPHIi[2][ind])/2 - shm_dFdPHIi[0][ind])
                         + Mob_SV*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[1][ind])) 
        - shm_K_other[ind]*(Mob_GL*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[0][ind])
                            + Mob_GL*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[1][ind]));
        
     	/* C-H equation of motion */	    
        shm_dPHIdtCH[0][ind] = - shm_K_D[ind]*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[1][ind])
        - shm_K_D[ind]*(shm_dFdPHIi[0][ind] - shm_dFdPHIi[2][ind]);
        
        shm_dPHIdtCH[1][ind] = - shm_K_D[ind]*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[0][ind])
        - shm_K_D[ind]*(shm_dFdPHIi[1][ind] - shm_dFdPHIi[2][ind]);
        
        shm_dPHIdtCH[2][ind] = - shm_K_D[ind]*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[0][ind])
        - shm_K_D[ind]*(shm_dFdPHIi[2][ind] - shm_dFdPHIi[1][ind]);    
		
    }                        
    
    /* calculate the max value of dphidt */
    
    double Gmiddle = 0;
if (dynamics_type != 1) 
{
#pragma omp parallel
    {
    	double priv_max = 0;
        
#pragma omp for nowait     
    	for (ind=0; ind<_NX*_NY*_NZ; ind++)   
    	{                
        	if ( priv_max <= fmax(fmax(fabs(shm_dPHIdt0[0][ind]), fabs(shm_dPHIdt0[1][ind])),fabs(shm_dPHIdt0[2][ind])) )
        	    priv_max = fmax(fmax(fabs(shm_dPHIdt0[0][ind]), fabs(shm_dPHIdt0[1][ind])),fabs(shm_dPHIdt0[2][ind]));
    	}
#pragma omp critical
        {
        	if ( priv_max >= Gmiddle)  Gmiddle = priv_max;
        }
    }
}
else if (dynamics_type == 1)
{
#pragma omp parallel
    {
    	double priv_max = 0;
        
#pragma omp for nowait     
    	for (ind=0; ind<_NX*_NY*_NZ; ind++)   
    	{                
        	if ( priv_max <= fmax(fmax(fabs(shm_dPHIdtCH[0][ind]), fabs(shm_dPHIdtCH[1][ind])),fabs(shm_dPHIdtCH[2][ind])) )
        	    priv_max = fmax(fmax(fabs(shm_dPHIdtCH[0][ind]), fabs(shm_dPHIdtCH[1][ind])),fabs(shm_dPHIdtCH[2][ind]));
    	}
#pragma omp critical
        {
        	if ( priv_max >= Gmiddle)  Gmiddle = priv_max;
        }
    }
}
    _G = Gmiddle;
            
    Phase_evolution();    
    
# if 0
    INFO_Printf("_F = %20.12e\n", _F);
    INFO_Printf("_G = %20.12e\n", _G);
    INFO_Printf("dmu1 = %20.12e _mprime[0] = %20.12e dmm = %20.12e\n", dmu1,_mprime[0],dmm);
    INFO_Printf("_dPHIdt0[0] = %20.12e _dPHIdt[0] = %20.12e \n", _dPHIdt0[0][0],_dPHIdt[0][0]);
    //INFO_Printf("_dPHIdx[0][0] = %20.12e \n", _dPHIdx[0][0]);
    write_array(_mprime,"_mprime.dat");
    write_array(_PHI[0],"PHI_0.dat");
    write_array(_PHI[1],"PHI_1.dat");
    write_array(_PHI[2],"PHI_2.dat");
    write_array(_d2PHI[0],"d2PHI_0.dat");
    write_array(_d2PHI[1],"d2PHI_1.dat");
    write_array(_d2PHI[2],"d2PHI_2.dat");
    write_array(_dPHIdt0[0],"dPHIdt0_0.dat");
    write_array(_dPHIdt0[1],"dPHIdt0_1.dat");
    write_array(_dPHIdt0[2],"dPHIdt0_2.dat");
    write_array(_dPHIdt[0],"dPHIdt_0.dat");
    write_array(_dPHIdx[0],"dPHIdx_0.dat"); 
    write_array(_dPHIdy[0],"dPHIdy_0.dat");
    write_array(_dPHIdz[0],"dPHIdz_0.dat");
    write_array(_dPHIdx[1],"dPHIdx_1.dat"); 
    write_array(_dPHIdy[1],"dPHIdy_1.dat");
    write_array(_dPHIdz[1],"dPHIdz_1.dat");
    write_array(_dPHIdx[2],"dPHIdx_2.dat"); 
    write_array(_dPHIdy[2],"dPHIdy_2.dat");
    write_array(_dPHIdz[2],"dPHIdz_2.dat");
    write_array(_K_S,"K_S.dat");
    write_array(_K_L,"K_L.dat");
    write_array(_K_LS,"K_LS.dat");
    write_array(_K_LV,"K_LV.dat");
    write_array(_K_SV,"K_SV.dat");
    write_array(_K_other,"K_other.dat");
    write_array(_dFdPHIi[0],"dFdPHIi_0.dat");  
    write_array(_dFdPHIi[1],"dFdPHIi_1.dat"); 
    write_array(_dFdPHIi[2],"dFdPHIi_2.dat"); 
    sleep(sleepseconds);
# endif
    
    
    return _F;        
    
}

void PhaseFieldFrame::Phase_evolution() 
{
    int n, p, i, j, k, ip, im, jp, jm, kp, km, ind;
    double h2,dmm, dmu1;
    double A_11, A_12, A_13, A_14, A_21, A_22, A_23, A_24, A_31, A_32, A_33, A_34, A_41, A_42, A_43, A_44;
    double B_1, B_2, B_3, B_4;
    double dmuL, dmuV, dfX, dfY, divider;
    double vol_incr, x_incr, y_incr;
    
#ifdef _OPENMP
    /* only local arrays can be shared by openmp threads */
    double *shm_PHI[MAXNUMFIELDS]; 
    double *shm_dPHIdt0[MAXNUMFIELDS], *shm_dPHIdt[MAXNUMFIELDS], *shm_dPHIdtCH[MAXNUMFIELDS], *shm_d2dPHIdt0[MAXNUMFIELDS];
    double *shm_K_S, *shm_K_L, *shm_K_LS, *shm_K_LV, *shm_K_SV, *shm_K_other, *shm_K_D;
    double *shm_matrix_x, *shm_matrix_y;
    double *shm_mprime;
    double *shm_dGLdmuL, *shm_dGLdmuV, *shm_dGVdmuL, *shm_dGVdmuV, *shm_dGSdmuL, *shm_dGSdmuV;
    
    for(n=0;n<num_fields;n++)
    {
        shm_PHI[n]   = _PHI[n];  
        shm_dPHIdt0[n] = _dPHIdt0[n];
        shm_dPHIdt[n]  = _dPHIdt[n];
        shm_dPHIdtCH[n]  = _dPHIdtCH[n];
        shm_d2dPHIdt0[n] = _d2dPHIdt0[n];
    }
    
    shm_K_S = _K_S;
    shm_K_L = _K_L;
    shm_K_LS = _K_LS;
    shm_K_LV = _K_LV;
    shm_K_SV = _K_SV;
    shm_K_other = _K_other;
    shm_K_D = _K_D;
    shm_matrix_x = matrix_x;
    shm_matrix_y = matrix_y;
    shm_mprime = _mprime;
    shm_dGLdmuL = _dGLdmuL;
    shm_dGLdmuV = _dGLdmuV;
    shm_dGVdmuL = _dGVdmuL;
    shm_dGVdmuV = _dGVdmuV;
    shm_dGSdmuL = _dGSdmuL;
    shm_dGSdmuV = _dGSdmuV;
    
#else
    
    /* make the code compilable when not using -openmp */
#define shm_PHI     _PHI
#define shm_dPHIdt   _dPHIdt
#define shm_dPHIdt0  _dPHIdt0
#define shm_d2dPHIdt0 _d2dPHIdt0
#define shm_dPHIdtCH _dPHIdtCH
#define shm_K_L _K_S
#define shm_K_LS _K_LS
#define shm_K_LV _K_LV
#define shm_K_SV _K_SV
#define shm_K_other _K_other
#define shm_K_D _K_D
#define shm_matrix_x matrix_x
#define shm_matrix_y matrix_y
#define shm_mprime _mprime     
#define shm_dGLdmuL _dGLdmuL
#define shm_dGLdmuV _dGLdmuV
#define shm_dGVdmuL _dGVdmuL
#define shm_dGVdmuV _dGVdmuV
#define shm_dGSdmuL _dGSdmuL
#define shm_dGSdmuV _dGSdmuV 
#endif

h2 = gridsize*gridsize;    
/* adjust time step at every step */
timestep = dtmax/fmax(_G*dtmax/dphimax,1);
    
vol_incr = vol_incr0/timestep;
x_incr = x_incr0*liquid_volume/timestep;
y_incr = y_incr0*liquid_volume/timestep;    
    
    
    
if (dynamics_type != 2 && dynamics_type != 3 && dynamics_type != 4 && dynamics_type != 5 && dynamics_type != 6 && dynamics_type != 7 && dynamics_type != 8) 
    {
        
        if (dynamics_type == 1) /* Cahn-Hillard equation */
        {
        /*   #pragma omp parallel for default(none) \
            private(n,i,j,k,ip,im,jp,jm,kp,km) \
            shared(h2,shm_d2dPHIdt0,shm_dPHIdtCH)  */ 
            for(i=0;i<_NX;i++)
            {
                ip = (i+1)%_NX;
                im = (i-1+_NX)%_NX;
                for (j=0;j<_NY;j++)
                {
                    jp = (j+1)%_NY;
                    jm = (j-1+_NY)%_NY;
                    
                    for (k=0;k<_NZ;k++)
                    {
                        kp = (k+1)%_NZ;
                        km = (k-1+_NZ)%_NZ;
                        
                        for(n=0;n<num_fields;n++)
                        {
                            shm_d2dPHIdt0[n][i*_NY*_NZ+j*_NZ+k] = 
                            (shm_dPHIdtCH[n][ip*_NY*_NZ+j*_NZ+k] + shm_dPHIdtCH[n][im*_NY*_NZ+j*_NZ+k]
                             +shm_dPHIdtCH[n][i*_NY*_NZ+jp*_NZ+k] + shm_dPHIdtCH[n][i*_NY*_NZ+jm*_NZ+k]
                             +shm_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+kp] + shm_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+km]
                             -shm_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+k] * 6.0 ) / h2;
                        }
                    }
                }
            }
        }
                      
        if (dynamics_type == 0) /* no additional constraints */
        {   
            #pragma omp parallel for default(none) \
            private(ind) \
            shared(shm_dPHIdt,shm_dPHIdt0)
            /* 0: Ginzburg-Landau equation: not conserved */
            for (ind=0; ind<_NX*_NY*_NZ; ind++)
            {    
                shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind];
                shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind];
                shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind];
            }
        }
        else if (dynamics_type == 1) /* no additional constraints */
        {
   /*         #pragma omp parallel for default(none) \
            private(ind) \
            shared(shm_dPHIdt,shm_d2dPHIdt0) b */
            /* 1: Cahn-Hilliard equation: conserved */
            for (ind=0; ind<_NX*_NY*_NZ; ind++)
            {
                shm_dPHIdt[0][ind] = -shm_d2dPHIdt0[0][ind];
                shm_dPHIdt[1][ind] = -shm_d2dPHIdt0[1][ind];
                shm_dPHIdt[2][ind] = -shm_d2dPHIdt0[2][ind];
            }
        }
    }
    
        else if (dynamics_type == 2) 
    {   
        /* Update _mprime */
        double loc_dmm;
        dmm = 0; 
        #pragma omp parallel for reduction(+:dmm) \
        private(loc_dmm,ind) \
        shared(shm_mprime,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other,shm_PHI)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            loc_dmm = 0;
            shm_mprime[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            loc_dmm = SQR(shm_mprime[ind]);
            dmm = dmm + loc_dmm;
        }
        
        /* Update DM */	
        double loc_dmu1;
        dmu1 = 0;
        #pragma omp parallel for reduction(+:dmu1) \
        private(loc_dmu1,ind) \
        shared(shm_mprime,shm_dPHIdt0,dmm)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)    
        {    
            loc_dmu1 = (shm_mprime[ind] * shm_dPHIdt0[0][ind]) / dmm;
            dmu1 = dmu1 - loc_dmu1;
        }
        
        #pragma omp parallel for default(none) \
        private(ind) \
        shared(dmu1,shm_PHI,shm_dPHIdt,shm_dPHIdt0,shm_mprime,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other)
        /* Update gradient due to the change of m */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmu1* shm_mprime[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmu1
            * (shm_K_LS[ind]*(-Mob_GL/2+Mob_LS)+shm_K_LV[ind]*Mob_GL+shm_K_SV[ind]*Mob_GL+shm_K_other[ind]*Mob_GL)
            *(5.0-5.0*SQR(tanh(10.0*shm_PHI[0][ind]-5.0)));
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmu1
            * (shm_K_LS[ind]*Mob_GL+shm_K_LV[ind]*(-Mob_GL/2+Mob_LV)+shm_K_SV[ind]*Mob_GL+shm_K_other[ind]*Mob_GL)
            *(5.0-5.0*SQR(tanh(10.0*shm_PHI[0][ind]-5.0)));                     
        }
        /* update chemical potential */
        _MU[0] = _MU[1] + dmu1;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    
    else if (dynamics_type == 3)
    {       
        /* Update _mprime */
        double loc_dmm;
        dmm = 0; 
        #pragma omp parallel for reduction(+:dmm) \
        private(loc_dmm,ind) \
        shared(shm_mprime,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other,shm_PHI)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            loc_dmm = 0;
            shm_mprime[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            loc_dmm = SQR(shm_mprime[ind]);
            dmm = dmm + loc_dmm;
        }
        
        /* Update DM */	
        double loc_dmu1;
        dmu1 = 0;
        #pragma omp parallel for reduction(+:dmu1) \
        private(loc_dmu1,ind) \
        shared(shm_mprime,shm_dPHIdt0,dmm)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)    
        {    
            loc_dmu1 = (shm_mprime[ind] * shm_dPHIdt0[0][ind]) / dmm;
            dmu1 = dmu1 - loc_dmu1;
        }
        
        #pragma omp parallel for default(none) \
        private(ind) \
        shared(dmu1,shm_PHI,shm_dPHIdt,shm_dPHIdt0,shm_mprime,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other)
        /* Update gradient due to the change of m */
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmu1* shm_mprime[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmu1
            * (shm_K_LS[ind]*(-Mob_GL/2+Mob_LS)+shm_K_LV[ind]*Mob_GL+shm_K_SV[ind]*Mob_GL+shm_K_other[ind]*Mob_GL)
            *(5.0-5.0*SQR(tanh(10.0*shm_PHI[0][ind]-5.0)));
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmu1
            * (shm_K_LS[ind]*Mob_GL+shm_K_LV[ind]*(-Mob_GL/2+Mob_LV)+shm_K_SV[ind]*Mob_GL+shm_K_other[ind]*Mob_GL)
            *(5.0-5.0*SQR(tanh(10.0*shm_PHI[0][ind]-5.0)));                     
        }
        /* update chemical potential */
        _MU[0] = _MU[1] + dmu1;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
        
        #pragma omp parallel for default(none) \
        private(n,i,j,k,ip,im,jp,jm,kp,km) \
        shared(h2,shm_d2dPHIdt0,shm_dPHIdtCH)
        for(i=0;i<_NX;i++)
        {
            ip = (i+1)%_NX;
            im = (i-1+_NX)%_NX;
            for (j=0;j<_NY;j++)
            {
                jp = (j+1)%_NY;
                jm = (j-1+_NY)%_NY;
                
                for (k=0;k<_NZ;k++)
                {
                    kp = (k+1)%_NZ;
                    km = (k-1+_NZ)%_NZ;
                    
                    for(n=0;n<num_fields;n++)
                    {
                        shm_d2dPHIdt0[n][i*_NY*_NZ+j*_NZ+k] = 
                        (shm_dPHIdtCH[n][ip*_NY*_NZ+j*_NZ+k] + shm_dPHIdtCH[n][im*_NY*_NZ+j*_NZ+k]
                         +shm_dPHIdtCH[n][i*_NY*_NZ+jp*_NZ+k] + shm_dPHIdtCH[n][i*_NY*_NZ+jm*_NZ+k]
                         +shm_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+kp] + shm_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+km]
                         -shm_dPHIdtCH[n][i*_NY*_NZ+j*_NZ+k] * 6.0 ) / h2;
                    }
                }
            }
        }
        
        #pragma omp parallel for default(none) \
        private(ind) \
        shared(shm_dPHIdt,shm_d2dPHIdt0)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            shm_dPHIdt[0][ind] = shm_dPHIdt[0][ind] - shm_d2dPHIdt0[0][ind];                      
            shm_dPHIdt[1][ind] = shm_dPHIdt[1][ind] - shm_d2dPHIdt0[1][ind];                      
            shm_dPHIdt[2][ind] = shm_dPHIdt[2][ind] - shm_d2dPHIdt0[2][ind];  
        }
    }
    
    else if (dynamics_type == 4) /* constrain all phase volume */
    {
       _MU[1] = 0;
        #pragma omp parallel for default(none) \
        private(ind) \
        shared(shm_PHI,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other) \
        shared(shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            shm_dGLdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGLdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGVdmuL[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));       
            
            shm_dGVdmuV[ind] = - shm_K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGSdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGSdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }  

        /* initialize the parameters */
        A_11=0, A_12=0, A_21=0, A_22=0, B_1=0, B_2=0, dmuL=0, dmuV=0;
        #pragma omp parallel for reduction(+:A_11,A_12,A_21) \
        private(ind) \
        shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_12 += shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_21 += shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }

        #pragma omp parallel for reduction(+:A_22,B_1,B_2) \
        private(ind) \
        shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {        
            A_22 += shm_dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0))); 
            B_1  += -(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
            B_2  += -(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))*shm_dPHIdt0[2][ind];
        }
       
       B_2 += vol_incr;
        
       dmuL = (A_12*B_2 - A_22*B_1)/(A_12*A_21-A_11*A_22);
       dmuV = -(A_11*B_2 - A_21*B_1)/(A_12*A_21-A_11*A_22);
        
        #pragma omp parallel for default(none) \
        private(ind) \
        shared(shm_dPHIdt,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
        shared(dmuL,dmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmuL*shm_dGLdmuL[ind] + dmuV*shm_dGLdmuV[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmuL*shm_dGSdmuL[ind] + dmuV*shm_dGSdmuV[ind];
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmuL*shm_dGVdmuL[ind] + dmuV*shm_dGVdmuV[ind];
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    
    else if (dynamics_type == 5) // constrain y-center of mass and all phase volumes
    {   
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_PHI,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other) \
shared(shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            shm_dGLdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGLdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGVdmuL[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));       
            
            shm_dGVdmuV[ind] = - shm_K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGSdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGSdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }
        
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
        B_1=0; B_2=0; B_3=0; dmuL=0; dmuV=0; dfY=0; divider=0; _MU[1]=0;
        
#pragma omp parallel for reduction(+:A_11,A_12,A_13) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_y)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_12 += shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_13 += shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
        
#pragma omp parallel for reduction(+:A_21,A_22,A_23) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_y)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_21 += shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            A_22 += shm_dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0))); 
            A_23 += shm_matrix_y[ind]*shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }
        
#pragma omp parallel for reduction(+:A_31,A_32,A_33) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_y)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_31 += shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_32 += shm_matrix_y[ind]*shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_33 += SQR(shm_matrix_y[ind])*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
#pragma omp parallel for reduction(+:B_1,B_2,B_3) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_y)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            B_1 += -(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
            B_2 += -(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))*shm_dPHIdt0[2][ind];
            B_3 += -shm_matrix_y[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
        }        
        
        B_2 += vol_incr;
        B_3 += y_incr;
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
        
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_dPHIdt,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_y) \
shared(dmuL,dmuV,dfY)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmuL*shm_dGLdmuL[ind] + dmuV*shm_dGLdmuV[ind] + dfY*shm_dGLdmuL[ind]*shm_matrix_y[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmuL*shm_dGSdmuL[ind] + dmuV*shm_dGSdmuV[ind] + dfY*shm_dGSdmuL[ind]*shm_matrix_y[ind];
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmuL*shm_dGVdmuL[ind] + dmuV*shm_dGVdmuV[ind] + dfY*shm_dGVdmuL[ind]*shm_matrix_y[ind];
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_y = 0.0 + dfY;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
  
    else if (dynamics_type == 6) // constrain x-center of mass and all phase volumes
    {   
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_PHI,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other) \
shared(shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            shm_dGLdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGLdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGVdmuL[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));       
            
            shm_dGVdmuV[ind] = - shm_K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGSdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGSdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
        B_1=0; B_2=0; B_3=0; dmuL=0; dmuV=0; dfX=0; divider=0; _MU[1]=0;
#pragma omp parallel for reduction(+:A_11,A_12,A_13) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x)        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_12 += shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_13 += shm_matrix_x[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
         }
#pragma omp parallel for reduction(+:A_21,A_22,A_23) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x)        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_21 += shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            A_22 += shm_dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0))); 
            A_23 += shm_matrix_x[ind]*shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }
#pragma omp parallel for reduction(+:A_31,A_32,A_33) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x)        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_31 += shm_matrix_x[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_32 += shm_matrix_x[ind]*shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_33 += SQR(shm_matrix_x[ind])*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
#pragma omp parallel for reduction(+:B_1,B_2,B_3) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x)        
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            B_1 += -(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
            B_2 += -(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))*shm_dPHIdt0[2][ind];
            B_3 += -shm_matrix_x[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
        }
        
        B_2 += vol_incr;
        B_3 += x_incr;
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfX  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
        
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_dPHIdt,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x) \
shared(dmuL,dmuV,dfX)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmuL*shm_dGLdmuL[ind] + dmuV*shm_dGLdmuV[ind] + dfX*shm_dGLdmuL[ind]*shm_matrix_x[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmuL*shm_dGSdmuL[ind] + dmuV*shm_dGSdmuV[ind] + dfX*shm_dGSdmuL[ind]*shm_matrix_x[ind];
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmuL*shm_dGVdmuL[ind] + dmuV*shm_dGVdmuV[ind] + dfX*shm_dGVdmuL[ind]*shm_matrix_x[ind];            
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_x = 0.0 + dfX;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
   
    else if (dynamics_type == 7) // constrain center of mass and all phase volumes 
    {   
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_PHI,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other) \
shared(shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            shm_dGLdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGLdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGVdmuL[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));       
            
            shm_dGVdmuV[ind] = - shm_K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGSdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGSdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }
        
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_14=0; A_21=0; A_22=0; A_23=0; A_24=0; A_31=0; A_32=0; A_33=0; A_34=0; A_41=0; A_42=0; A_43=0; A_44=0;
        B_1=0; B_2=0; B_3=0; B_4=0; dmuL=0; dmuV=0; dfX=0; dfY=0; divider=0; _MU[1]=0;
#pragma omp parallel for reduction(+:A_11,A_12,A_13) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_12 += shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_13 += shm_matrix_x[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
#pragma omp parallel for reduction(+:A_14,A_21,A_22) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_14 += shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_21 += shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            A_22 += shm_dGVdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0))); 
        }
#pragma omp parallel for reduction(+:A_23,A_24,A_31) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_23 += shm_matrix_x[ind]*shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            A_24 += shm_matrix_y[ind]*shm_dGVdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            A_31 += shm_matrix_x[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
#pragma omp parallel for reduction(+:A_32,A_33,A_34) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_32 += shm_matrix_x[ind]*shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_33 += SQR(shm_matrix_x[ind])*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_34 += shm_matrix_x[ind]*shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
#pragma omp parallel for reduction(+:A_41,A_42,A_43) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_41 += shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_42 += shm_matrix_y[ind]*shm_dGLdmuV[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_43 += shm_matrix_x[ind]*shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
         }
#pragma omp parallel for reduction(+:A_44,B_1,B_2) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {     
            A_44 += SQR(shm_matrix_y[ind])*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            B_1 += -(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
            B_2 += -(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))*shm_dPHIdt0[2][ind];
        }

#pragma omp parallel for reduction(+:B_3,B_4) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            B_3 += -shm_matrix_x[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
            B_4 += -shm_matrix_y[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
        }
        
        B_2 += vol_incr;
        B_3 += x_incr;
        B_4 += y_incr;
        
        divider = A_11*A_22*A_33*A_44-A_11*A_22*A_34*A_43-A_11*A_23*A_32*A_44+A_11*A_23*A_34*A_42+A_11*A_24*A_32*A_43
        - A_11*A_24*A_33*A_42-A_12*A_21*A_33*A_44+A_12*A_21*A_34*A_43+A_12*A_23*A_31*A_44-A_12*A_23*A_34*A_41
        - A_12*A_24*A_31*A_43+A_12*A_24*A_33*A_41+A_13*A_21*A_32*A_44-A_13*A_21*A_34*A_42-A_13*A_22*A_31*A_44
        + A_13*A_22*A_34*A_41+A_13*A_24*A_31*A_42-A_13*A_24*A_32*A_41-A_14*A_21*A_32*A_43+A_14*A_21*A_33*A_42
        + A_14*A_22*A_31*A_43-A_14*A_22*A_33*A_41-A_14*A_23*A_31*A_42+A_14*A_23*A_32*A_41;
        dmuL = -(A_12*A_23*A_34*B_4-A_12*A_23*A_44*B_3-A_12*A_24*A_33*B_4+A_12*A_24*A_43*B_3
                 +A_12*A_33*A_44*B_2-A_12*A_34*A_43*B_2-A_13*A_22*A_34*B_4+A_13*A_22*A_44*B_3
                 +A_13*A_24*A_32*B_4-A_13*A_24*A_42*B_3-A_13*A_32*A_44*B_2+A_13*A_34*A_42*B_2
                 +A_14*A_22*A_33*B_4-A_14*A_22*A_43*B_3-A_14*A_23*A_32*B_4+A_14*A_23*A_42*B_3
                 +A_14*A_32*A_43*B_2-A_14*A_33*A_42*B_2-A_22*A_33*A_44*B_1+A_22*A_34*A_43*B_1
                 +A_23*A_32*A_44*B_1-A_23*A_34*A_42*B_1-A_24*A_32*A_43*B_1+A_24*A_33*A_42*B_1)/divider;
        dmuV =  (A_11*A_23*A_34*B_4-A_11*A_23*A_44*B_3-A_11*A_24*A_33*B_4+A_11*A_24*A_43*B_3
                 +A_11*A_33*A_44*B_2-A_11*A_34*A_43*B_2-A_13*A_21*A_34*B_4+A_13*A_21*A_44*B_3
                 +A_13*A_24*A_31*B_4-A_13*A_24*A_41*B_3-A_13*A_31*A_44*B_2+A_13*A_34*A_41*B_2
                 +A_14*A_21*A_33*B_4-A_14*A_21*A_43*B_3-A_14*A_23*A_31*B_4+A_14*A_23*A_41*B_3
                 +A_14*A_31*A_43*B_2-A_14*A_33*A_41*B_2-A_21*A_33*A_44*B_1+A_21*A_34*A_43*B_1
                 +A_23*A_31*A_44*B_1-A_23*A_34*A_41*B_1-A_24*A_31*A_43*B_1+A_24*A_33*A_41*B_1)/divider;
        dfX  = -(A_11*A_22*A_34*B_4-A_11*A_22*A_44*B_3-A_11*A_24*A_32*B_4+A_11*A_24*A_42*B_3
                 +A_11*A_32*A_44*B_2-A_11*A_34*A_42*B_2-A_12*A_21*A_34*B_4+A_12*A_21*A_44*B_3
                 +A_12*A_24*A_31*B_4-A_12*A_24*A_41*B_3-A_12*A_31*A_44*B_2+A_12*A_34*A_41*B_2
                 +A_14*A_21*A_32*B_4-A_14*A_21*A_42*B_3-A_14*A_22*A_31*B_4+A_14*A_22*A_41*B_3
                 +A_14*A_31*A_42*B_2-A_14*A_32*A_41*B_2-A_21*A_32*A_44*B_1+A_21*A_34*A_42*B_1
                 +A_22*A_31*A_44*B_1-A_22*A_34*A_41*B_1-A_24*A_31*A_42*B_1+A_24*A_32*A_41*B_1)/divider;
        dfY  =  (A_11*A_22*A_33*B_4-A_11*A_22*A_43*B_3-A_11*A_23*A_32*B_4+A_11*A_23*A_42*B_3
                 +A_11*A_32*A_43*B_2-A_11*A_33*A_42*B_2-A_12*A_21*A_33*B_4+A_12*A_21*A_43*B_3
                 +A_12*A_23*A_31*B_4-A_12*A_23*A_41*B_3-A_12*A_31*A_43*B_2+A_12*A_33*A_41*B_2
                 +A_13*A_21*A_32*B_4-A_13*A_21*A_42*B_3-A_13*A_22*A_31*B_4+A_13*A_22*A_41*B_3
                 +A_13*A_31*A_42*B_2-A_13*A_32*A_41*B_2-A_21*A_32*A_43*B_1+A_21*A_33*A_42*B_1
                 +A_22*A_31*A_43*B_1-A_22*A_33*A_41*B_1-A_23*A_31*A_42*B_1+A_23*A_32*A_41*B_1)/divider;
        
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_dPHIdt,shm_dPHIdt0,shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV) \
shared(shm_matrix_x,shm_matrix_y) \
shared(dmuL,dmuV,dfX,dfY)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmuL*shm_dGLdmuL[ind] + dmuV*shm_dGLdmuV[ind] 
            + dfX*shm_dGLdmuL[ind]*shm_matrix_x[ind] + dfY*shm_dGLdmuL[ind]*shm_matrix_y[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmuL*shm_dGSdmuL[ind] + dmuV*shm_dGSdmuV[ind]
            + dfX*shm_dGSdmuL[ind]*shm_matrix_x[ind] + dfY*shm_dGSdmuL[ind]*shm_matrix_y[ind];
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmuL*shm_dGVdmuL[ind] + dmuV*shm_dGVdmuV[ind]
            + dfX*shm_dGVdmuL[ind]*shm_matrix_x[ind] + dfY*shm_dGVdmuL[ind]*shm_matrix_y[ind];
            
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_x = 0.0 + dfX;
        Fext_y = 0.0 + dfY;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        }
    }
    else if (dynamics_type == 8)  // constrain center of mass and droplet volume
    {   
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_PHI,shm_K_LS,shm_K_SV,shm_K_LV,shm_K_other) \
shared(shm_dGLdmuL,shm_dGLdmuV,shm_dGVdmuL,shm_dGVdmuV,shm_dGSdmuL,shm_dGSdmuV)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {
            shm_dGLdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGLdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGVdmuL[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));       
            
            shm_dGVdmuV[ind] = - shm_K_LS[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
            
            shm_dGSdmuL[ind] = - shm_K_LS[ind]*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_SV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            
            shm_dGSdmuV[ind] = - shm_K_LS[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_SV[ind]*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_LV[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)))
            - shm_K_other[ind]*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*shm_PHI[2][ind]-5.0)));
        }
        
        /* initialize the parameters */
        A_11=0; A_12=0; A_13=0; A_21=0; A_22=0; A_23=0; A_31=0; A_32=0; A_33=0;
        B_1=0; B_2=0; B_3=0; dmuL=0; dfX=0; dfY=0; divider=0; _MU[1]=0;
        
#pragma omp parallel for reduction(+:A_11,A_12,A_13) \
private(ind) \
shared(shm_PHI,shm_dGLdmuL,shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
            A_11 += shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_12 += shm_matrix_x[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_13 += shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));            
        }
      
#pragma omp parallel for reduction(+:A_21,A_22,A_23) \
private(ind) \
shared(shm_PHI,shm_dGLdmuL,shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {              
            A_21 += shm_matrix_x[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_22 += SQR(shm_matrix_x[ind])*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
            A_23 += shm_matrix_x[ind]*shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
        }
 
#pragma omp parallel for reduction(+:A_31,A_32,A_33) \
private(ind) \
shared(shm_PHI,shm_dGLdmuL,shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {   
              A_31 += shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
              A_32 += shm_matrix_x[ind]*shm_matrix_y[ind]*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));
              A_33 += SQR(shm_matrix_y[ind])*shm_dGLdmuL[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)));            
        } 
        
#pragma omp parallel for reduction(+:B_1,B_2,B_3) \
private(ind) \
shared(shm_PHI,shm_dPHIdt0,shm_dGLdmuL,shm_dGVdmuL,shm_dGSdmuL) \
shared(shm_matrix_x,shm_matrix_y) 
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {              
            B_1 += -(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
            B_2 += -shm_matrix_x[ind]*(5.0-5.0*SQR(tanh(10*shm_PHI[0][ind]-5.0)))*shm_dPHIdt0[0][ind];
            B_3 += -matrix_y[ind]*(5.0-5.0*SQR(tanh(10*_PHI[0][ind]-5.0)))*_dPHIdt0[0][ind];
        } 
      
        B_2 += x_incr;
        B_3 += y_incr;
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dfX  = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;

       
        
#pragma omp parallel for default(none) \
private(ind) \
shared(shm_dPHIdt,shm_dPHIdt0,shm_dGLdmuL,shm_dGVdmuL,shm_dGSdmuL) \
shared(shm_matrix_x,shm_matrix_y) \
shared(dmuL,dfX,dfY)
        for (ind=0; ind<_NX*_NY*_NZ; ind++)
        {           
            shm_dPHIdt[0][ind] = shm_dPHIdt0[0][ind] + dmuL*shm_dGLdmuL[ind] 
            + dfX*shm_dGLdmuL[ind]*shm_matrix_x[ind] + dfY*shm_dGLdmuL[ind]*shm_matrix_y[ind];
            
            shm_dPHIdt[1][ind] = shm_dPHIdt0[1][ind] + dmuL*shm_dGSdmuL[ind] 
            + dfX*shm_dGSdmuL[ind]*shm_matrix_x[ind] + dfY*shm_dGSdmuL[ind]*shm_matrix_y[ind];
            
            shm_dPHIdt[2][ind] = shm_dPHIdt0[2][ind] + dmuL*shm_dGVdmuL[ind] 
            + dfX*shm_dGVdmuL[ind]*shm_matrix_x[ind] + dfY*shm_dGVdmuL[ind]*shm_matrix_y[ind];
        }
        
        /* update chemical potential */
        _MU[0] = _MU[1] + dmuL;
        Fext_x = 0.0 + dfX;
        Fext_y = 0.0 + dfY;
        
        for (i=0;i<3;i++) {
            for (j=0;j<3;j++) {
                _M[i][j] = _MU[i]-_MU[j];
            }
        } 
    }    
  
    else
    {
        ERROR("unknown dynamics_type = "<<dynamics_type);
    }


}

