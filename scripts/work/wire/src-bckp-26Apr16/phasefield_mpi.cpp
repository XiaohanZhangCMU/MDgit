/* traditional (single) phase field model */
double PhaseFieldFrame::FreeEnergy_single_mpi()
{
	_F = 0;
        _Fraw = 0;

	if (gridsize <= 0)
	{
		FATAL("FreeEnergy: gridesize = " << gridsize);
	}

	// XXX for _NZ
	const int z_begin = _idx->get_range_x().begin();
	const int z_end   = _idx->get_range_x().end();

	// XXX for _NY
	const int y_begin = _idx->get_range_y().begin();
	const int y_end   = _idx->get_range_y().end();

	// XXX for _NX
	const int x_begin = _idx->get_range_z().begin();
	const int x_end   = _idx->get_range_z().end();

	for (int n = 0; n < num_fields; n++)
	{
		(s_PHI[n])->set_boundary_periodic();
	}

	/* compute spatial gradients */
	double h2 = gridsize * gridsize;
	for (int n = 0; n < num_fields; n++)
	{
		double local_F = 0.0;

#pragma omp parallel for reduction(+:local_F)
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					double phi = (*s_PHI[n])(k, j, i);
					double d2phi = ((*s_PHI[n])(k,     j,     i + 1) + (*s_PHI[n])(k,     j,     i - 1)
					             +  (*s_PHI[n])(k,     j + 1, i    ) + (*s_PHI[n])(k,     j - 1, i    )
					             +  (*s_PHI[n])(k + 1, j,     i    ) + (*s_PHI[n])(k - 1, j,     i    )
					             -  (*s_PHI[n])(k,     j,     i    ) * 6.0) / h2;

					double dphidx = ((*s_PHI[n])(k,     j,     i + 1) - (*s_PHI[n])(k,     j,     i - 1)) / (2.0 * gridsize);
					double dphidy = ((*s_PHI[n])(k,     j + 1, i    ) - (*s_PHI[n])(k,     j - 1, i    )) / (2.0 * gridsize);
					double dphidz = ((*s_PHI[n])(k + 1, j,     i    ) - (*s_PHI[n])(k - 1, j,     i    )) / (2.0 * gridsize);

					(*s_dPHIdx[n])(k, j, i) = dphidx;
					(*s_dPHIdy[n])(k, j, i) = dphidy;
					(*s_dPHIdz[n])(k, j, i) = dphidz;

					local_F += (phi * phi * (1 - phi) * (1 - phi)) * _U[n][n]; 
					(*s_dFdPHI[n])(k, j, i) = (2.0 * phi * (1 - phi) * (1 - 2.0 * phi)) * _U[n][n];

					/* isotropic interface energy */
					if ((_EPS1[n][n] == 0) && (_EPS2[n][n] == 0) && (_EPS3[n][n] == 0))
					{
						double EPS_prefactor = 0.5*SQR(_EPS[n][n]);

						local_F += (dphidx * dphidx + dphidy * dphidy + dphidz * dphidz) * EPS_prefactor;
						(*s_dFdPHI[n])(k, j, i) += - 2.0 * d2phi * EPS_prefactor;
					}
					else
					/* cubic anisotropic interface energy */
					{
						double dphi_SQR   = SQR(dphidx) + SQR(dphidy) + SQR(dphidz);
						double absdphi = sqrt(dphi_SQR + 1e-24);

#define _R _ROT_MATRIX 
						double nx_cryst = _R[0][0] * (dphidx / absdphi) + _R[0][1] * (dphidy / absdphi) + _R[0][2] * (dphidz / absdphi);
						double ny_cryst = _R[1][0] * (dphidx / absdphi) + _R[1][1] * (dphidy / absdphi) + _R[1][2] * (dphidz / absdphi);
						double nz_cryst = _R[2][0] * (dphidx / absdphi) + _R[2][1] * (dphidy / absdphi) + _R[2][2] * (dphidz / absdphi);

						double dnx_lab_dphidx = 1./absdphi - dphidx * dphidx / CUBE(absdphi);
						double dny_lab_dphidx =            - dphidy * dphidx / CUBE(absdphi);
						double dnz_lab_dphidx =            - dphidz * dphidx / CUBE(absdphi);

						double dnx_lab_dphidy =            - dphidx * dphidy / CUBE(absdphi);
						double dny_lab_dphidy = 1./absdphi - dphidy * dphidy / CUBE(absdphi);
						double dnz_lab_dphidy =            - dphidz * dphidy / CUBE(absdphi);

						double dnx_lab_dphidz =            - dphidx * dphidz / CUBE(absdphi);
						double dny_lab_dphidz =            - dphidy * dphidz / CUBE(absdphi);
						double dnz_lab_dphidz = 1./absdphi - dphidz * dphidz / CUBE(absdphi);

						double triplesum = SQR(nx_cryst) * SQR(ny_cryst)
						                 + SQR(ny_cryst) * SQR(nz_cryst)
						                 + SQR(nz_cryst) * SQR(nx_cryst);
						double EPS_loc   = _EPS[n][n]
						                 + _EPS1[n][n] * triplesum
						                 + _EPS2[n][n] * SQR(nx_cryst * ny_cryst * nz_cryst)
						                 + _EPS3[n][n] * SQR(triplesum);
						(*s_EPS_loc[n])(k, j, i) = EPS_loc;

						local_F += 0.5 * SQR(EPS_loc) * dphi_SQR;
						/* dFdPHI += -d/dx ( SQR(EPS_loc) * dphidx + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidx) )
						             -d/dy ( SQR(EPS_loc) * dphidy + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidy) )
						             -d/dz ( SQR(EPS_loc) * dphidz + dphi_SQR * EPS_loc * d(EPS_loc)/d(dphidz) ); */

						double dEPS_dnx_cryst = 2.0 * _EPS1[n][n] * nx_cryst * (SQR(ny_cryst) + SQR(nz_cryst))
						                      + 2.0 * _EPS2[n][n] * nx_cryst *  SQR(ny_cryst) * SQR(nz_cryst)
						                      + 4.0 * _EPS3[n][n] * nx_cryst * (SQR(ny_cryst) + SQR(nz_cryst)) * triplesum;
						double dEPS_dny_cryst = 2.0 * _EPS1[n][n] * ny_cryst * (SQR(nz_cryst) + SQR(nx_cryst))
						                      + 2.0 * _EPS2[n][n] * ny_cryst *  SQR(nz_cryst) * SQR(nx_cryst)
						                      + 4.0 * _EPS3[n][n] * ny_cryst * (SQR(nz_cryst) + SQR(nx_cryst)) * triplesum;
						double dEPS_dnz_cryst = 2.0 * _EPS1[n][n] * nz_cryst * (SQR(nx_cryst) + SQR(ny_cryst))
						                      + 2.0 * _EPS2[n][n] * nz_cryst *  SQR(nx_cryst) * SQR(ny_cryst)
						                      + 4.0 * _EPS3[n][n] * nz_cryst * (SQR(nx_cryst) + SQR(ny_cryst)) * triplesum;
  
						double dEPS_dnx_lab = _R[0][0] * dEPS_dnx_cryst
						                    + _R[1][0] * dEPS_dny_cryst
						                    + _R[2][0] * dEPS_dnz_cryst;
						double dEPS_dny_lab = _R[0][1] * dEPS_dnx_cryst
						                    + _R[1][1] * dEPS_dny_cryst
						                    + _R[2][1] * dEPS_dnz_cryst;
						double dEPS_dnz_lab = _R[0][2] * dEPS_dnx_cryst
						                    + _R[1][2] * dEPS_dny_cryst
						                    + _R[2][2] * dEPS_dnz_cryst;

						//tmp_x[n](k,j,i) = SQR(EPS_loc)*dphidx + dphi_SQR*EPS_loc*(dEPS_dnx_lab*dnx_lab_dphidx + dEPS_dny_lab*dny_lab_dphidx + dEPS_dnz_lab*dnz_lab_dphidx);
						//tmp_y[n][ind] = SQR(EPS_loc)*dphidy + dphi_SQR*EPS_loc*(dEPS_dnx_lab*dnx_lab_dphidy + dEPS_dny_lab*dny_lab_dphidy + dEPS_dnz_lab*dnz_lab_dphidy);
						//tmp_z[n][ind] = SQR(EPS_loc)*dphidz + dphi_SQR*EPS_loc*dEPS_dphidz;

						/* new scheme for computing variational derivative, 2012/04/27  */
						(*s_tmp_x[n])(k, j, i) = dphi_SQR * EPS_loc * (dEPS_dnx_lab * dnx_lab_dphidx + dEPS_dny_lab * dny_lab_dphidx + dEPS_dnz_lab * dnz_lab_dphidx);
						(*s_tmp_y[n])(k, j, i) = dphi_SQR * EPS_loc * (dEPS_dnx_lab * dnx_lab_dphidy + dEPS_dny_lab * dny_lab_dphidy + dEPS_dnz_lab * dnz_lab_dphidy);
						(*s_tmp_z[n])(k, j, i) = dphi_SQR * EPS_loc * (dEPS_dnx_lab * dnx_lab_dphidz + dEPS_dny_lab * dny_lab_dphidz + dEPS_dnz_lab * dnz_lab_dphidz);

						(*s_dFdPHI[n])(k, j, i) += - 1.0 * SQR(EPS_loc) * d2phi;
					}
				}
			}
		}

#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &local_F);
#endif

		_F += local_F;

		/* cubic anisotropic interface energy */
		if ((_EPS1[n][n] != 0) || (_EPS2[n][n] != 0) || (_EPS3[n][n] != 0))
		{
			(s_tmp_x[n])->set_boundary_periodic();
			(s_tmp_y[n])->set_boundary_periodic();
			(s_tmp_z[n])->set_boundary_periodic();
			(s_EPS_loc[n])->set_boundary_periodic();

#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_dFdPHI[n])(k, j, i) += - ((*s_tmp_x[n])(k,     j,     i + 1) - (*s_tmp_x[n])(k,     j,     i - 1)) / (2.0 * gridsize)
						                           - ((*s_tmp_y[n])(k,     j + 1, i    ) - (*s_tmp_y[n])(k,     j - 1, i    )) / (2.0 * gridsize)
						                           - ((*s_tmp_z[n])(k + 1, j,     i    ) - (*s_tmp_z[n])(k - 1, j,     i    )) / (2.0 * gridsize);

						/* new scheme for computing variational derivative, 2012/04/27 */
						(*s_dFdPHI[n])(k, j, i) += - (SQR((*s_EPS_loc[n])(k,     j,     i + 1)) - SQR((*s_EPS_loc[n])(k,     j,     i - 1))) / (2.0 * gridsize) * (*s_dPHIdx[n])(k, j, i)
						                           - (SQR((*s_EPS_loc[n])(k,     j + 1, i    )) - SQR((*s_EPS_loc[n])(k,     j - 1, i    ))) / (2.0 * gridsize) * (*s_dPHIdy[n])(k ,j, i)
						                           - (SQR((*s_EPS_loc[n])(k + 1, j,     i    )) - SQR((*s_EPS_loc[n])(k - 1, j,     i    ))) / (2.0 * gridsize) * (*s_dPHIdz[n])(k, j, i);
					}
				}
			}
		}

		_F *= CUBE(gridsize);
                _Fraw = _F;

		/* Cahn-Hillard equation */
		if (dynamics_type == 1)
		{
			(s_dFdPHI[n])->set_boundary_periodic();

#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_d2dFdPHI[n])(k, j, i) = ((*s_dFdPHI[n])(k,     j,     i + 1) + (*s_dFdPHI[n])(k,     j,     i - 1)
						                          +  (*s_dFdPHI[n])(k,     j + 1, i    ) + (*s_dFdPHI[n])(k,     j - 1, i    )
						                          +  (*s_dFdPHI[n])(k + 1, j,     i    ) + (*s_dFdPHI[n])(k - 1, j,     i    )
						                          -  (*s_dFdPHI[n])(k,     j,     i    ) * 6.0) / h2;
					}
				}
			}
		}

		/* Cahn-Hillard equation */
		if ((dynamics_type == 0) || (dynamics_type == 2))
		{
#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_dPHIdt0[n])(k, j, i) = -1.0 * Mob_GL * (*s_dFdPHI[n])(k, j, i);
					}
				}
			}
		}
		/* Ginzburg-Landau equation */
		else if (dynamics_type == 1)
		{
#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_dPHIdt0[n])(k, j, i) = -1.0 * Mob_D * (*s_d2dFdPHI[n])(k, j, i);
					}
				}
			}
		}

		/* no additional constraints */
		if ((dynamics_type == 0) || (dynamics_type == 1))
		{
			/* 0: Ginzburg-Landau equation: not conserved */
			/* 1: Cahn-Hilliard equation: conserved */
#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_dPHIdt[n])(k, j, i) = (*s_dPHIdt0[n])(k, j, i);
					}
				}
			}
		}
		/* Ginzburg-Landau equation but with constant volume constraint */
		else if (dynamics_type == 2) 
		{
			double avg_dPHIdt0 = 0;

// FIXME OpenMP makes numerical diff
//#pragma omp parallel for reduction(+:avg_dPHIdt0)
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						avg_dPHIdt0 += (*s_dPHIdt0[n])(k, j, i);
					}
				}
			}

#ifdef _STK_MPI
			_node->reduce(STK_SUM, STK_DOUBLE, &avg_dPHIdt0);
#endif

			avg_dPHIdt0 /= (_NX * _NY * _NZ);

#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_dPHIdt[n])(k, j, i) = (*s_dPHIdt0[n])(k, j, i) - avg_dPHIdt0;
					}
				}
			}
		}
		else
		{
			ERROR("unknown dynamics_type = " << dynamics_type);
		}
	} /* end of for(n=0;n<num_fields;n++) */

	/* calculate the max value of dphidt */
	double Gmiddle = 0;
	for (int n = 0; n < num_fields; n++)
	{
#pragma omp parallel for reduction(max:Gmiddle)
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					double temp_max = fabs((*s_dPHIdt[n])(k, j, i));
					if (Gmiddle <= temp_max)
					{
						Gmiddle = temp_max;
					}
				}
			}
		}
	}

#ifdef _STK_MPI
	_node->reduce(STK_MAX, STK_DOUBLE, &Gmiddle);
#endif

	_G = Gmiddle;
    
    /* adjust time step at every step */
    timestep = dtmax/fmax(_G*dtmax/dphimax,1);
    
	return _F;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

/* multi-phase field model (MPF) */
double PhaseFieldFrame::FreeEnergy_multi_mpi()
{
	double h2 = gridsize * gridsize;

	_F = 0;
    _Fraw = 0; 
        
	if (gridsize <= 0)
	{
		FATAL("FreeEnergy: gridesize = " << gridsize);
	}

	// XXX for _NZ
	const int z_begin = _idx->get_range_x().begin();
	const int z_end   = _idx->get_range_x().end();

	// XXX for _NY
	const int y_begin = _idx->get_range_y().begin();
	const int y_end   = _idx->get_range_y().end();

	// XXX for _NX
	const int x_begin = _idx->get_range_z().begin();
	const int x_end   = _idx->get_range_z().end();

	for (int n = 0; n < num_fields; n++)
	{
		(s_PHI[n])->set_boundary_periodic();
	}

	/* Zero out time derivative of phase fields */
	for (int n = 0; n < num_fields; n++)
	{
#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					(*s_dPHIdt0[n])(k, j, i) = 0;
					(*s_dFdPHI[n])(k, j, i) = 0;
					(*s_dPHIdtCH[n])(k, j, i) = 0;
				}
			}
		}
	}

	/* Reset liquid chemical potential, will be updated by constraint */
	if (dynamics_type == 2 || dynamics_type == 3 || dynamics_type == 8)
	{
		_M[0][1] = 0;
		_M[0][2] = _M[1][2];
		_M[1][0] = 0;
		_M[2][0] = _M[2][1];
		_MU[0]   = _MU[1];
	}
     if (dynamics_type == 4 || dynamics_type == 5 || dynamics_type == 6 || dynamics_type == 7)
    {
        _MU[0] = _MU[1];
        _MU[2] = _MU[1];
    }
    
    
	/* compute spatial gradients */
	for(int n = 0; n < num_fields; n++)
	{
#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					(*s_d2PHI[n])(k, j, i) = ((*s_PHI[n])(k,     j,     i + 1) + (*s_PHI[n])(k,     j,     i - 1)
					                       +  (*s_PHI[n])(k,     j + 1, i    ) + (*s_PHI[n])(k,     j - 1, i    )
					                       +  (*s_PHI[n])(k + 1, j,     i    ) + (*s_PHI[n])(k - 1, j,     i    )
					                       -  (*s_PHI[n])(k,     j,     i    ) * 6.0) / h2;

					(*s_dPHIdx[n])(k, j, i) = ((*s_PHI[n])(k,     j,     i + 1) - (*s_PHI[n])(k,     j,     i - 1)) / (2.0 * gridsize);
					(*s_dPHIdy[n])(k, j, i) = ((*s_PHI[n])(k,     j + 1, i    ) - (*s_PHI[n])(k,     j - 1, i    )) / (2.0 * gridsize);
					(*s_dPHIdz[n])(k, j, i) = ((*s_PHI[n])(k + 1, j,     i    ) - (*s_PHI[n])(k - 1, j,     i    )) / (2.0 * gridsize);
				}
			}
		}
	}

	/* Compute Free energy and time derivative of phase fields */

	/* assign the interface mobility and assigan the gradient matrix for external force term*/

#pragma omp parallel for
	for (int i = x_begin; i < x_end; i++)
	{
		for (int j = y_begin; j < y_end; j++)
		{
			for (int k = z_begin; k < z_end; k++)
			{
				(*s_K_D)(k, j, i) = Mob_D;
			}
		}
	}

#pragma omp parallel for
	for (int i = x_begin; i < x_end; i++)
	{
		for (int j = y_begin; j < y_end; j++)
		{
			for (int k = z_begin; k < z_end; k++)
			{
				/* define the s-v interface where it's near the solid boundary
				   but away from the liquid boundary */
				double phi_1 = (*s_PHI[1])(k,     j,     i    ) - 0.5;
				if (((phi_1 * ((*s_PHI[1])(k,     j,     i + 1) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j,     i - 1) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j,     i + 2) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j,     i - 2) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j,     i + 3) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j,     i - 3) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j + 1, i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j - 1, i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j + 2, i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j - 2, i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j + 3, i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k,     j - 3, i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k + 1, j,     i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k - 1, j,     i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k + 2, j,     i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k - 2, j,     i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k + 3, j,     i    ) - 0.5) <= 0)
				  || (phi_1 * ((*s_PHI[1])(k - 3, j,     i    ) - 0.5) <= 0)))
				{
					(*s_K_S)(k, j, i) = 1;
				}
				else
				{
					(*s_K_S)(k, j, i) = 0;
				}

				double phi_0 = (*s_PHI[0])(k,     j,     i    ) - 0.5;
				if (((phi_0 * ((*s_PHI[0])(k,     j,     i + 1) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j,     i - 1) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j,     i + 2) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j,     i - 2) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j,     i + 3) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j,     i - 3) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j + 1, i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j - 1, i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j + 2, i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j - 2, i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j + 3, i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k,     j - 3, i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k + 1, j,     i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k - 1, j,     i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k + 2, j,     i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k - 2, j,     i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k + 3, j,     i    ) - 0.5) <= 0)
				  || (phi_0 * ((*s_PHI[0])(k - 3, j,     i    ) - 0.5) <= 0)))
				{ 
					(*s_K_L)(k, j, i) = 1;
				}
				else
				{
					(*s_K_L)(k, j, i) = 0;
				}

				if (((*s_K_L)(k, j, i) == 1) && ((*s_K_S)(k, j, i) == 1))
				{
					(*s_K_LS)(k, j, i) = 1;
				}
				else
				{
					(*s_K_LS)(k, j, i) = 0;
				}
				if (((*s_K_L)(k, j, i) == 0) && ((*s_K_S)(k, j, i) == 1))
				{
					(*s_K_SV)(k, j, i) = 1;
				}
				else
				{
					(*s_K_SV)(k, j, i) = 0;
				}

				if (((*s_K_L)(k, j, i) == 1) && ((*s_K_S)(k, j, i) == 0))
				{
					(*s_K_LV)(k, j, i) = 1;
				}
				else
				{
					(*s_K_LV)(k, j, i) = 0;
				}

				(*s_K_other)(k, j, i) = 1 - (*s_K_LV)(k, j, i) - (*s_K_SV)(k, j, i) - (*s_K_LS)(k, j, i);
			}
		}
	}

	for(int n = 0; n < num_fields; n++)
	{
		double local_F = 0.0;
                double local_Fraw = 0.0;
#pragma omp parallel for reduction(+:local_F)
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					double grad_term = SQR((*s_dPHIdx[n])(k, j, i)) + SQR((*s_dPHIdy[n])(k, j, i)) + SQR((*s_dPHIdz[n])(k, j, i));

					double phi_0 = (*s_PHI[0])(k, j, i);
					double phi_n = (*s_PHI[n])(k, j, i);

					double f = SQR(phi_n) * SQR(1 - phi_n);
					//double fp = phi_n;
					double fp = tanh((phi_n - 0.5) / 0.10) / 2 + 0.5;
					local_F += f * _U[n][n] + fp * _MU[n];
                                        local_Fraw += f * _U[n][n];
					// the influence of the external force (only for liquid)
					if ( n == 0 )
					{
						local_F += Fext_x * fp * (*s_matrix_x)(k, j, i) + Fext_y * fp * (*s_matrix_y)(k, j, i);
					}

					// Penalty term for three phase co-existing    
					local_F += Penalty_3phase / num_fields * SQR(phi_0) * SQR((*s_PHI[1])(k, j, i)) * SQR((*s_PHI[2])(k, j, i));

                    // Penalty term for deviation from original NW shape
                    local_F += Penalty_NW/num_fields*SQR((*s_PHI[1])(k, j, i)-(*s_NW_orig)(k,j,i));
                    
					/* isotropic interface energy */
					if (((_EPS1[n][n] == 0) && (_EPS2[n][n] == 0) && (_EPS3[n][n] == 0)))
					{
						double EPS_loc = 1.0;
						double EPS_prefactor = SQR(EPS_loc);

						local_F += _EPS[n][n] * EPS_prefactor * grad_term;
                                                local_Fraw += _EPS[n][n] * EPS_prefactor * grad_term;
						//local_F += f * _U[n][n] + fp * _MU[n];

						/* the equation of motion should follow the Steinbach 1999 formulation */
						(*s_dFdPHI[n])(k, j, i) = _U[n][n] * (4.0 * phi_n * phi_n * phi_n + 2.0 * phi_n - 6.0 * phi_n * phi_n)
						                        + _MU[n]   * (5.0 - 5.0 * SQR(tanh(10 * phi_n - 5.0)))
						                        - 2.0 * _EPS[n][n] * EPS_prefactor * (*s_d2PHI[n])(k, j, i);

						/* _dFdPHIi[n][ind] = _U[n][n] * (4.0 * _PHI[n][ind] * _PHI[n][ind] * _PHI[n][ind]
						  + 2.0 * _PHI[n][ind] - 6.0 * _PHI[n][ind] * _PHI[n][ind]) 
						  + _MU[n] - 2.0 * _EPS[n][n] * EPS_prefactor * _d2PHI[n][ind]; */

						// with the external force 
						if (n == 0)
						{
							(*s_dFdPHI[0])(k, j, i) += Fext_x * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0))) * (*s_matrix_x)(k, j, i)
							                         + Fext_y * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0))) * (*s_matrix_y)(k, j, i);
						}
					}
					/* anisotropic interface energy */
					else
					{
						//if (  !( ((n==1)&&(p==0)) || ((n==0)&&(p==1)) || ((n==1)&&(p==2)) || ((n==2)&&(p==1)) ))
						//{
						//   FATAL("FreeEnergy_multi: something is wrong");
						//}
						//if (  !( ((n==1)&&(p==2)) || ((n==2)&&(p==1)) ) )
						//{
						//   FATAL("FreeEnergy_multi: something is wrong");
						//}

						/* the EPS term will be computed and added below */

						/* surface normal defined by phase gradient for all the three phase */
						double dphidx = (*s_dPHIdx[n])(k, j, i);
						double dphidy = (*s_dPHIdy[n])(k, j, i);
						double dphidz = (*s_dPHIdz[n])(k, j, i);

						double absdphi = sqrt((SQR(dphidx) + SQR(dphidy) + SQR(dphidz)) + 1e-24);

						double nx_lab = dphidx / absdphi;
						double ny_lab = dphidy / absdphi;
						double nz_lab = dphidz / absdphi;

						double nx_cryst = _R[0][0] * nx_lab + _R[0][1] * ny_lab + _R[0][2] * nz_lab;
						double ny_cryst = _R[1][0] * nx_lab + _R[1][1] * ny_lab + _R[1][2] * nz_lab;
						double nz_cryst = _R[2][0] * nx_lab + _R[2][1] * ny_lab + _R[2][2] * nz_lab;

						double eps1 = _EPS1[n][n]; 
						double eps2 = _EPS2[n][n]; 
						double eps3 = _EPS3[n][n];

						double triplesum = SQR(nx_cryst) * SQR(ny_cryst) + SQR(ny_cryst) * SQR(nz_cryst) + SQR(nz_cryst) * SQR(nx_cryst);

						/* Notice the difference from the Matlab implementation */
						double EPS_loc   = 1.0 + eps1 * triplesum + eps2 * SQR(nx_cryst * ny_cryst * nz_cryst) + eps3 * SQR(triplesum);
						(*s_EPS_loc[n])(k, j, i) = EPS_loc;

						/* Free energy for anisotropic surface energy */
						//f = SQR(_PHI[n][ind])*SQR(1-_PHI[n][ind]);
						//fp = _PHI[n][ind];
						//fp = tanh((_PHI[n][ind]-0.5)/0.10)/2 + 0.5;
						local_F += _EPS[n][n] * SQR(EPS_loc) * grad_term;;
                                                local_Fraw += _EPS[n][n] * SQR(EPS_loc) * grad_term;;
						//local_F += f * _U[n][n] + fp * _MU[n] + _EPS[n][n] * SQR(EPS_loc) * grad_term;

						double dnx_lab_dphidx = 1. / absdphi - dphidx * dphidx / CUBE(absdphi);
						double dny_lab_dphidx =              - dphidy * dphidx / CUBE(absdphi);
						double dnz_lab_dphidx =              - dphidz * dphidx / CUBE(absdphi);

						double dnx_lab_dphidy =              - dphidx * dphidy / CUBE(absdphi);
						double dny_lab_dphidy = 1. / absdphi - dphidy * dphidy / CUBE(absdphi);
						double dnz_lab_dphidy =              - dphidz * dphidy / CUBE(absdphi);

						double dnx_lab_dphidz =              - dphidx * dphidz / CUBE(absdphi);
						double dny_lab_dphidz =              - dphidy * dphidz / CUBE(absdphi);
						double dnz_lab_dphidz = 1. / absdphi - dphidz * dphidz / CUBE(absdphi);

						double dEPS_dnx_cryst = 2. * eps1 * nx_cryst * (SQR(ny_cryst) + SQR(nz_cryst))
						                      + 2. * eps2 * nx_cryst *  SQR(ny_cryst) * SQR(nz_cryst)
						                      + 4. * eps3 * nx_cryst * (SQR(ny_cryst) + SQR(nz_cryst)) * triplesum;
						double dEPS_dny_cryst = 2. * eps1 * ny_cryst * (SQR(nz_cryst) + SQR(nx_cryst))
						                      + 2. * eps2 * ny_cryst *  SQR(nz_cryst) * SQR(nx_cryst)
						                      + 4. * eps3 * ny_cryst * (SQR(nz_cryst) + SQR(nx_cryst)) * triplesum;
						double dEPS_dnz_cryst = 2. * eps1 * nz_cryst * (SQR(nx_cryst) + SQR(ny_cryst))
						                      + 2. * eps2 * nz_cryst *  SQR(nx_cryst) * SQR(ny_cryst)
						                      + 4. * eps3 * nz_cryst * (SQR(nx_cryst) + SQR(ny_cryst)) * triplesum;

						double dEPS_dnx_lab = _R[0][0] * dEPS_dnx_cryst + _R[1][0] * dEPS_dny_cryst + _R[2][0] * dEPS_dnz_cryst;
						double dEPS_dny_lab = _R[0][1] * dEPS_dnx_cryst + _R[1][1] * dEPS_dny_cryst + _R[2][1] * dEPS_dnz_cryst;
						double dEPS_dnz_lab = _R[0][2] * dEPS_dnx_cryst + _R[1][2] * dEPS_dny_cryst + _R[2][2] * dEPS_dnz_cryst;

						double dEPS_dphidx_loc, dEPS_dphidy_loc, dEPS_dphidz_loc;
					
							dEPS_dphidx_loc = dEPS_dnx_lab * dnx_lab_dphidx + dEPS_dny_lab * dny_lab_dphidx + dEPS_dnz_lab * dnz_lab_dphidx;
							dEPS_dphidy_loc = dEPS_dnx_lab * dnx_lab_dphidy + dEPS_dny_lab * dny_lab_dphidy + dEPS_dnz_lab * dnz_lab_dphidy;
							dEPS_dphidz_loc = dEPS_dnx_lab * dnx_lab_dphidz + dEPS_dny_lab * dny_lab_dphidz + dEPS_dnz_lab * dnz_lab_dphidz;
    
                        double grad_loc_SQR = SQR((*s_dPHIdx[n])(k, j, i)) + SQR((*s_dPHIdy[n])(k, j, i)) + SQR((*s_dPHIdz[n])(k, j, i));

						//tmp_x[n][ind] = SQR(EPS_loc)*_dPHIdx[n][ind] + grad_loc_SQR*EPS_loc*dEPS_dphidx_loc;
						//tmp_y[n][ind] = SQR(EPS_loc)*_dPHIdy[n][ind] + grad_loc_SQR*EPS_loc*dEPS_dphidy_loc;
						//tmp_z[n][ind] = SQR(EPS_loc)*_dPHIdz[n][ind] + grad_loc_SQR*EPS_loc*dEPS_dphidz_loc;

						/* new scheme for computing variational derivative, 2012/04/27  */
						(*s_tmp_x[n])(k, j, i) = grad_loc_SQR * EPS_loc * dEPS_dphidx_loc;
						(*s_tmp_y[n])(k, j, i) = grad_loc_SQR * EPS_loc * dEPS_dphidy_loc;
						(*s_tmp_z[n])(k, j, i) = grad_loc_SQR * EPS_loc * dEPS_dphidz_loc;

						(*s_dFdPHI[n])(k, j, i) += _U[n][n] * (4 * phi_n * phi_n * phi_n + 2 * phi_n - 6 * phi_n * phi_n)
						                          + _MU[n]   * (5.0 - 5.0 * SQR(tanh(10 * phi_n - 5.0)));

						/* _dFdPHIi[n][ind] += _U[n][n]*(4*_PHI[n][ind]*_PHI[n][ind]*_PHI[n][ind]
						                     + 2*_PHI[n][ind] - 6*_PHI[n][ind]*_PHI[n][ind])
						                     + _MU[n]; */

						(*s_dFdPHI[n])(k, j, i) += (-2.0) * _EPS[n][n] * SQR(EPS_loc) * (*s_d2PHI[n])(k, j, i);

						// with the external force term
						if (n == 0)
						{
							(*s_dFdPHI[0])(k, j, i) += Fext_x * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0))) * (*s_matrix_x)(k, j, i)
							                         + Fext_y * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0))) * (*s_matrix_y)(k, j, i);
						}
					}
				}
			}
		}

#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &local_F);
        _node->reduce(STK_SUM, STK_DOUBLE, &local_Fraw);
#endif

 		_F += local_F;
        _Fraw += local_Fraw;               

		/* cubic anisotropic interface energy */
		if ((_EPS1[n][n] != 0) || (_EPS2[n][n] != 0) || (_EPS3[n][n] != 0))
		{
			(s_tmp_x[n])->set_boundary_periodic();
			(s_tmp_y[n])->set_boundary_periodic();
			(s_tmp_z[n])->set_boundary_periodic();
			(s_EPS_loc[n])->set_boundary_periodic();

#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_dFdPHI[n])(k, j, i) += _EPS[n][n] * (-2.0) * ((*s_tmp_x[n])(k,     j,     i + 1) - (*s_tmp_x[n])(k,     j,     i - 1)) / (2.0 * gridsize)
						                         + _EPS[n][n] * (-2.0) * ((*s_tmp_y[n])(k,     j + 1, i    ) - (*s_tmp_y[n])(k,     j - 1, i    )) / (2.0 * gridsize)
						                         + _EPS[n][n] * (-2.0) * ((*s_tmp_z[n])(k + 1, j,     i    ) - (*s_tmp_z[n])(k - 1, j,     i    )) / (2.0 * gridsize);

						/* new scheme for computing variational derivative, 2012/04/27 */
						double eps_temp = _EPS[n][n] * (-4.0) * (*s_EPS_loc[n])(k, j, i);
						(*s_dFdPHI[n])(k, j, i) += eps_temp * ((*s_EPS_loc[n])(k,     j,     i + 1) - (*s_EPS_loc[n])(k,     j,     i - 1)) / (2.0 * gridsize) * (*s_dPHIdx[n])(k, j, i)
						                         + eps_temp * ((*s_EPS_loc[n])(k,     j + 1, i    ) - (*s_EPS_loc[n])(k,     j - 1, i    )) / (2.0 * gridsize) * (*s_dPHIdy[n])(k, j, i)
						                         + eps_temp * ((*s_EPS_loc[n])(k + 1, j,     i    ) - (*s_EPS_loc[n])(k - 1, j,     i    )) / (2.0 * gridsize) * (*s_dPHIdz[n])(k, j, i);
					}
				}
			}
		}
	}

	// variational derivative due to the triple junction penalty term
#pragma omp parallel for
	for (int i = x_begin; i < x_end; i++)
	{
		for (int j = y_begin; j < y_end; j++)
		{
			for (int k = z_begin; k < z_end; k++)
			{
				(*s_dFdPHI[0])(k, j, i) += 2 * Penalty_3phase * (*s_PHI[0])(k, j, i) * SQR((*s_PHI[1])(k, j, i)) * SQR((*s_PHI[2])(k, j, i));
				(*s_dFdPHI[1])(k, j, i) += 2 * Penalty_3phase * (*s_PHI[1])(k, j, i) * SQR((*s_PHI[0])(k, j, i)) * SQR((*s_PHI[2])(k, j, i));
				(*s_dFdPHI[2])(k, j, i) += 2 * Penalty_3phase * (*s_PHI[2])(k, j, i) * SQR((*s_PHI[0])(k, j, i)) * SQR((*s_PHI[1])(k, j, i));
                
                /* variational derivative due to the NW shape penalty term */
                (*s_dFdPHI[1])(k, j, i) += 2*Penalty_NW*((*s_PHI[1])(k, j, i) - (*s_NW_orig)(k, j, i));
                
			}
		}
	}

	// equation of motion
#pragma omp parallel for
	for (int i = x_begin; i < x_end; i++)
	{
		for (int j = y_begin; j < y_end; j++)
		{
			for (int k = z_begin; k < z_end; k++)
			{
				double dfdphi_0 = (*s_dFdPHI[0])(k, j, i);
				double dfdphi_1 = (*s_dFdPHI[1])(k, j, i);
				double dfdphi_2 = (*s_dFdPHI[2])(k, j, i);

				/* Steinback 1999 formulation */
				(*s_dPHIdt0[0])(k, j, i) = - (*s_K_SV)(k, j, i) * (Mob_GL *  (dfdphi_0 - dfdphi_1) + Mob_GL * (dfdphi_0 - dfdphi_2))
				                           - (*s_K_LV)(k, j, i) * (Mob_GL * ((dfdphi_0 + dfdphi_2) / 2 - dfdphi_1) + Mob_LV * (dfdphi_0 - dfdphi_2))
				                           - (*s_K_LS)(k, j, i) * (Mob_GL * ((dfdphi_0 + dfdphi_1) / 2 - dfdphi_2) + Mob_LS * (dfdphi_0 - dfdphi_1))
				                           - (*s_K_other)(k, j, i) * (Mob_GL * (dfdphi_0 - dfdphi_1) + Mob_GL * (dfdphi_0 - dfdphi_2));

				(*s_dPHIdt0[1])(k, j, i) = - (*s_K_LV)(k, j, i) * (Mob_GL *  (dfdphi_1 - dfdphi_2) + Mob_GL * (dfdphi_1 - dfdphi_0))
				                           - (*s_K_SV)(k, j, i) * (Mob_GL * ((dfdphi_1 + dfdphi_2) / 2 - dfdphi_0) + Mob_SV * (dfdphi_1 - dfdphi_2))
				                           - (*s_K_LS)(k, j, i) * (Mob_GL * ((dfdphi_0 + dfdphi_1) / 2 - dfdphi_2) + Mob_LS * (dfdphi_1 - dfdphi_0))
				                           - (*s_K_other)(k, j, i) * (Mob_GL * (dfdphi_1 - dfdphi_0) + Mob_GL * (dfdphi_1 - dfdphi_2));

				(*s_dPHIdt0[2])(k, j, i) = - (*s_K_LS)(k, j, i) * (Mob_GL *  (dfdphi_2 - dfdphi_1) + Mob_GL * (dfdphi_2 - dfdphi_0))
				                           - (*s_K_LV)(k, j, i) * (Mob_GL * ((dfdphi_0 + dfdphi_2) / 2 - dfdphi_1) + Mob_LV * (dfdphi_2 - dfdphi_0))
				                           - (*s_K_SV)(k, j, i) * (Mob_GL * ((dfdphi_1 + dfdphi_2) / 2 - dfdphi_0) + Mob_SV * (dfdphi_2 - dfdphi_1))
				                           - (*s_K_other)(k, j, i) * (Mob_GL * (dfdphi_2 - dfdphi_0) + Mob_GL * (dfdphi_2 - dfdphi_1));

				(*s_dPHIdtCH[0])(k, j, i) = -(*s_K_D)(k, j, i) * (dfdphi_0 - dfdphi_1) - (*s_K_D)(k, j, i) * (dfdphi_0 - dfdphi_2);
				(*s_dPHIdtCH[1])(k, j, i) = -(*s_K_D)(k, j, i) * (dfdphi_1 - dfdphi_0) - (*s_K_D)(k, j, i) * (dfdphi_1 - dfdphi_2);
				(*s_dPHIdtCH[2])(k, j, i) = -(*s_K_D)(k, j, i) * (dfdphi_2 - dfdphi_0) - (*s_K_D)(k, j, i) * (dfdphi_2 - dfdphi_1);
			}
		}
	}

	_F *= CUBE(gridsize);
    _Fraw *= CUBE(gridsize);
    
	/* calculate the max value of dphidt */
	double Gmiddle = 0;
    if (dynamics_type != 1)
    {
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {
                    double temp_max = fmax(fmax(fabs((*s_dPHIdt0[0])(k, j, i)),
                                                fabs((*s_dPHIdt0[1])(k, j, i))),
                                           fabs((*s_dPHIdt0[2])(k, j, i)));
                    if (Gmiddle <= temp_max)
                    {
                        Gmiddle = temp_max;
                    }
                }
            }
        }  
    }
    else if (dynamics_type == 1)
    {
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {
                    double temp_max = fmax(fmax(fabs((*s_dPHIdtCH[0])(k, j, i)),
                                                fabs((*s_dPHIdtCH[1])(k, j, i))),
                                           fabs((*s_dPHIdtCH[2])(k, j, i)));
                    if (Gmiddle <= temp_max)
                    {
                        Gmiddle = temp_max;
                    }
                }
            }
        }  
    }
    
#ifdef _STK_MPI
	_node->reduce(STK_MAX, STK_DOUBLE, &Gmiddle);
#endif
    
	_G = Gmiddle;

    /* adjust time step at every step */
    timestep = dtmax/fmax(_G*dtmax/dphimax,1);
    
    double vol_incr = vol_incr0/timestep;
    double x_incr = x_incr0*liquid_volume/timestep;
    double y_incr = y_incr0*liquid_volume/timestep;
    
    if (dynamics_type != 2 && dynamics_type != 3 && dynamics_type != 4 && dynamics_type != 5 && dynamics_type != 6 && dynamics_type != 7 && dynamics_type != 8) 
    {
		for (int n = 0; n < num_fields; n++)
		{
			/* Cahn-Hillard equation */
			if (dynamics_type == 1)
			{
				(s_dPHIdtCH[n])->set_boundary_periodic();

#pragma omp parallel for
				for (int i = x_begin; i < x_end; i++)
				{
					for (int j = y_begin; j < y_end; j++)
					{
						for (int k = z_begin; k < z_end; k++)
						{
							(*s_d2dPHIdt0[n])(k, j, i) = ((*s_dPHIdtCH[n])(k,     j,     i + 1) + (*s_dPHIdtCH[n])(k,     j,     i - 1)
							                           +  (*s_dPHIdtCH[n])(k,     j + 1, i    ) + (*s_dPHIdtCH[n])(k,     j - 1, i    )
							                           +  (*s_dPHIdtCH[n])(k + 1, j,     i    ) + (*s_dPHIdtCH[n])(k - 1, j,     i    )
							                           -  (*s_dPHIdtCH[n])(k,     j,     i    ) * 6.0 ) / h2;
						}
					}
				}
			}

			/* no additional constraints */
			if (dynamics_type == 0)
			{
				/* 0: Ginzburg-Landau equation: not conserved */
#pragma omp parallel for
				for (int i = x_begin; i < x_end; i++)
					for (int j = y_begin; j < y_end; j++)
						for (int k = z_begin; k < z_end; k++)
							(*s_dPHIdt[n])(k, j, i) = (*s_dPHIdt0[n])(k, j, i);
			}
			/* no additional constraints */
			else if (dynamics_type == 1)
			{
				/* 1: Cahn-Hilliard equation: conserved */
#pragma omp parallel for
				for (int i = x_begin; i < x_end; i++)
					for (int j = y_begin; j < y_end; j++)
						for (int k = z_begin; k < z_end; k++)
							(*s_dPHIdt[n])(k, j, i) = -(*s_d2dPHIdt0[n])(k, j, i);
			}
		}
	}
	else if (dynamics_type == 2)
	{
		/* specify initial m value for the test case (need to be done before F is calculated) */
		//_M[0][0] = 0; _M[1][1] = 0; _M[2][2] = 0;
		//_M[0][1] = 0; _M[1][0] = 0; _M[0][2] = _M[1][2]; _M[2][0] = _M[2][1];

		/* Update _mprime */
#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					double phi_0 = (*s_PHI[0])(k, j, i);
					//_mprime[i*_NY*_NZ+j*_NZ+k] = -2*_K[i*_NY*_NZ+j*_NZ+k];
					(*s_mprime)(k, j, i) = - (*s_K_LS)(k, j, i) * (Mob_GL / 2 + Mob_LS) * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)))
					                       - (*s_K_SV)(k, j, i) * (Mob_GL     + Mob_GL) * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)))
					                       - (*s_K_LV)(k, j, i) * (Mob_GL / 2 + Mob_LV) * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)))
					                       - (*s_K_other)(k, j, i) * (Mob_GL + Mob_GL)  * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)));
					//_mprime[ind] = (-2.0)*_K[ind];
				}
			}
		}

		/* Initialize DM */	
		// _DM[0][0] = _DM[0][1] = _DM[0][2] = 0;
		// _DM[1][0] = _DM[1][1] = _DM[1][2] = 0;
		// _DM[2][0] = _DM[2][1] = _DM[2][2] = 0;

		/* Update DM */	
		double dmm = 0;
// FIXME OpenMP makes numerical diff
//#pragma omp parallel for reduction(+:dmm)
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
					dmm += SQR((*s_mprime)(k, j, i));

#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &dmm);
#endif

		double dmu1 = 0;
// FIXME OpenMP makes numerical diff
//#pragma omp parallel for reduction(+:dmu1)
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
					dmu1 -= ((*s_mprime)(k, j, i) * (*s_dPHIdt0[0])(k, j, i)) / dmm;

#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &dmu1);
#endif

		//_DM[1][0] = -_DM[0][1];
		//_DM[0][2] =  _DM[0][1];
		//_DM[2][0] = -_DM[0][1];

		/* Update gradient due to the change of m */
#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					(*s_dPHIdt[0])(k, j, i) = (*s_dPHIdt0[0])(k, j, i) + dmu1 * (*s_mprime)(k, j, i);
					(*s_dPHIdt[1])(k, j, i) = (*s_dPHIdt0[1])(k, j, i) + dmu1
					                        * ((*s_K_LS)(k, j, i) * (-Mob_GL / 2 + Mob_LS) + (*s_K_LV)(k, j, i) * Mob_GL + (*s_K_SV)(k, j, i) * Mob_GL + (*s_K_other)(k, j, i) * Mob_GL)
					                        * (5.0 - 5.0 * SQR(tanh(10.0 * (*s_PHI[0])(k, j, i) - 5.0)));
					(*s_dPHIdt[2])(k, j, i) = (*s_dPHIdt0[2])(k, j, i) + dmu1
					                        * ((*s_K_LS)(k, j, i) * Mob_GL + (*s_K_LV)(k, j, i) * (-Mob_GL / 2 + Mob_LV) + (*s_K_SV)(k, j, i) * Mob_GL + (*s_K_other)(k, j, i) * Mob_GL)
					                        * (5.0 - 5.0 * SQR(tanh(10.0 * (*s_PHI[0])(k, j, i) - 5.0)));

					/* _dPHIdt[1][i*_NY*_NZ+j*_NZ+k] = _dPHIdt0[1][i*_NY*_NZ+j*_NZ+k] + dmu1
					                                 * (_K[i*_NY*_NZ+j*_NZ+k]);
					   _dPHIdt[2][i*_NY*_NZ+j*_NZ+k] = _dPHIdt0[2][i*_NY*_NZ+j*_NZ+k] + dmu1
					                                 * (_K[i*_NY*_NZ+j*_NZ+k]); */
				}
			}
		}

		/* update chemical potential */
		_MU[0] = _MU[1] + dmu1;
#pragma omp parallel for
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];
	}
	/* G-L+C-H with liquid volume constriant */
	else if (dynamics_type == 3)
	{
#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					double phi_0 = (*s_PHI[0])(k, j, i);
					//_mprime[i*_NY*_NZ+j*_NZ+k] = -2*_K[i*_NY*_NZ+j*_NZ+k];
					(*s_mprime)(k, j, i) = - (*s_K_LS)(k, j, i) * (Mob_GL / 2 + Mob_LS) * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)))
					                       - (*s_K_SV)(k, j, i) * (Mob_GL     + Mob_GL) * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)))
					                       - (*s_K_LV)(k, j, i) * (Mob_GL / 2 + Mob_LV) * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)))
					                       - (*s_K_other)(k, j, i) * (Mob_GL + Mob_GL)  * (5.0 - 5.0 * SQR(tanh(10 * phi_0 - 5.0)));
					//_mprime[ind] = (-2.0)*_K[ind];
				}
			}
		}

		/* Initialize DM */	
		// _DM[0][0] = _DM[0][1] = _DM[0][2] = 0;
		// _DM[1][0] = _DM[1][1] = _DM[1][2] = 0;
		// _DM[2][0] = _DM[2][1] = _DM[2][2] = 0;

		/* Update DM */	
		double dmm = 0;
// FIXME OpenMP makes numerical diff
//#pragma omp parallel for reduction(+:dmm)
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
					dmm += SQR((*s_mprime)(k, j, i));

#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &dmm);
#endif

		double dmu1 = 0;
// FIXME OpenMP makes numerical diff
//#pragma omp parallel for reduction(+:dmu1)
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
					dmu1 -= ((*s_mprime)(k, j, i) * (*s_dPHIdt0[0])(k, j, i)) / dmm;

#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &dmu1);
#endif

		//_DM[1][0] = -_DM[0][1]; 
		//_DM[0][2] =  _DM[0][1]; 
		//_DM[2][0] = -_DM[0][1];
        
	        /* Update gradient due to the change of m */
#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					(*s_dPHIdt[0])(k, j, i) = (*s_dPHIdt0[0])(k, j, i) + dmu1 * (*s_mprime)(k, j, i);
					(*s_dPHIdt[1])(k, j, i) = (*s_dPHIdt0[1])(k, j, i) + dmu1
					                        * ((*s_K_LS)(k, j, i) * (-Mob_GL / 2 + Mob_LS) + (*s_K_LV)(k, j, i) * Mob_GL + (*s_K_SV)(k, j, i) * Mob_GL + (*s_K_other)(k, j, i) * Mob_GL)
					                        * (5.0 - 5.0 * SQR(tanh(10.0 * (*s_PHI[0])(k, j, i) - 5.0)));
					(*s_dPHIdt[2])(k, j, i) = (*s_dPHIdt0[2])(k, j, i) + dmu1
					                        * ((*s_K_LS)(k, j, i) * Mob_GL + (*s_K_LV)(k, j, i) * (-Mob_GL / 2 + Mob_LV) + (*s_K_SV)(k, j, i) * Mob_GL + (*s_K_other)(k, j, i) * Mob_GL)
					                        * (5.0 - 5.0 * SQR(tanh(10.0 * (*s_PHI[0])(k, j, i) - 5.0))); 

						/* _dPHIdt[1][i*_NY*_NZ+j*_NZ+k] = _dPHIdt0[1][i*_NY*_NZ+j*_NZ+k] + dmu1
						                                 * (_K[i*_NY*_NZ+j*_NZ+k]);
						   _dPHIdt[2][i*_NY*_NZ+j*_NZ+k] = _dPHIdt0[2][i*_NY*_NZ+j*_NZ+k] + dmu1
						                                 * (_K[i*_NY*_NZ+j*_NZ+k]); */
				}
			}
		}

		/* update chemical potential */
		_MU[0] = _MU[1] + dmu1;
#pragma omp parallel for
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];

		for (int n = 0; n < num_fields; n++)
		{
			(s_dPHIdtCH[n])->set_boundary_periodic();

#pragma omp parallel for
			for (int i = x_begin; i < x_end; i++)
			{
				for (int j = y_begin; j < y_end; j++)
				{
					for (int k = z_begin; k < z_end; k++)
					{
						(*s_d2dPHIdt0[n])(k, j, i) = ((*s_dPHIdtCH[n])(k,     j,     i + 1) + (*s_dPHIdtCH[n])(k,     j,     i - 1)
						                           +  (*s_dPHIdtCH[n])(k,     j + 1, i    ) + (*s_dPHIdtCH[n])(k,     j - 1, i    )
						                           +  (*s_dPHIdtCH[n])(k + 1, j,     i    ) + (*s_dPHIdtCH[n])(k - 1, j,     i    )
						                           -  (*s_dPHIdtCH[n])(k,     j,     i    ) * 6.0 ) / h2;
					}
				}
			}
		}

#pragma omp parallel for
		for (int i = x_begin; i < x_end; i++)
		{
			for (int j = y_begin; j < y_end; j++)
			{
				for (int k = z_begin; k < z_end; k++)
				{
					(*s_dPHIdt[0])(k, j, i) -= (*s_d2dPHIdt0[0])(k, j, i);
					(*s_dPHIdt[1])(k, j, i) -= (*s_d2dPHIdt0[1])(k, j, i);
					(*s_dPHIdt[2])(k, j, i) -= (*s_d2dPHIdt0[2])(k, j, i);
				}
			}
		}
	}
    
    else if (dynamics_type == 4)
    {   
        _MU[1] = 0;
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {   
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    (*s_dGLdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGLdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGVdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));       
                    
                    (*s_dGVdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGSdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGSdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                }
            }
        }  
        

        double A_11=0; double A_12=0; double A_21=0; double A_22=0; 
        double B_1=0; double B_2=0; 
    
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    A_11 += (*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_12 += (*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_21 += (*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_22 += (*s_dGVdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0))); 
                    B_1  += -(5.0-5.0*SQR(tanh(10*phi_0-5.0)))* (*s_dPHIdt0[0])(k,j,i);
                    B_2  += -(5.0-5.0*SQR(tanh(10*phi_2-5.0)))* (*s_dPHIdt0[2])(k,j,i);
                }    
#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &A_11);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_12);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_21);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_22);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_1);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_2);
#endif       
        
        B_2 += vol_incr;
        
        double dmuL = (A_12*B_2 - A_22*B_1)/(A_12*A_21-A_11*A_22);
        double dmuV = -(A_11*B_2 - A_21*B_1)/(A_12*A_21-A_11*A_22);
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {           
                    (*s_dPHIdt[0])(k,j,i) = (*s_dPHIdt0[0])(k,j,i) + dmuL*(*s_dGLdmuL)(k,j,i) + dmuV*(*s_dGLdmuV)(k,j,i);
            
                    (*s_dPHIdt[1])(k,j,i) = (*s_dPHIdt0[1])(k,j,i) + dmuL*(*s_dGSdmuL)(k,j,i) + dmuV*(*s_dGSdmuV)(k,j,i);
            
                    (*s_dPHIdt[2])(k,j,i) = (*s_dPHIdt0[2])(k,j,i) + dmuL*(*s_dGVdmuL)(k,j,i) + dmuV*(*s_dGVdmuV)(k,j,i);
            
                }
        
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        
        for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];

    }
    
    else if (dynamics_type == 5)
    {   
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {   
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    (*s_dGLdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGLdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGVdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));       
                    
                    (*s_dGVdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGSdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGSdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                }
            }
        }  
        
        
        double A_11=0; double A_12=0; double A_13=0; 
        double A_21=0; double A_22=0; double A_23=0;
        double A_31=0; double A_32=0; double A_33=0;
        double B_1=0; double B_2=0; double B_3=0;
        double dmuL=0; double dmuV=0; double dfY=0;
        double divider=0;
        _MU[1]=0;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    double mat_y = (*s_matrix_y)(k, j, i);
                    A_11 += (*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_12 += (*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_13 += mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_21 += (*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_22 += (*s_dGVdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0))); 
                    A_23 += mat_y*(*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_31 += mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_32 += mat_y*(*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_33 += SQR(mat_y)*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    B_1  += -(5.0-5.0*SQR(tanh(10*phi_0-5.0)))* (*s_dPHIdt0[0])(k,j,i);
                    B_2  += -(5.0-5.0*SQR(tanh(10*phi_2-5.0)))* (*s_dPHIdt0[2])(k,j,i);
                    B_3  += -mat_y*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))*(*s_dPHIdt0[0])(k,j,i);
                }    
#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &A_11);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_12);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_13);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_21);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_22);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_23);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_31);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_32);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_33);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_1);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_2);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_3);
#endif  
        
        B_2 += vol_incr;
        B_3 += y_incr;
        
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {           
                    (*s_dPHIdt[0])(k,j,i) = (*s_dPHIdt0[0])(k,j,i) + dmuL*(*s_dGLdmuL)(k,j,i) + dmuV*(*s_dGLdmuV)(k,j,i) + dfY*(*s_dGLdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                    (*s_dPHIdt[1])(k,j,i) = (*s_dPHIdt0[1])(k,j,i) + dmuL*(*s_dGSdmuL)(k,j,i) + dmuV*(*s_dGSdmuV)(k,j,i) + dfY*(*s_dGSdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                    (*s_dPHIdt[2])(k,j,i) = (*s_dPHIdt0[2])(k,j,i) + dmuL*(*s_dGVdmuL)(k,j,i) + dmuV*(*s_dGVdmuV)(k,j,i) + dfY*(*s_dGVdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                }
        
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_y = 0.0 + dfY;
        
        for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];
        
    }
    
    else if (dynamics_type == 6)
    {   
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {   
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    (*s_dGLdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGLdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGVdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));       
                    
                    (*s_dGVdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGSdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGSdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                }
            }
        }  
        
        
        double A_11=0; double A_12=0; double A_13=0; 
        double A_21=0; double A_22=0; double A_23=0;
        double A_31=0; double A_32=0; double A_33=0;
        double B_1=0; double B_2=0; double B_3=0;
        double dmuL=0; double dmuV=0; double dfX=0; double divider=0;
        _MU[1]=0;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    double mat_x = (*s_matrix_x)(k,j,i);
                    A_11 += (*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_12 += (*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_13 += mat_x*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_21 += (*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_22 += (*s_dGVdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0))); 
                    A_23 += mat_x*(*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_31 += mat_x*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_32 += mat_x*(*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_33 += SQR(mat_x)*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    B_1  += -(5.0-5.0*SQR(tanh(10*phi_0-5.0)))* (*s_dPHIdt0[0])(k,j,i);
                    B_2  += -(5.0-5.0*SQR(tanh(10*phi_2-5.0)))* (*s_dPHIdt0[2])(k,j,i);
                    B_3  += -mat_x*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))*(*s_dPHIdt0[0])(k,j,i);
                }    
#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &A_11);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_12);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_13);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_21);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_22);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_23);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_31);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_32);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_33);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_1);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_2);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_3);
#endif  
        
        B_2 += vol_incr;
        B_3 += x_incr;
                
        divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        dmuV = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        dfX  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {           
                    (*s_dPHIdt[0])(k,j,i) = (*s_dPHIdt0[0])(k,j,i) + dmuL*(*s_dGLdmuL)(k,j,i) + dmuV*(*s_dGLdmuV)(k,j,i) + dfX*(*s_dGLdmuL)(k,j,i)*(*s_matrix_x)(k,j,i);
                    
                    (*s_dPHIdt[1])(k,j,i) = (*s_dPHIdt0[1])(k,j,i) + dmuL*(*s_dGSdmuL)(k,j,i) + dmuV*(*s_dGSdmuV)(k,j,i) + dfX*(*s_dGSdmuL)(k,j,i)*(*s_matrix_x)(k,j,i);
                    
                    (*s_dPHIdt[2])(k,j,i) = (*s_dPHIdt0[2])(k,j,i) + dmuL*(*s_dGVdmuL)(k,j,i) + dmuV*(*s_dGVdmuV)(k,j,i) + dfX*(*s_dGVdmuL)(k,j,i)*(*s_matrix_x)(k,j,i);
                    
                }
        
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_x = 0.0 + dfX;
        
        for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];
        
    }
    
    else if (dynamics_type == 7)
    {   
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {   
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    (*s_dGLdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGLdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGVdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));       
                    
                    (*s_dGVdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGSdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGSdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                }
            }
        }  
        
        
        double A_11=0; double A_12=0; double A_13=0; double A_14=0; 
        double A_21=0; double A_22=0; double A_23=0; double A_24=0;
        double A_31=0; double A_32=0; double A_33=0; double A_34=0;
        double A_41=0; double A_42=0; double A_43=0; double A_44=0;
        double B_1=0; double B_2=0; double B_3=0; double B_4=0;
	double divider=0; double dmuL=0; double dmuV=0; double dfX=0; double dfY=0;
        _MU[1]=0;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    double mat_x = (*s_matrix_x)(k,j,i);
                    double mat_y = (*s_matrix_y)(k,j,i);
                    A_11 += (*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_12 += (*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_13 += mat_x*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_14 += mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_21 += (*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_22 += (*s_dGVdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0))); 
                    A_23 += mat_x*(*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_24 += mat_y*(*s_dGVdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    A_31 += mat_x*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_32 += mat_x*(*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_33 += SQR(mat_x)*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_34 += mat_x*mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_41 += mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_42 += mat_y*(*s_dGLdmuV)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_43 += mat_x*mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_44 += SQR(mat_y)*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    B_1  += -(5.0-5.0*SQR(tanh(10*phi_0-5.0)))* (*s_dPHIdt0[0])(k,j,i);
                    B_2  += -(5.0-5.0*SQR(tanh(10*phi_2-5.0)))* (*s_dPHIdt0[2])(k,j,i);
                    B_3  += -mat_x*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))*(*s_dPHIdt0[0])(k,j,i);
                    B_4  += -mat_y*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))*(*s_dPHIdt0[0])(k,j,i);
                }    
#ifdef _STK_MPI
	_node->reduce(STK_SUM, STK_DOUBLE, &A_11);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_12);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_13);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_14);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_21);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_22);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_23);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_24);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_31);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_32);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_33);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_34);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_41);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_42);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_43);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_44);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_1);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_2);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_3);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_4);
#endif       
        
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
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {           
                    (*s_dPHIdt[0])(k,j,i) = (*s_dPHIdt0[0])(k,j,i) + dmuL*(*s_dGLdmuL)(k,j,i) + dmuV*(*s_dGLdmuV)(k,j,i) + dfX*(*s_dGLdmuL)(k,j,i)*(*s_matrix_x)(k,j,i) + dfY*(*s_dGLdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                    (*s_dPHIdt[1])(k,j,i) = (*s_dPHIdt0[1])(k,j,i) + dmuL*(*s_dGSdmuL)(k,j,i) + dmuV*(*s_dGSdmuV)(k,j,i) + dfX*(*s_dGSdmuL)(k,j,i)*(*s_matrix_x)(k,j,i) + dfY*(*s_dGSdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                    (*s_dPHIdt[2])(k,j,i) = (*s_dPHIdt0[2])(k,j,i) + dmuL*(*s_dGVdmuL)(k,j,i) + dmuV*(*s_dGVdmuV)(k,j,i) + dfX*(*s_dGVdmuL)(k,j,i)*(*s_matrix_x)(k,j,i) + dfY*(*s_dGVdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                }
        
        _MU[0] = _MU[1] + dmuL;
        _MU[2] = _MU[1] + dmuV;
        Fext_x = 0.0 + dfX;
        Fext_y = 0.0 + dfY;
        
        for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];
        
    }
    
    else if (dynamics_type == 8)
    {   
        for (int i = x_begin; i < x_end; i++)
        {
            for (int j = y_begin; j < y_end; j++)
            {
                for (int k = z_begin; k < z_end; k++)
                {   
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double phi_2 = (*s_PHI[2])(k, j, i);
                    (*s_dGLdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2+Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGLdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGVdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2-Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));       
                    
                    (*s_dGVdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2+Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(Mob_GL/2+Mob_LV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(Mob_GL+Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                    
                    (*s_dGSdmuL)(k,j,i) = - (*s_K_LS)(k,j,i)*(Mob_GL/2-Mob_LS)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_SV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    
                    (*s_dGSdmuV)(k,j,i) = - (*s_K_LS)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_SV)(k,j,i)*(Mob_GL/2-Mob_SV)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_LV)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)))
                    - (*s_K_other)(k,j,i)*(-Mob_GL)*(5.0-5.0*SQR(tanh(10*phi_2-5.0)));
                }
            }
        }  
        
        
        double A_11=0; double A_12=0; double A_13=0;
        double A_21=0; double A_22=0; double A_23=0;
        double A_31=0; double A_32=0; double A_33=0;
        double B_1=0; double B_2=0; double B_3=0; 
        _MU[1]=0;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {
                    double phi_0 = (*s_PHI[0])(k, j, i);
                    double mat_x = (*s_matrix_x)(k,j,i);
                    double mat_y = (*s_matrix_y)(k,j,i);
                    A_11 += (*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_12 += mat_x*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_13 += mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_21 += mat_x*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_22 += SQR(mat_x)*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_23 += mat_x*mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_31 += mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_32 += mat_x*mat_y*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    A_33 += SQR(mat_y)*(*s_dGLdmuL)(k,j,i)*(5.0-5.0*SQR(tanh(10*phi_0-5.0)));
                    B_1  += -(5.0-5.0*SQR(tanh(10*phi_0-5.0)))* (*s_dPHIdt0[0])(k,j,i);
                    B_2  += -mat_x*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))*(*s_dPHIdt0[0])(k,j,i);
                    B_3  += -mat_y*(5.0-5.0*SQR(tanh(10*phi_0-5.0)))*(*s_dPHIdt0[0])(k,j,i);
                }    
#ifdef _STK_MPI
		_node->reduce(STK_SUM, STK_DOUBLE, &A_11);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_12);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_13);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_21);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_22);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_23);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_31);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_32);
        _node->reduce(STK_SUM, STK_DOUBLE, &A_33);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_1);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_2);
        _node->reduce(STK_SUM, STK_DOUBLE, &B_3);
#endif       
        
        double divider = A_11*A_22*A_33-A_11*A_23*A_32-A_12*A_21*A_33+A_12*A_23*A_31+A_13*A_21*A_32-A_13*A_22*A_31;
        double dmuL =  (A_12*A_23*B_3-A_12*A_33*B_2-A_13*A_22*B_3+A_13*A_32*B_2+A_22*A_33*B_1-A_23*A_32*B_1)/divider;
        double dfX  = -(A_11*A_23*B_3-A_11*A_33*B_2-A_13*A_21*B_3+A_13*A_31*B_2+A_21*A_33*B_1-A_23*A_31*B_1)/divider;
        double dfY  =  (A_11*A_22*B_3-A_11*A_32*B_2-A_12*A_21*B_3+A_12*A_31*B_2+A_21*A_32*B_1-A_22*A_31*B_1)/divider;
        
		for (int i = x_begin; i < x_end; i++)
			for (int j = y_begin; j < y_end; j++)
				for (int k = z_begin; k < z_end; k++)
                {           
                    (*s_dPHIdt[0])(k,j,i) = (*s_dPHIdt0[0])(k,j,i) + dmuL*(*s_dGLdmuL)(k,j,i) + dfX*(*s_dGLdmuL)(k,j,i)*(*s_matrix_x)(k,j,i) + dfY*(*s_dGLdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                    (*s_dPHIdt[1])(k,j,i) = (*s_dPHIdt0[1])(k,j,i) + dmuL*(*s_dGSdmuL)(k,j,i) + dfX*(*s_dGSdmuL)(k,j,i)*(*s_matrix_x)(k,j,i) + dfY*(*s_dGSdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                    (*s_dPHIdt[2])(k,j,i) = (*s_dPHIdt0[2])(k,j,i) + dmuL*(*s_dGVdmuL)(k,j,i) + dfX*(*s_dGVdmuL)(k,j,i)*(*s_matrix_x)(k,j,i) + dfY*(*s_dGVdmuL)(k,j,i)*(*s_matrix_y)(k,j,i);
                    
                }
        
        _MU[0] = _MU[1] + dmuL;
        Fext_x = 0.0 + dfX;
        Fext_y = 0.0 + dfY;
        
        for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_M[i][j] = _MU[i]-_MU[j];
        
    }

	else
	{
		ERROR("unknown dynamics_type = " << dynamics_type);
	}

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
