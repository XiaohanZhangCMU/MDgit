/* A test version by William, started 12/06/2013 */
void MDPARALLELFrame::stringrelax_parallel_2()
{
    /* parallel string method
     * assuming Broadcast_Atoms have been called
     * every CPU has a copy of all atoms
     */
    int n, i, ipt, j, nin, nout;
    int moveleftend, moverightend, EmaxDomain, yesclimbimage, islengthconstant;
    double s, r2, fr, fm2;
    double lavg0=0, TanMag2, E0, Emax;
    double *Ec, *Ec_global, *Fm, *Fm_global, *Ft, *Ft_global, *dR, *dR_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    FILE *fp;
    int *ind_left, *ind_left_global, *ind_right, *ind_right_global;
    int cgrelaxsteps;
 
    INFO_Printf("stringrelax[%d]: String Relax Parallel 2\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }
    
    s = _nebinterp[myDomain];
    hinv = _H.inv();
    n = constrainatoms[0];
    
    /* Allocate data array along path */
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    _Ec=Ec;
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Fm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force normal to path */
    Fm_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    Ft=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* force along path */
    Ft_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */

    ind_left=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the left */
    ind_left_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    ind_right=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the right */
    ind_right_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    memset(ind_left,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_left_global,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right_global,0,sizeof(int)*(_CHAINLENGTH+1));

    fp = NULL;
    if(myDomain==0) /* Master opens file */
    {
        fp=fopen("stringeng.out","w");
    }

    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5];
    cgrelaxsteps = nebspec[6];  if (cgrelaxsteps<=0) cgrelaxsteps = 10;
    if (yesclimbimage)
        INFO_Printf("String: The climbing image method will be applied after %d steps.\n",equilsteps);

  /* Done by Keonwook?  EmaxDomain seems to be reset before it is used anyway...Commented out by WK 2013-12-06
    // Initialize EmaxDomain depending on (moverightend) and (moveleftend) 
    if(moverightend && !(moveleftend))
        EmaxDomain=0;
    else if(moveleftend && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
    else if(!(moveleftend) && !(moverightend))
        EmaxDomain=_CHAINLENGTH;
  */

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;

    /* start String Relax iteration */
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */

        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
	    /* Original method only allowed for the relaxation of a small fraction
		of the total atoms, while the rest would maintain an interpolated
		position between A and B.  Thus, ever atom is interpolated above,
		while the atoms that are participating in the relaxation are given
		their previous positions from the real coordinates (_Rc0) below */
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
        }

        /* relax surrounding atoms */
        if(nebspec[0]==1) relax(); // NOT WORKING !!

        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
        if(curstep==step0) MDFrame::clearR0();
	/* Get energy and forces */
        MDFrame::call_potential();
        memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
        Ec[myDomain]=_EPOT;
        MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
        /* Now everybody has the entire Ec array */

        /* receive _Rc1, _Rc2 from neighbors */
        nin = nout = 0;
        if(myDomain>0) /* receive _Rc1 from left node, send _Rc0 to left node */
        {
            MPI_Irecv(_Rc1,n*3,MPI_DOUBLE,myDomain-1,
                      MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
            nin ++;
            MPI_Isend(_Rc0,n*3,MPI_DOUBLE,myDomain-1,
                      MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
            nout ++;
        }
        if(myDomain<_CHAINLENGTH) /* receive _Rc2 from right node, send _Rc0 to right node */
        {
            MPI_Irecv(_Rc2,n*3,MPI_DOUBLE,myDomain+1,
                      MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
            nin ++;
            MPI_Isend(_Rc0,n*3,MPI_DOUBLE,myDomain+1,
                      MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
            nout ++;
        }
        MPI_Waitall(nout, outRequests,outStatus);
        MPI_Waitall(nin,  inRequests, inStatus );

        /* _Fc0 is the force on constrained atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc0[i]=_F[ipt];
        }
        
        /* Find max-energy replica for climbing image method */
        Emax=Ec[0]; EmaxDomain=0;
        for(j=1;j<=_CHAINLENGTH;j++)
            if(Ec[j]>Emax)
            {
                Emax=Ec[j]; EmaxDomain=j;
            }
        if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
            INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);
    
        /* calculate tangential vector */
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the left end */
        //if(moveleftend && myDomain==0) Better to have some version of the tangent vector than none
	if(myDomain==0)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc2[i]-_Rc0[i];

            TanMag2 = 0;
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
        /* calculate tangential vector at the right end */
        //if(moverightend && myDomain==_CHAINLENGTH) Better to have some version of the tangent vector than none
	if(myDomain==_CHAINLENGTH)
        {
            for(i=0;i<n;i++) // n: number of constrained atoms
                _Tan[i] = _Rc0[i]-_Rc1[i];
            
            TanMag2 = 0;        
            for(i=0;i<n;i++) // compute magnitude squared of tangential vector
                TanMag2 += _Tan[i].norm2();
        }
                 
        /* orthogonalize forces */
        memset(Fm,0,sizeof(double)*(_CHAINLENGTH+1));
        memset(Ft,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<_CHAINLENGTH))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
            {
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                if (yesclimbimage && (curstep>=equilsteps) && myDomain==EmaxDomain)
		/* We've already subtracted out the tangential component, but for the climbing
		    image we want it to climb to the highest point.  This is probably located
		    in the direction opposite where the force was pushing the configuration.
		    Thus, by subtracting the force a second time, we get a force in the opposite
		    direction as the original force, pushing the configuration up the energy hill. */
                    _Fc0[i]-=(_Tan[i]*fr);
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        /* orthogonalize forces at the end */
        if((moverightend && myDomain==_CHAINLENGTH)
           || (moveleftend && myDomain==0))
        {
            fr=0;            
            for(i=0;i<n;i++)
                fr+=dot(_Fc0[i],_Tan[i]);
            Ft[myDomain] = fr/sqrt(TanMag2);
            fr /= TanMag2; // normalizing
            
            fm2=0;
            for(i=0;i<n;i++) /* orthogonalizing force against tangential vector */
	    /* We are not really orthogonalizing here, as the force itself is
		maintained, while the magnitude is adjusted.  This may be due to the
		fact that the tangents at the end are already not very good and
		allowing the end configuration to free fall down the energy slope
		is not necessarily a bad thing. */
            {
                tmpvec = _Fc0[i];
                _Fc0[i]-=(_Tan[i]*fr);
                fm2 += _Fc0[i].norm2();
                _Fc0[i] = tmpvec;
            }
            Fm[myDomain] = sqrt(fm2); /* store residual force magnitude */
        }
        MPI_Allreduce(Fm,Fm_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Fm,Fm_global,sizeof(double)*(_CHAINLENGTH+1));
        MPI_Allreduce(Ft,Ft_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(Ft,Ft_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate dR */
        memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
        if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
        {
            /* abs(R_j - R_{j-1}): The distance between a configuration "D" and the one
		immediately to the left of "D"*/
            r2=0;
            for(i=0;i<n;i++)
            {
                dr=_Rc0[i]-_Rc1[i];
                r2+=dr.norm2();
            }
            dR[myDomain]=sqrt(r2);
        }
        MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

        /* calculate the initial averaged chain length, lavg0.  Done on first step */
        if(curstep==step0)
        {
            r2=0;
            for(j=1;j<=_CHAINLENGTH;j++)
                r2+=dR[j];
            
            lavg0=r2/_CHAINLENGTH;
        }
        
        if(myDomain==0) /* Master print file */
        {
            if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
            {
                INFO_Printf("curstep = %d\n",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"%d ",curstep);
                for(j=0;j<=_CHAINLENGTH;j++)
                    fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ", j, Ec[j]-Ec[0], Fm[j], Ft[j], dR[j]);
                fprintf(fp,"\n");
                fflush(fp); 
            }
        }

        if(curstep<(step0 + totalsteps))
        {   /* move chains along force */
          /* do not allow other replicas to move if energy of right end is too high */
          if( ((Ec[_CHAINLENGTH]-Ec[0]) < (Emax-Ec[0])*0.5) || (EmaxDomain < 0.5*_CHAINLENGTH) ) {

            if((myDomain>0)&&(myDomain<_CHAINLENGTH))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }
            /* move the ends along force */
            if((moverightend && myDomain==_CHAINLENGTH)
               || (moveleftend && myDomain==0))
                for(i=0;i<n;i++)
                {
		    //ipt = constrainatoms[i+1];
                    //if (fixed[ipt]==1) continue;
                    _Rc0[i]+=_Fc0[i]*_TIMESTEP;
                }

         	//if(curstep%nebspec[1]==0) reparametrization_with_trimming();
           }
           else {

            /* only move right end, no trimming */
            if (moverightend && myDomain==_CHAINLENGTH)
            { /* do CGRelax */
               /* copy _Rc0 to _SR */
               for(i=0;i<n;i++)
               {
                  ipt=constrainatoms[i+1];
                  _SR[ipt]=hinv*_Rc0[i];
               }
               conj_ftol = 1e-4; conj_itmax = 10; conj_fevalmax = 10;
               conj_fixbox = 1;
               relax();
               SHtoR();
               /* copy _R to _Rc0 */
               for(i=0;i<n;i++)
               {
                   ipt=constrainatoms[i+1];
                   _Rc0[i]=_R[ipt];
               }
               //Ec[myDomain]=_EPOT;
               INFO_Printf("Ec[%d]-Ec[0] = %20.12e\n",myDomain,_EPOT-Ec[0]);
            }

	    //if(curstep%nebspec[1]==0) reparametrization_with_trimming();

           }

#if 0
            for(i=0;i<n;i++)
            {
               ipt=constrainatoms[i+1];
               _SR[ipt]=hinv*_Rc0[i];
            }
            MDFrame::call_potential();
            INFO_Printf("[%d] _EPOT = %20.12e  Ec[%d] = %20.12e",myDomain,_EPOT,myDomain,Ec[myDomain]);

            /* need to modifiy this condition to avoid one more call_potential */
            if (_EPOT > Ec[myDomain]) /* if energy increases after steepest descent step, do CGRelax */
#endif
          if( (myDomain > 0) && (myDomain<_CHAINLENGTH) )
          {
            if ( (Ec[myDomain]-(Ec[myDomain-1]+Ec[myDomain+1])*0.5) > 0.5 ) 
            { /* do CGRelax */
               INFO_Printf("[%d] _EPOT = %20.12e  Ec[%d] = %20.12e  Need relax",myDomain,_EPOT,myDomain,Ec[myDomain]);
               /* copy _Rc0 to _SR */
               for(i=0;i<n;i++)
               {
                  ipt=constrainatoms[i+1];
                  _SR[ipt]=hinv*_Rc0[i];
               }
               conj_ftol = 1e-4; conj_itmax = cgrelaxsteps; conj_fevalmax = cgrelaxsteps;
               conj_fixbox = 1;
               relax();
               SHtoR();
               /* copy _R to _Rc0 */
               for(i=0;i<n;i++)
               {
                   ipt=constrainatoms[i+1];
                   _Rc0[i]=_R[ipt];
               }
               //Ec[myDomain]=_EPOT;
               INFO_Printf("Ec[%d]-Ec[0] = %20.12e\n",myDomain,_EPOT-Ec[0]);
            }
          }
         
	   if(curstep%nebspec[1]==0) reparametrization_with_trimming();

	   //if( curstep%nebspec[1]==0)  { /* pasting the content of reparameterization with trimming */
           //} /* end of reparameterization_with_trimming block */
	    
        } // if(curstep<(step0 + totalsteps))
    }

    free(Ec);
    free(Ec_global);
    free(Fm);
    free(Fm_global);
    free(Ft);
    free(Ft_global);
    free(dR);
    free(dR_global);
    free(ind_left);
    free(ind_left_global);
    free(ind_right);
    free(ind_right_global);

    if(myDomain==0) fclose(fp);

    INFO_Printf("stringrelax[%d]: exit\n",myDomain);

}

void MDPARALLELFrame::reparametrization_with_trimming() 
{
    int n, i, ipt, j, nin, nout, right_index;
    int moveleftend, moverightend, EmaxDomain, left_index;
    int yesclimbimage, islengthconstant;
    int lowleft, lowright, leftendDomain, rightendDomain;
    double s, r2, alpha, new_position, new_totalR;
    double Emax;
    double *Ec, *Ec_global, *dR, *dR_global;
    int *ind_left, *ind_left_global, *ind_right, *ind_right_global;
    Vector3 ds0, ds, dr, dR1, dR2, tmpvec;
    Matrix33 hinv;
    
    INFO_Printf("RP[%d]: Reparametrization With Trimming\n", myDomain);

    if( numDomains != (_CHAINLENGTH+1) )
    {
        ERROR("number of cpus("<<numDomains<<") does not match chainlength("<<_CHAINLENGTH<<")+1");
        return;
    }

//INFO_Printf("RP0\n");
    s = _nebinterp[myDomain];
//INFO_Printf("RP0.1\n");
    hinv = _H.inv();
//INFO_Printf("RP0.2\n");
    n = constrainatoms[0];
//INFO_Printf("RP0.3\n");
    
    /* Allocate data array along path */
    //Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* energy along path */
    Ec=_Ec;
    Ec_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    dR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* segment length in Angstrom */
    dR_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
    ind_left=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the left */
    ind_left_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    ind_right=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* configuration to the right */
    ind_right_global=(int *)malloc(sizeof(int)*(_CHAINLENGTH+1)); /* for communication */
    memset(ind_left,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_left_global,0,sizeof(int)*(_CHAINLENGTH+1));
    memset(ind_right_global,0,sizeof(int)*(_CHAINLENGTH+1));

#ifdef _DEBUG
    Ap=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* alpha along path */
    Ap_global=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1)); /* for communication */
//INFO_Printf("RP0.4\n");
#endif


#ifdef _DEBUG   //////////NEW
    /* fix constrained atoms for relaxing surrounding atoms */
    if(nebspec[0]==1)
    {
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }
#endif

//INFO_Printf("RP0.5\n");
    /* if moveleftend = 0, fix the left end replica
                      = 1, free the left end */
    moveleftend = nebspec[2]; moverightend = nebspec[3];
    yesclimbimage = nebspec[4]; islengthconstant = nebspec[5];

#if 0  /* no longer needed since Ec = _Ec,  Wei 2/8/2014 */
#ifdef _DEBUG  //////////////NEW
//INFO_Printf("RP0.6\n");
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];
//    MDFrame::clearR0(); MDFrame::call_potential(); E0=_EPOT;
//INFO_Printf("RP0.7: Energy: %g, Domain: %d\n", E0, myDomain);
        /* get forces, s is reaction coordinate */
        interpCN(s); /* interpolate all atoms for a given s */
//INFO_Printf("RP0.8\n");
        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
	}
//INFO_Printf("RP0.9\n");
        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
//	MDFrame::clearR0(); 
#endif

	MDFrame::call_potential();
//INFO_Printf("RP0.10\n");

    memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));

#ifdef _DEBUG
    memset(Ap,0,sizeof(double)*(_CHAINLENGTH+1));
#endif

    Ec[myDomain]=_EPOT;

#ifdef _DEBUG
INFO_Printf("Energy = %g\n", _EPOT);
#endif

    MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));

#endif /* no longer needed since Ec = _Ec,  Wei 2/8/2014 */

    /* Now everybody has the entire Ec array */

    /* Find max-energy replica for climbing image method */
    Emax=Ec[0]; EmaxDomain=0;
    for(j=1;j<=_CHAINLENGTH;j++)
        if(Ec[j]>Emax)
        {
            Emax=Ec[j]; EmaxDomain=j;
        }
    if(EmaxDomain==0 || EmaxDomain==_CHAINLENGTH)
        INFO_Printf("Warning: Max-energy domain = %d\n",EmaxDomain);
    // Get the lowest energy configurations to the left and right of the
    // maximum energy configuration.
    lowleft = 0; lowright = _CHAINLENGTH;
    for(j=0;j<=EmaxDomain;j++) {
	//if(Ec[j]<Ec[lowleft]) lowleft = j;
        //prevent very shallow region near state A
	if(Ec[j]<Ec[lowleft] + 0.01 ) lowleft = j; /* Wei 2/8/2014 */
    }
    for(j=EmaxDomain;j<=_CHAINLENGTH;j++) {
	if(Ec[j]<Ec[lowright]) lowright = j;
    }
    // We don't necessarily want to redistribute between the two lowest configurations,
    // however.  We want to try to keep the left and right ends near the same energy.
    leftendDomain = lowleft; rightendDomain = lowright;
    while(Ec[leftendDomain+1]-Ec[0] < 0) leftendDomain++;
    while(Ec[rightendDomain-1]-Ec[0] < 0) rightendDomain--;

/*	for(i = 0; i <= _CHAINLENGTH; i++) {
//		ind_left[i] = 0;
		ind_left_global[i] = 0;
//		ind_right[i] = 0;
		ind_right_global[i] = 0;
	}
*/

//INFO_Printf("RP1\n");


    /* receive _Rc1, _Rc2 from neighbors */
    nin = nout = 0;
    if(myDomain>0) /* receive _Rc1 from left node, send _Rc0 to left node */
    {
	MPI_Irecv(_Rc1,n*3,MPI_DOUBLE,myDomain-1,
		  MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
	nin ++;
	MPI_Isend(_Rc0,n*3,MPI_DOUBLE,myDomain-1,
	          MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	nout ++;
    }
    if(myDomain<_CHAINLENGTH) /* receive _Rc2 from right node, send _Rc0 to right node */
    {
	MPI_Irecv(_Rc2,n*3,MPI_DOUBLE,myDomain+1,
		  MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
	nin ++;
	MPI_Isend(_Rc0,n*3,MPI_DOUBLE,myDomain+1,
		  MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	nout ++;
    }
    MPI_Waitall(nout, outRequests,outStatus);
    MPI_Waitall(nin,  inRequests, inStatus );


    /* calculate dR */
    memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
    if((myDomain>0)&&(myDomain<=_CHAINLENGTH))
    {
        /* abs(R_j - R_{j-1}): The distance between a configuration "D" and the one
	    immediately to the left of "D" */
        r2=0;
        for(i=0;i<n;i++)
        {
            dr=_Rc0[i]-_Rc1[i];
            r2+=dr.norm2();
        }
        dR[myDomain]=sqrt(r2);
    }
    MPI_Allreduce(dR,dR_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(dR,dR_global,sizeof(double)*(_CHAINLENGTH+1));

#ifdef _DEBUG
INFO_Printf("RP2\n");
    for(i = 0; i <= _CHAINLENGTH; i++) {
	INFO_Printf("Configuration %d Energy: %g dR: %g \n",i,Ec[i]-Ec[0],dR[i]);
    }
INFO_Printf("moverightend = %d, moveleftend = %d, Energy Difference = %g\n", moverightend,moveleftend,Ec[0]-Ec[_CHAINLENGTH]);
#endif

    /* Redistribution! */
//INFO_Printf("RP2.3\n");
    new_totalR = 0.0;
    for(j=leftendDomain+1;j<=rightendDomain;j++) {
	new_totalR+=dR[j];

#ifdef _DEBUG
	INFO_Printf("Config %d dR: %g, TotalR: %g\n",j,dR[j],new_totalR);
#endif

    }
    new_position = myDomain*new_totalR/_CHAINLENGTH;
    left_index = leftendDomain;
    // The second condition on this while loop should be unnecessary...
    while(new_position > dR[left_index+1] && left_index < _CHAINLENGTH) {
	left_index++;
	new_position -= dR[left_index];
    }


 
//INFO_Printf("RP3\n");
    if (left_index == _CHAINLENGTH) {
	left_index = _CHAINLENGTH-1; 
	right_index = _CHAINLENGTH;
	new_position += dR[right_index];
    } else { right_index = left_index+1; }
    ind_left[myDomain] = left_index; ind_right[myDomain] = right_index;
    MPI_Allreduce(ind_left,ind_left_global,(_CHAINLENGTH+1),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(ind_left,ind_left_global,sizeof(int)*(_CHAINLENGTH+1));
    MPI_Allreduce(ind_right,ind_right_global,(_CHAINLENGTH+1),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    memcpy(ind_right,ind_right_global,sizeof(int)*(_CHAINLENGTH+1));
    /* Now everyone knows what everyone needs */

#ifdef _DEBUG
INFO_Printf("RP3.1: Config 23: Left: %d; Right: %d\n",ind_left[23],ind_right[23]);
#endif

//	for(i = 0; i <= _CHAINLENGTH; i++) {
//		INFO_Printf("Config %d Left: %d; Right: %d\n",i,ind_left[i],ind_right[i]);
//	}
    /* receive _Rc1, _Rc2 from neighbors */
    nin = nout = 0;
//INFO_Printf("RP3.11\n");
    /* receive _Rc1 from left index, send _Rc0 to those who need it as a left index */
    MPI_Irecv(_Rc1,n*3,MPI_DOUBLE,ind_left[myDomain],
		MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
//INFO_Printf("RP3.12\n");
    nin ++;
//INFO_Printf("RP3.13\n");
    for(i = 0; i <= _CHAINLENGTH; i++) {
	if(ind_left[i] == myDomain) { 
//INFO_Printf("RP3.14\n");
	    MPI_Isend(_Rc0,n*3,MPI_DOUBLE,i,
			MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	    nout ++;
	}
    }
//INFO_Printf("RP3.2\n", myDomain);
    /* receive _Rc2 from right index, send _Rc0 to those who need it as a right index */
    MPI_Irecv(_Rc2,n*3,MPI_DOUBLE,ind_right[myDomain],
		MSG_NEB_ATOM,MPI_COMM_WORLD,&inRequests[nin]);
    nin ++;
    for(i = 0; i <= _CHAINLENGTH; i++) {
	if(ind_right[i] == myDomain) { 
	    MPI_Isend(_Rc0,n*3,MPI_DOUBLE,i,
			MSG_NEB_ATOM,MPI_COMM_WORLD,&outRequests[nout]);
	    nout ++;
	}
    }
//INFO_Printf("RP3.3\n", myDomain);
    MPI_Waitall(nout, outRequests,outStatus);
    MPI_Waitall(nin,  inRequests, inStatus );
    /* Now we need to move everything */
//    new_position += dR[left_index];
    alpha = new_position/dR[right_index];

#ifdef _DEBUG
INFO_Printf("RP4: alpha = %g\n", alpha);
    Ap[myDomain] = alpha;
    MPI_Allreduce(Ap,Ap_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ap,Ap_global,sizeof(double)*(_CHAINLENGTH+1));
	for(i = 0; i <= _CHAINLENGTH; i++) {
		INFO_Printf("Config %d Left: %d; Right: %d; Alpha: %g \n",i,ind_left[i],ind_right[i],Ap[i]);
	}
#endif



    for(i=0;i<n;i++) { 
	_Rc0[i] = _Rc1[i]*(1.0-alpha)+_Rc2[i]*alpha;
    }

	/* This only works when we are applying this to every atom */
        /* copy _Rc0 to _SR */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc0[i];
            //_SR[ipt].subint();                    /* 2010/04/09 Keonwook Kang */
	}

#ifdef _DEBUG
        /* clear neighbor list => enforce neighbor list updated at the very 1st step */
	MDFrame::clearR0(); MDFrame::call_potential();
    memset(Ec,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(Ec_global,0,sizeof(double)*(_CHAINLENGTH+1));
    Ec[myDomain]=_EPOT;
INFO_Printf("Energy = %g\n", _EPOT);
    MPI_Allreduce(Ec,Ec_global,(_CHAINLENGTH+1),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    memcpy(Ec,Ec_global,sizeof(double)*(_CHAINLENGTH+1));
    /* Now everybody has the entire Ec array */
    for(i = 0; i <= _CHAINLENGTH; i++) {
	INFO_Printf("Configuration %d Energy: %g\n",i,Ec[i]-Ec[0]);
    }
INFO_Printf("RP5\n");
#endif



#ifdef _DEBUG
/* call_potential first */
/* put this at the end of the function to print out eng to screen */
       
    if(myDomain==0) /* Master print to screen */
    {
	for(j=0;j<=_CHAINLENGTH;j++) INFO_Printf("%8d %25.15e\n", j, Ec[j]-Ec[0]);
    }

    /* Everyone print to screen */
    INFO_Printf("CPU[%d] needs config %d (left) and %d (right)\n", myDomain, ind_left[myDomain], ind_right[myDomain]);
#endif



    //free(Ec);  /* Wei 2/8/2014 */
    free(Ec_global);
    free(dR);
    free(dR_global);
    free(ind_left);
    free(ind_left_global);
    free(ind_right);
    free(ind_right_global);
}


