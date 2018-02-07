/* special functions to perform Monte Carlo simulation
   in which the number of atoms can change 
   
   Activate by: GEN+=-D_MC_CHANGE_NUM_ATOMS in compile
*/

int MDFrame::MCstep_fracatom()
{ /* allowing fractional atom
   * atom _NP-1 is fractional: MC_fracatom
   */
    int ret;
    
    ret = MCstep_moveatom_fracatom();

    if(ensemble_type[1]=='P')
    {
        if(drand48()<(1.0/MC_dV_freq))
            ret = MCstep_dV_fracatom();
    }

    if(ensemble_type[0]=='u')
    {
        if(drand48()<(1.0/MC_dN_freq))
            ret = MCstep_dN_fracatom();
    }

    return ret;
}


int MDFrame::MCstep_dN()
{
    int atom_id, i;
    double dE,tmp,xi;

    //INFO("MCstep_dN curstep="<<curstep);

    xi = drand48() - 0.5;
    //xi = 1; /* only try insert atom */
    
    if (xi < 0)
    { /* remove one atom */
        atom_id = (int) floor(drand48()*_NP);
        potential_energyonly();
        dE = -_EPOT_RMV[atom_id];

        tmp = (-dE - MC_mu_ext)/(KB*_TDES) + log(_NP*MC_Lambda3/_H.det());

	INFO_Printf("MC trial remove: dE = %20.12e tmp = %20.12e\n",dE,tmp);
        
        if(tmp>0)
            MC_accept=1;
        else
        {
            if(tmp<-40) MC_accept=0;
            else
            {
                if( drand48() < exp(tmp) )
                    MC_accept=1;
                else
                    MC_accept=0;
            }
        }    
        if(MC_accept)
        { /* remove atom_id */
            _EPOT += dE;
            _EPOT0=_EPOT;
            MC_accept_tot++;

            for(i=atom_id+1;i<_NP;i++)
            {
                _SR[i-1]=_SR[i];
                fixed[i-1]=fixed[i];
                species[i-1]=species[i];
                _EPOT_IND[i-1]=_EPOT_IND[i];
                _EPOT_RMV[i-1]=_EPOT_RMV[i];
                _TOPOL[i-1]=_TOPOL[i];
            }
            _NP--;
        }
    }
    else
    { /* insert one atom */
        atom_id = _NP;

        _SR[atom_id].x = drand48()-0.5;
        _SR[atom_id].y = drand48()-0.5;
        _SR[atom_id].z = drand48()-0.5;
        fixed[atom_id] = 0;
        
        _NP ++;
        
        /* brute force approach */
        NbrList_reconstruct_use_link_list();
        /* smart algorithm */
        //NbrList_reconstruct_use_link_list(particle_id); /* append list */
        
        potential_energyonly();

        dE = _EPOT - _EPOT0;
        
        tmp = (-dE + MC_mu_ext)/(KB*_TDES) - log(_NP*MC_Lambda3/_H.det());
        INFO_Printf("V = %20.12e -log(NP*Lambda3/V) = %20.12e\n",_H.det(), -log(_NP*MC_Lambda3/_H.det()));
	INFO_Printf("MC trial insert: dE = %20.12e (%20.12e) tmp = %20.12e\n",dE,_EPOT_RMV[atom_id],tmp);
        
        if(tmp>0)
            MC_accept=1;
        else
        {
            if(tmp<-40) MC_accept=0;
            else
            {
                if( drand48() < exp(tmp) )
                    MC_accept=1;
                else
                    MC_accept=0;
            }
        }
    
        if(MC_accept)
        {
            _EPOT0=_EPOT;
            MC_accept_tot++;
        }
        else
        {
            _EPOT=_EPOT0;
            _NP--;
            /* brute force approach */
            NbrList_reconstruct_use_link_list();
            /* smart algorithm */
            //NbrList_remove_atom(particle_id); /* remove entry from list */
        }
    }
    
    return MC_accept;
}        


int MDFrame::MCstep_moveatom_fracatom()
{ /* allowing fractional atom
   * atom _NP-1 is fractional: MC_fracatom
   */
    Vector3 dr, ds, dr0, ds0;
    Matrix33 hinv;
    int i;
    double tmp;

    //INFO("MCstep_moveatom curstep="<<curstep);
    
    hinv=_H.inv();
    
    /* displace atom */
    if(MC_atom>=0)
    {
        FATAL("MCstep_moveatom_fracatom only works for MC_atom = -1, i.e. move all atoms at once");
    }
    else
    {/* move all atoms at once */
        dr0.clear();
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]==-1) continue;
            dr.x=(drand48()-0.5)*_TIMESTEP;
            dr.y=(drand48()-0.5)*_TIMESTEP;
            dr.z=(drand48()-0.5)*_TIMESTEP;
            dr0+=dr;
            ds=hinv*dr;
            _VSR[i]=ds;
        }
        /* fix center of mass */
        dr0/=_NP;
        ds0=hinv*dr0;
        for(i=0;i<_NP;i++)
        {
            _VSR[i]-=ds0;
            _SR[i]+=_VSR[i];
        }
        potential_energyonly();
        _EPOT-=_EPOT_RMV[_NP-1]*(1-MC_fracatom);
    }
    
    if(_EPOT<_EPOT0)
        MC_accept=1;
    else
    {
        tmp=(_EPOT0-_EPOT)/(KB*_TDES);
        if(tmp<-40) MC_accept=0;
        else
        {
            if(drand48()<exp(tmp))
                MC_accept=1;
            else
                MC_accept=0;
        }
    }
    
    if(MC_accept)
    {
        _EPOT0=_EPOT;
        MC_accept_tot++;
    }
    else
    {
        if(MC_atom>=0)
        {
            _SR[MC_atom]-=ds;
        }
        else
        {
            for(i=0;i<_NP;i++)
            {
                _SR[i]-=_VSR[i];
            }
            _EPOT=_EPOT0;
        }
    }
    return MC_accept;
}

int MDFrame::MCstep_dV_fracatom()
{
    double tmp, V, dV;
    
    //INFO("MCstep_dV curstep="<<curstep);
    
    V = _H.det();

    dV = MC_dVmax*(2*drand48()-1);
    _H0 = _H;

    tmp = pow((V+dV)/V,1.0/3.0);
    _H *= tmp;

    /* brute force approach */
    NbrList_reconstruct_use_link_list();
    potential_energyonly();
    _EPOT-=_EPOT_RMV[_NP-1]*(1-MC_fracatom);

    tmp=(_EPOT0-_EPOT - MC_P_ext*dV/160.2e3)/(KB*_TDES) + _NP*log((V+dV)/V);
    
    if(tmp>0)
        MC_accept=1;
    else
    {
        if(tmp<-40) MC_accept=0;
        else
        {
            if(drand48()<exp(tmp))
                MC_accept=1;
            else
                MC_accept=0;
        }
    }
    
    if(MC_accept)
    {
        _EPOT0=_EPOT;
        MC_accept_tot++;
    }
    else
    {
        _H = _H0;
        _EPOT = _EPOT0;
        /* brute force approach */
        NbrList_reconstruct_use_link_list();
    }
    return MC_accept;
}

int MDFrame::MCstep_dN_fracatom()
{
    int atom_id;
    double dE,tmp,xi;

    //INFO("MCstep_dN curstep="<<curstep);

    xi = drand48() - 0.5;
    
    if (xi < 0)
    { /* gradually remove one atom */
        atom_id = _NP-1;
        MC_fracatom -= MC_dfrac;
        potential_energyonly();
        _EPOT-=_EPOT_RMV[_NP-1]*(1-MC_fracatom);
        dE = _EPOT - _EPOT0;

        if(fabs(dE+_EPOT_RMV[atom_id]*MC_dfrac)>1e-10)
            FATAL("MC trial remove: dE = "<<dE<<" EPOT_RMV["<<atom_id<<"] = "<<_EPOT_RMV[atom_id]
                  <<" MC_dfrac = "<<MC_dfrac);
        
        tmp = (-dE - MC_mu_ext*MC_dfrac)/(KB*_TDES) + log(_NP*MC_Lambda3/_H.det())*MC_dfrac;

	INFO_Printf("MC trial remove: dE = %20.12e tmp = %20.12e frac=%f\n",dE,tmp,MC_fracatom);
        
        if(tmp>0)
            MC_accept=1;
        else
        {
            if(tmp<-40) MC_accept=0;
            else
            {
                if( drand48() < exp(tmp) )
                    MC_accept=1;
                else
                    MC_accept=0;
            }
        }    
        if(MC_accept)
        { /* remove atom_id */
            _EPOT0=_EPOT;
            MC_accept_tot++;

            if(MC_fracatom==0)
            {
                _NP --;
                /* select a random atom exchange with NP-1 (use NP as buffer)*/
                atom_id = (int) floor(drand48()*_NP);

                _SR[(_NP)]=_SR[atom_id]; fixed[_NP]=fixed[atom_id]; species[_NP]=species[atom_id]; _EPOT_IND[_NP]=_EPOT_IND[atom_id]; _EPOT_RMV[_NP]=_EPOT_RMV[atom_id]; _TOPOL[_NP]=_TOPOL[atom_id];
                _SR[(atom_id-1)]=_SR[_NP-1]; fixed[atom_id-1]=fixed[_NP-1]; species[atom_id-1]=species[_NP-1]; _EPOT_IND[atom_id-1]=_EPOT_IND[_NP-1]; _EPOT_RMV[atom_id-1]=_EPOT_RMV[_NP-1]; _TOPOL[atom_id-1]=_TOPOL[_NP-1];
                _SR[(_NP-1)]=_SR[_NP]; fixed[_NP-1]=fixed[_NP]; species[_NP-1]=species[_NP]; _EPOT_IND[_NP-1]=_EPOT_IND[_NP]; _EPOT_RMV[_NP-1]=_EPOT_RMV[_NP]; _TOPOL[_NP-1]=_TOPOL[_NP];
                    
                MC_fracatom = 1;
            }
            _EPOT = _EPOT0;
            INFO_Printf("MC remove accepted: NP = %d frac = %f\n",_NP,MC_fracatom);            
        }
        else
        {
            MC_fracatom += MC_dfrac;
            _EPOT = _EPOT0;
        }
    }
    else
    { /* gradually insert one atom */
        if(MC_fracatom==1)
        { /* insert new atom at fraction (MC_dfrac) */        
            atom_id = _NP;
            _SR[atom_id].x = drand48()-0.5;
            _SR[atom_id].y = drand48()-0.5;
            _SR[atom_id].z = drand48()-0.5;
            fixed[atom_id] = 0;        
            _NP ++;
            MC_fracatom = MC_dfrac;
            /* brute force approach */
            NbrList_reconstruct_use_link_list();
            /* smart algorithm */
            //NbrList_reconstruct_use_link_list(particle_id); /* append list */
            potential_energyonly();
            _EPOT-=_EPOT_RMV[_NP-1]*(1-MC_fracatom);
            dE = _EPOT - _EPOT0;

            if(fabs((dE-_EPOT_RMV[atom_id]*MC_dfrac)/dE)>1e-10)
                FATAL("MC trial insert: dE = "<<dE<<" EPOT_RMV["<<atom_id<<"] = "<<_EPOT_RMV[atom_id]
                      <<" MC_dfrac = "<<MC_dfrac);
        }
        else
        {
            atom_id = _NP-1;
            MC_fracatom += MC_dfrac;
            potential_energyonly();
            _EPOT-=_EPOT_RMV[atom_id]*(1-MC_fracatom);
            dE = _EPOT - _EPOT0;
        
            if(fabs(dE-_EPOT_RMV[atom_id]*MC_dfrac)>1e-10)
                FATAL("MC trial insert: dE = "<<dE<<" EPOT_RMV["<<atom_id<<"] = "<<_EPOT_RMV[atom_id]
                      <<" MC_dfrac = "<<MC_dfrac);
        }

        /* exact form of this expression needs to be double checked! */
        tmp = (-dE + MC_mu_ext*MC_dfrac)/(KB*_TDES) - log(_NP*MC_Lambda3/_H.det())*MC_dfrac;
        INFO_Printf("V = %20.12e -log(NP*Lambda3/V) = %20.12e\n",_H.det(), -log(_NP*MC_Lambda3/_H.det()));
	INFO_Printf("MC trial insert: dE = %20.12e (%20.12e) tmp = %20.12e frac=%f\n",dE,_EPOT_RMV[atom_id],tmp,MC_fracatom);
        
        if(tmp>0)
            MC_accept=1;
        else
        {
            if(tmp<-40) MC_accept=0;
            else
            {
                if( drand48() < exp(tmp) )
                    MC_accept=1;
                else
                    MC_accept=0;
            }
        }
    
        if(MC_accept)
        {
            INFO_Printf("MC insert accepted: NP = %d frac = %f\n",_NP,MC_fracatom);
            _EPOT0=_EPOT;
            MC_accept_tot++;
        }
        else
        {
            _EPOT=_EPOT0;
            if(MC_fracatom==MC_dfrac)
            {
                _NP--;
                MC_fracatom=1;
                /* brute force approach */
                NbrList_reconstruct_use_link_list();
                /* smart algorithm */
                //NbrList_remove_atom(particle_id); /* remove entry from list */
            }
            else
            {
                MC_fracatom -= MC_dfrac;
            }
        }
    }
    
    return MC_accept;
}        


