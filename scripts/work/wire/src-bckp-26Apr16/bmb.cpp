/*
  bmb.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Mon Aug 27 11:49:40 2007

  FUNCTION  :  BMB potential (for YSZ)
  
  Test cases: scripts/work/ysz/zro2.tcl
*/

#include "bmb.h"

void BMBFrame::initparser()
{
    MDPARALLELFrame::initparser();
}

int BMBFrame::exec(const char *name)
{
    if(MDPARALLELFrame::exec(name)==0) return 0;

    bindcommand(name,"dope_Yttria",dope_Yttria());    
    return -1;
}

void BMBFrame::Alloc()
{
    MDPARALLELFrame::Alloc();
}

void BMBFrame::initvars()
{
    MDPARALLELFrame::initvars();
}

void BMBFrame::calcprop()
{
    MDPARALLELFrame::calcprop();
}

void BMBFrame::born_meyer_buckingham()
{
    int i,j,ipt,jpt,isp,jsp,k;
    double U,r,r2,ri6, cx1,cx2,x,xcube;
    Vector3 sij, rij, fij;
    
    DUMP(HIG"Lennard Jones"NOR);
        
    refreshneighborlist();

    BMB_Rc = Ewald_Rcut;
    
    _EPOT=0;

    for(i=0;i<_NP;i++)
    {_F[i].clear(); _EPOT_IND[i]=0;}
    _VIRIAL.clear();

    for(ipt=0;ipt<_NP;ipt++)
    {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[ipt]==-1) continue;
        
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            if(ipt>jpt) continue;

            sij=_SR[jpt]-_SR[ipt];
            sij.subint();
            rij=_H*sij;
            r2=rij.norm2();
            r=sqrt(r2);
            if(r<=BMB_Rc)
            {
              ri6=1./(r2*r2*r2);

              // calculate the cutoff variables
              if(r>(BMB_Rc-1.0))
              {
                  x=r-BMB_Rc+1.0;
                  xcube=CUBE(x);
                  cx1=SQR(1.0-xcube*x);
                  cx2=8.0*(1.0-xcube*x)*xcube;
              }
              else
              {
                  cx1=1.0;
                  cx2=0.0;
              }
              
              // get the species indices
              isp=species[ipt]; // species of ipt-atom
              jsp=species[jpt]; // species of jpt-atom
              k=isp+jsp; // interaction index
              
              // calculate potential
              U=(-BMB_POT[k][0]*ri6+BMB_POT[k][1]*exp(-BMB_POT[k][2]*r))*cx1;
              
              // calculate force
              fij=rij*(((-6.0*BMB_POT[k][0]*ri6/r+BMB_POT[k][1]*BMB_POT[k][2]
                          *exp(-BMB_POT[k][2]*r))*cx1+U*cx2)/r);
                            
              _F[ipt]-=fij;
              _F[jpt]+=fij;
              _EPOT_IND[ipt]+=U*0.5;
              _EPOT_IND[jpt]+=U*0.5;
              _EPOT+=U;
              _VIRIAL.addnvv(1.,fij,rij);
            }
        }
    }
}

void BMBFrame::born_meyer_buckingham_Ewald()
{
    int i;
    
    born_meyer_buckingham();

    /* Classical Ewald */
    if(Ewald_CE_or_PME==0){CE();}
    /* Paricle-Mesh Ewald */
    else{PME();}

    /* Sum up energy, force and stress */
    _EPOT += _EPOT_Ewald;
    _VIRIAL += _VIRIAL_Ewald;
    for(i=0;i<_NP;i++)
    {
        _F[i]+=_F_Ewald[i];
        _EPOT_IND[i]+=_EPOT_IND_Ewald[i];
    }

    /* check the results */
    /*
    INFO("==============================================================================");
    INFO(" _EPOT = "<<_EPOT);
    INFO("------------------------------------------------------------------------------");
    INFO(" _EPOT_BMB = "<<_EPOT-_EPOT_Ewald);
    INFO("------------------------------------------------------------------------------");
    INFO(" _EPOT_Ewald = "<<_EPOT_Ewald);
    INFO("------------------------------------------------------------------------------");
    INFO(" _EPOT_Ewald_Real = "<<_EPOT_Ewald_Real);
    INFO(" _EPOT_Ewald_Rec = "<<_EPOT_Ewald_Rec);
    INFO(" _EPOT_Ewald_Self = "<<_EPOT_Ewald_Self);
    INFO("------------------------------------------------------------------------------");
    INFO(" _VIRIAL_Ewald = "<<_VIRIAL_Ewald);
    INFO("------------------------------------------------------------------------------");
    INFO(" _VIRIAL_Ewald_Real = "<<_VIRIAL_Ewald_Real);
    INFO(" _VIRIAL_Ewald_Rec = "<<_VIRIAL_Ewald_Rec);
    INFO("==============================================================================");
    */
}

void BMBFrame::potential()
{
    born_meyer_buckingham_Ewald();
}

void BMBFrame::dope_Yttria()
{
    /* dope perfect ZrO2 crystal with Y2O3 */    
    int Y_count,V_count,Max_Y,Max_V;
    int N_cat,N_an,N_cat_set,N_an_set;
    int I_cat_1,I_cat_2,I_an_1,I_an_2;
    double yttriaconc, F_vac, rand_temp;

    /* input Yttria concentration */
    yttriaconc = input[0];

    if(yttriaconc==0.0){WARNING("DOPEYTTRIA : YTTRIA CONCENTRATION IS NOT SPECIFIED"); return;}
    N_cat = _NP/12*4; // number of cations
    N_an = _NP/12*8; // number of anions
    F_vac = yttriaconc/(2.0+2.0*yttriaconc); // vacancy fraction (N_vac / N_an)
    Max_V = (int) floor(F_vac*((double)N_an)+0.5); // number of vacancy
    Max_Y = Max_V * 2; // number of Yttrium
    Y_count=0; // initialize Yttrium number
    V_count=0; // initialize vacancy number
    INFO("------------------------------------------------------------------------------");
    INFO("------------------------------------------------------------------------------");
    INFO("                             Doping Yttria [PL]");
    INFO("        Substitute Zirconium with Yttrium and Create Oxygen Vacancies");
    INFO("------------------------------------------------------------------------------");
    INFO("    Total Number of Atomic Sites  = "<<(_NP+Max_V));
    INFO("    Total Number of Cation Sites  = "<<N_cat);
    INFO("    Total Number of Anion Sites   = "<<N_an);
    INFO("    Vacancy Fraction              = "<<F_vac);
    INFO("------------------------------------------------------------------------------");
    INFO("    Total Number of Zirconiums    = "<<(N_cat-Max_Y));
    INFO("    Total Number of Oxygens       = "<<N_an-Max_V);
    INFO("    Total Number of Yttriums      = "<<Max_Y);
    INFO("    Total Number of Vacancies     = "<<Max_V);
    INFO("------------------------------------------------------------------------------");
    INFO("------------------------------------------------------------------------------");

    // dope Yttrium into cation sites
    while(Y_count<Max_Y)
    {
        // pick a random number between 0 and 1
        rand_temp = ((double)drand48());
        if(rand_temp!=1.0)
        {
            // get the atom index of a random cation
            I_cat_1 = (int) floor(((double)N_cat)*rand_temp);
            N_an_set = I_cat_1/4; // implicitly floor
            I_cat_2 = I_cat_1+8*N_an_set; // final cation index

            // check the species of picked atom
            if(species[I_cat_2]==0) // if target is Zr
            {
                species[I_cat_2]=3; // set the species index to 3
                Y_count++; // increase the Yttrium count
            }
        }
    }
    // free all atoms
    freeallatoms();

    // create vacancy
    while(V_count<Max_V)
    {
        // pick a random number between 0 and 1
        rand_temp = ((double)drand48());
        if(rand_temp!=1.0)
        {
            // get the atom index of a random anion
            I_an_1 = (int) floor(((double)N_an)*rand_temp);
            N_cat_set = I_an_1/8+1;
            I_an_2 = I_an_1+4*N_cat_set;

            if(species[I_an_2]==1) // if target is oxygen
            {
                fixed[I_an_2]=1; // fix picked oxygen atom
                V_count++; // count number of picked atoms
            }
        }
    }
    
    // create oxygen vacancies
    if (input[1] == -1 )
        markremovefixedatoms();
    else
    {
        removefixedatoms();
    }
}

#ifdef _TEST

/* Main Program Begins */
class BMBFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

