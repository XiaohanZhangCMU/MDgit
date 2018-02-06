#include "wire.h"

//modify linkfibers and makefibers so that one thread and one cell is created
void WIREFrame::linkfibers()
{
  int fiber_type, nfiber, chainlen;
  int ipt, i, j;
  double bond_len;

  fiber_type = (int) input[0]; /* 0: loop (over PBC), 1: segment */
  nfiber     = (int) input[1];
  chainlen   = (int) input[2];
  bond_len   = input[3];

  for(ipt=0;ipt<_NP;ipt++) {
    /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
    if (species[ipt] != 0) 
      continue;

    group[ipt] = 0;
  }
  if (fiber_type == 0 )
    { /* loop type */
      ERROR("fiber_type "<<fiber_type<<" not implemented!");
      return;
    }
  ipt = 0;
  for(i=0;i<nfiber;i++)
    {
      for(j=0;j<chainlen;j++)
	{
	  if(j==0)
	    {
	      /* only one bond */
	      num_bonds[ipt] = 1;
	      bond_index[ipt*MAXNUMBOND + 0] = ipt+1;
               
	    } else if (j==(chainlen-1))
	    {   /* only one bond */
	      num_bonds[ipt]=1;
	      bond_index[ipt*MAXNUMBOND + 0] = ipt-1;
	    } else
	    {   /* two bonds */
	      num_bonds[ipt]=2;
	      bond_index[ipt*MAXNUMBOND + 0] = ipt-1;
	      bond_index[ipt*MAXNUMBOND + 1] = ipt+1;
	    }
	  group[ipt] = j;
	  ipt ++;
	}
    }
#if 0
  //link of cell
  num_bonds[ipt] = 1;
  bond_index[ipt*MAXNUMBOND+0] = ipt+1;
  group[ipt] = 0;
  
  ipt++;
  num_bonds[ipt] = 1;
  bond_index[ipt*MAXNUMBOND+0] = ipt-1;
  group[ipt] = 1;
#endif
}

void WIREFrame::makefibers()
{
  printf("make fiber of wire is called\n");
    int fiber_type, nfiber, chainlen, fix_bottom, flatten;
    int ipt, i, j;
    double bond_len,d_by_two,offset,orient_pref, bottom_coord;
    Vector3 dr, ds, ddr, dds;
    Matrix33 hinv;

    fiber_type = (int) input[0]; /* 0: loop (over PBC), 1: segment */
    nfiber     = (int) input[1];
    chainlen   = (int) input[2];
    bond_len   = input[3];

    // hardcode
    assert(nfiber == 1);
#if 0
    int ncells = 1;
    assert(ncells == 1);
#endif

    d_by_two=0.04;
    offset=0.2;
    
    _H.clear();
    _H[0][0] = input[4];
    _H[1][1] = input[5];
    _H[2][2] = input[6];

    orient_pref = input[7];
    fix_bottom = (int)input[8];
    bottom_coord = input[9];
    flatten = (int)input[10];
#if 0
    _NP = nfiber*chainlen + ncells*2;
#else
    _NP = nfiber*chainlen;
#endif
    Alloc();

    hinv = _H.inv();
    ipt = 0;
    for(i=0;i<nfiber;i++)
    {
        for(j=0;j<chainlen;j++)
        {
           if(j==0)
           {
               /* the thread initialized on the top*/
	     if (flatten == 1) _SR[ipt].set(0.075,0,0.5);
	     else              _SR[ipt].set(0.075,0,0.5);
	     if (fix_bottom) _SR[ipt].y = bottom_coord;
               /* choose a growth in z direction */
	     dr.set(0,0,-1);
	     dr.y += orient_pref;
	     if (flatten == 1) dr.y = 0.0;
	     dr/=dr.norm();  dr*=bond_len; ds = hinv*dr;
           } 
           else 
           {
               /* randomly perturb the growth direction */	     
	     ddr.set(drand48()-0.5,drand48()-0.5,drand48()-0.5);
	     ddr/=ddr.norm();  ddr*=bond_len*0.1; dds = hinv*ddr;
               /* add a new atom in the randomly chosen direction */
	     _SR[ipt] = _SR[ipt-1] + ds + dds;	     
           }
           species[ipt] = 0;
           ipt ++;
        }
    }

#if 0
    // This initialize a blood cell in the center
    if (flatten == 1) _SR[ipt].set(0.075,0.02,0.08);
    else              _SR[ipt].set(0.1,0,0.1);
    if (fix_bottom) _SR[ipt].y = bottom_coord;
    species[ipt] = 0;

    ipt ++;
    if (flatten == 1) _SR[ipt].set(0.075,-0.02,0.08);
    else              _SR[ipt].set(0.1,0,0.1);
    if (fix_bottom) _SR[ipt].y = bottom_coord;
    species[ipt] = 0;
#endif

//    linkfibers();
}


//check the energy is the sum of the two 
void WIREFrame::fiber_fem_interaction()
{
      /* no need of neighbor list */
  int i,ipt,iele,j,jpt, ia;
  Vector3 ds, dr, fint, s0, savg;
  double  Eint;

    //    int map23[2] = {0, 2};

    DUMP("FEM");

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;
    for(i=0;i<_NP;i++)
    {
      _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear();
    }
    _VIRIAL.clear();

    for(ia=0;ia<_NANCHOR;ia++)
    {
      ipt = atom_fem_links[ia*2];
      iele = atom_fem_links[ia*2+1];

      jpt=elements[iele*_NNODE_PER_ELEMENT+0];
      s0=_SR[jpt];

      savg.clear();      
      for(j=0;j<_NNODE_PER_ELEMENT;j++)
	{
	  jpt=elements[iele*_NNODE_PER_ELEMENT+j];
          ds = _SR[jpt] - s0; ds.subint();
	  
          ds *= atom_fem_weights[ia*_NNODE_PER_ELEMENT+j];
          savg += ds;
        }  
      savg += s0;

      /* distance between atom and fem_anchor_point */
      ds = _SR[ipt] - savg;  ds.subint();
      ds.y = 0;

      _H.multiply(ds, dr);
      dr.y = 0;

      Eint = 0.5*_K_ANCHOR*dr.norm2();
      fint = dr*_K_ANCHOR;

      _EPOT += Eint;
      _F[ipt] -= fint;
      
      for(j=0;j<_NNODE_PER_ELEMENT;j++)
	{
	  jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	  _F[jpt] += fint * atom_fem_weights[ia*_NNODE_PER_ELEMENT+j];
        }  
    }
}

void WIREFrame::potential()
{
  clcvars();
  clvars();  LJBONDFrame::potential();
  accumvars(0); 
  clvars(); 
  FEMFrame::potential(); 
  accumvars(1);
  if (beads_are_docked) { clvars();  fiber_fem_interaction(); accumvars(2); }

  swapvars();
}

//clear inherited variables
void WIREFrame::clvars()
{
  for(int i=0;i<_NP;i++) { 
    _F[i].clear();  _EPOT_IND[i] = 0; 
    _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear(); 
  } _VIRIAL.clear(); _EPOT = 0;
}

//clear cached variables
void WIREFrame::clcvars()
{
  for(int i=0;i<_NP;i++) { 
    cached_F[i].clear();  cached_EPOT_IND[i] = 0; 
    cached_EPOT_RMV[i]=0; cached_VIRIAL_IND[i].clear(); 
  } cached_VIRIAL.clear();   cached_EPOT = 0;
}

void WIREFrame::swapvars()
{
  for(int i=0;i<_NP;i++) { 
    _F[i] = cached_F[i];  _EPOT_IND[i] = cached_EPOT_IND[i]; 
    _EPOT_RMV[i] = cached_EPOT_RMV[i]; _VIRIAL_IND[i] = cached_VIRIAL_IND[i]; 
  } _VIRIAL = cached_VIRIAL; _EPOT =  cached_EPOT;
}

void WIREFrame::accumvars(int config_id)
{
  int start, end;
  // if (config_id == 1) {
  //   start = n_existed_atoms; end = _NP;
  // } else if (config_id == 2) {
  //   start = 0; end= _NP;
  // }  else {
  // start = 0; end = n_existed_atoms;
  // }


  start = 0; end = _NP;
  INFO_Printf("accum: %d, %d, %d, %d, %d\n",config_id, start,end,n_existed_atoms,_NP);
  for(int i=start;i<end;i++) { 
    cached_F[i] += _F[i]; cached_EPOT_IND[i] += _EPOT_IND[i]; 
    cached_EPOT_RMV[i] += _EPOT_RMV[i]; cached_VIRIAL_IND[i] += _VIRIAL_IND[i];
  } cached_VIRIAL += _VIRIAL; cached_EPOT += _EPOT;
}

/*
  helper1: select_bead(int). given a bead, judge whether is should be docked
  helper2: belongs_to(int). find the element that it belongs to. 
                            if the bead is on the interface of two or three elements, 
			    pick the first
  Both helpers can be modified. and need to be visually checked 
 */
void WIREFrame::dockbeads()
{
  int bead_start, bead_end;
  FILE *dockfile = fopen("dockfile","w");
  assert(dockfile != NULL);

  //INFO_Printf("nfemNodes = %d, n_existed_atoms = %d, _NP = %d\n", nfemNodes, n_existed_atoms, _NP);
  if ( nfemNodes == n_existed_atoms) { //fem is read in first
    bead_start = n_existed_atoms; bead_end = _NP;
  } 
  if ( nfemNodes == _NP - n_existed_atoms ) { //fiber is read in first
    bead_start = 0; bead_end = n_existed_atoms;
  }

  int cnt = 0;

  INFO_Printf("bead_start=%d, bead_end = %d, n_existed_atoms=%d, NP = %d\n", bead_start, bead_end, n_existed_atoms,_NP);
  for (int i = bead_start; i<bead_end; i++) {
    if ( select_bead (i, bead_end-bead_start) ) {
      //int eleid = belongs_to (i);
      int eleid; double weights[4];
      eleid = get_weight (i, weights);
      
      assert(eleid != -1);

      fprintf(dockfile, "%d  %d %20.15E %20.15E %20.15E %20.15E\n", i, eleid, weights[0], weights[1], weights[2], weights[3]);

      //temperary. should read from file
      atom_fem_links[2*cnt] = i;
      atom_fem_links[2*cnt+1] = eleid;

      atom_fem_weights[_NNODE_PER_ELEMENT*cnt]   = weights[0];
      atom_fem_weights[_NNODE_PER_ELEMENT*cnt+1] = weights[1];
      atom_fem_weights[_NNODE_PER_ELEMENT*cnt+2] = weights[2];
      atom_fem_weights[_NNODE_PER_ELEMENT*cnt+3] = weights[3];
      cnt ++;      
    }
  }  
  INFO_Printf("_NANCHOR = %d, cnt = %d\n", _NANCHOR, cnt);
  assert(cnt == _NANCHOR);
  fclose(dockfile);
  
  beads_are_docked = 1;
}


bool WIREFrame::select_bead( int bead_id, int nbeads)
{
  //  return !(bead_id % (nbeads/_NANCHOR) );
  return true;

  if(bead_id<nbeads-2) {
    INFO_Printf("beadId = %d, nbeads-2 = %d\n", bead_id, nbeads-2);
    return true;
  }
  else
    return false;
}

int WIREFrame::belongs_to (int bead_id)
{
  double xmax, xmin, zmax, zmin;
  Vector3 srj, sr;
  sr = _SR[bead_id];
  sr.subint();
  double bdx = sr.x;
  double bdz = sr.z;

  for (int iele = 0;iele<_NELE; iele++) {
    if (_NDIM == 2)  {
      xmax = -1e5; xmin = 1e5; zmax = -1e5; zmin = 1e5;
      for (int j = 0;j<_NNODE_PER_ELEMENT;j++) {
	int jpt = elements[iele*_NNODE_PER_ELEMENT+j];
	srj = _SR[jpt];
	srj.subint();
	xmax = srj.x >= xmax ? srj.x : xmax;
	xmin = srj.x <= xmin ? srj.x : xmin;
	zmax = srj.z >= zmax ? srj.z : zmax;
	zmin = srj.z <= zmin ? srj.z : zmin;
	//INFO_Printf("srj.x = %lf, srj.z = %lf, xmax = %lf, xmin = %lf, zmax = %lf, zmin = %lf, bdx = %lf, bdz = %lf\n", srj.x, srj.z, xmax,xmin,zmax,zmin,bdx,bdz);
      }
      if (bdx <= xmax && bdx >= xmin && bdz <= zmax && bdz >=zmin )
	return iele;
    }
  }
  return -1;
}

// this should be moved to the matlab code
int WIREFrame::get_weight (int bead_id, double* weight)
{
  double xmax, xmin, zmax, zmin;
  double xi, zi;
  Vector3 sr, srj;
  sr = _SR[bead_id];
  sr.subint();
  double bdx = sr.x;
  double bdz = sr.z;
  double x0, x1, z0, z3;


  for (int iele = 0;iele<_NELE; iele++) {
    if (_NDIM == 2)  {
      xmax = -1e5; xmin = 1e5; zmax = -1e5; zmin = 1e5;
      
      for (int j = 0;j<_NNODE_PER_ELEMENT;j++) {
	int jpt = elements[iele*_NNODE_PER_ELEMENT+j];
	srj = _SR[jpt];
	srj.subint();
	xmax = srj.x >= xmax ? srj.x : xmax;
	xmin = srj.x <= xmin ? srj.x : xmin;
	zmax = srj.z >= zmax ? srj.z : zmax;
	zmin = srj.z <= zmin ? srj.z : zmin;

	//INFO_Printf("srj.x = %lf, srj.z = %lf, xmax = %lf, xmin = %lf, zmax = %lf, zmin = %lf, bdx = %lf, bdz = %lf\n", srj.x, srj.z, xmax,xmin,zmax,zmin,bdx,bdz);

	if (j == 0) {
	  x0 = srj.x; z0 = srj.z;
	}
	if (j == 1) {
	  x1 = srj.x;
	}
	if (j == 3) {
	  z3 = srj.z;
	}
      }
           
      if (bdx <= xmax && bdx >= xmin && bdz <= zmax && bdz >=zmin ) {
        xi = (bdx-x0)/(x1-x0)*2.0 - 1.0;
        zi = (bdz-z0)/(z3-z0)*2.0 - 1.0;

	weight[0] = 1.0/4.0*(1-xi)*(1-zi);
	weight[1] = 1.0/4.0*(1+xi)*(1-zi);
	weight[2] = 1.0/4.0*(1+xi)*(1+zi);
	weight[3] = 1.0/4.0*(1-xi)*(1+zi);
	return iele;
      }
    }
  }
  return -1;
}

void WIREFrame::plot()
{
    int ipt,jpt,j;
    double L;
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    Vector3 sri, srj, dsrij, ri, rj, sij;

    MDPARALLELFrame::plot();

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=max(_H[0][0],_H[1][1]);
    L=max(L,_H[2][2])*.5;

    SHtoR();
    win->Lock();
    //win->Clear();

    /* Draw Bonds (from connectivity matrix) */
    for(ipt=0;ipt<_NP;ipt++)
    {
       for(j=0;j<num_bonds[ipt];j++)
       {
          jpt = bond_index[ipt*MAXNUMBOND + j];
          if(ipt>jpt) continue;

          sri=_SR[ipt];
          if(plot_map_pbc==1) sri.subint();
          ri = _H*sri;
            
          srj=_SR[jpt];
          dsrij = srj - sri;  dsrij.subint();
          srj = sri + dsrij;
          rj = _H*srj;

          if(plot_limits[0])
              if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
               ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
               ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
               ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                        continue;

          x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
          x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
          dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
          dx/=dr;dy/=dr;dz/=dr;

          win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                        y1+dy*atomradius[species[ipt]]/L,
                        z1+dz*atomradius[species[ipt]]/L,
                        x2-dx*atomradius[species[jpt]]/L,
                        y2-dy*atomradius[species[jpt]]/L,
                        z2-dz*atomradius[species[jpt]]/L,
                        colors[MAXCOLORS+1],bondradius/L,1);
       }
    }
    /* draw covalent bonds between atoms */
    // FEMFrame::plot();    


  _NDIM = 2;
  _NNODE_PER_ELEMENT = 4;
  if (elements != NULL)
    for(int iele=0;iele<_NELE;iele++)
    {
      if (_NDIM == 2)
        for(j=0;j<_NNODE_PER_ELEMENT;j++)
        {	  
             jpt=elements[iele*_NNODE_PER_ELEMENT+j];
             ipt=elements[iele*_NNODE_PER_ELEMENT+((j+1)%_NNODE_PER_ELEMENT)];

             sri=_SR[ipt];
             if(plot_map_pbc==1) sri.subint();
             ri = _H*sri;

             sij=_SR[jpt]-_SR[ipt];
             sij.subint();
             srj=sri+sij;
             rj = _H*srj;
                
             if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                     continue;

                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
                    x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
                    dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
                    dx/=dr;dy/=dr;dz/=dr;
                    win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                                  y1+dy*atomradius[species[ipt]]/L,
                                  z1+dz*atomradius[species[ipt]]/L,
                                  x2-dx*atomradius[species[jpt]]/L,
                                  y2-dy*atomradius[species[jpt]]/L,
                                  z2-dz*atomradius[species[jpt]]/L,
                                  colors[MAXCOLORS+1],bondradius/L,1);
        }
    }
    
    win->Unlock();
    win->Refresh();

}

void WIREFrame::Alloc()
{
  FEMFrame::Alloc();
  LJBONDFrame::Alloc();
  int size = _NP * allocmultiple;
  Realloc(cached_F, Vector3, size);
  Realloc(cached_EPOT_IND, double, size);
  Realloc(cached_EPOT_RMV, double, size);
  Realloc(cached_VIRIAL_IND, Matrix33, size);

  _NNODE_PER_ELEMENT = 4; // try to get rid of it
  Realloc(atom_fem_links, int, (_NANCHOR*2));
  Realloc(atom_fem_weights, double, (_NANCHOR*_NNODE_PER_ELEMENT));
}

void WIREFrame::initvars()
{
  FEMFrame::initvars();
  LJBONDFrame::initvars();
}

void WIREFrame::initparser()
{
  FEMFrame::initparser();
  //LJBONDFrame::initparser();

    int i, j;
    char s[100];
    MDPARALLELFrame::initparser();

    /* input */
    for(i=0;i<MAXSP;i++)
       for(j=i;j<MAXSP;j++)
       {
           sprintf(s,"C12_%d%d",i,j);
           bindvar(s,&(_ALJ[i][j]),DOUBLE);
           sprintf(s,"C6_%d%d",i,j);
           bindvar(s,&(_BLJ[i][j]),DOUBLE);
       }

    bindvar("BOND_R0",&BOND_R0,DOUBLE);
    bindvar("BOND_K", &BOND_K, DOUBLE); /* BOND_K = E * A / L */
    bindvar("BOND_B", &BOND_B, DOUBLE); /* BOND_B = E * I / L */
    bindvar("Rcut",&LJ_RC,DOUBLE);
    bindvar("usrfile",usrfile,STRING);

    bindvar("substrate_Y",  &_SUBSTRATE_Y,DOUBLE);
    bindvar("substrate_REP",&_SUBSTRATE_REP,DOUBLE);
    bindvar("substrate_ATR",&_SUBSTRATE_ATR,DOUBLE);
    bindvar("substrate_CON",&_SUBSTRATE_CON,DOUBLE);

    bindvar("NANCHOR",&_NANCHOR,INT);
    bindvar("K_ANCHOR",&_K_ANCHOR,DOUBLE);

}

int WIREFrame::exec(const char* name)
{
  if(FEMFrame::exec(name)==0) return 0;
  //  if(LJBONDFrame::exec(name)==0) return 0;
  bindcommand(name,"initLJ",LJBONDFrame::initLJ());
  bindcommand(name,"makefibers",makefibers());
  bindcommand(name,"linkfibers",linkfibers());
  bindcommand(name,"writeatomeyeusr",LJBONDFrame::writeatomeyeusrfile(usrfile));
  bindcommand(name,"dockbeads",dockbeads());

  return -1;
}

#ifdef _TEST

/* Main Program Begins */
class WIREFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
