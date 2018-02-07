#include "xfem.h"

/*
  This read2cn() is specific for fiber+FEMsubstrate
  It first dumps the already existed configuration, and combine the dumped file with the current "incnfile"  
 */
int XFEMFrame::read2cn()
{
  assert(_NP>0); 
  int ival1, ival2; 
  double dvalx, dvaly, dvalz;
  n_existed_atoms = _NP;

  strcpy(finalcnfile,  "config1.cn");
  writefinalcnfile(0,false);

  FILE *config1fp  = fopen(finalcnfile,"r");
  FILE *config2fp  = fopen(incnfile,"r");   
  FILE *config12fp = fopen("config12.cn","w");

  assert (config1fp!=NULL); 
  assert (config2fp!=NULL);
  assert (config12fp!=NULL);

  fscanf(config1fp, "%d", &ival1);
  fscanf(config2fp, "%d", &ival2);  

  nfemNodes = ival2;
  fprintf(config12fp, "%d\n" ,ival1+ival2);

  for (int nd = 0;nd<ival1;nd++) {
    fscanf(config1fp, "%lf %lf %lf" , &dvalx, &dvaly, &dvalz);
    fprintf(config12fp, "%23.18E %23.18E %23.18E\n" ,dvalx, dvaly, dvalz);
  }
  for (int nd = 0;nd<ival2;nd++) {
    fscanf(config2fp, "%lf %lf %lf" , &dvalx, &dvaly, &dvalz);
    fprintf(config12fp, "%23.18E %23.18E %23.18E\n" ,dvalx, dvaly, dvalz);
  }
  //box information is the same as the existed one
  for (int nd = 0;nd<3;nd++) {
    fscanf(config1fp, "%lf %lf %lf" , &dvalx, &dvaly, &dvalz);
    fprintf(config12fp, "%23.18E %23.18E %23.18E\n" ,dvalx, dvaly, dvalz);
  }

  fclose(config1fp);
  fclose(config2fp);
  fclose(config12fp);

  //perform a regular readcn on the combined cn file
  strcpy(incnfile,  "config12.cn");
  return (FEMFrame::readcn());
}

void XFEMFrame::initparser()
{
  FEMFrame::initparser();
  bindvar("fem_bdy_nodes_file",fem_bdy_nodes_file,STRING);
  bindvar("n_fem_nodes",&nfemNodes,INT);
}

int XFEMFrame::exec(const char *name)
{
  if(FEMFrame::exec(name)==0) return 0;

  bindcommand(name,"read_bdy_nodes",read_bdy_nodes(fem_bdy_nodes_file));
  bindcommand(name,"shift_fem_node_id",shift_fem_node_id((int) input[0]));
  bindcommand(name,"read2cn",read2cn());

#ifdef _PARALLEL
  bindcommand(name,"Broadcast_FEM_Param",Broadcast_FEM_Param());
#endif

  return -1;
}

int XFEMFrame::read_bdy_nodes(const char* fname)
{
  INFO_Printf("fname = %s\n", fname);
  INFO_Printf("bdy_file = %s\n", fem_bdy_nodes_file);
  //FILE *bdyfile = fopen(fname,"r"); assert(bdyfile != NULL); 
  FILE *bdyfile = fopen(fem_bdy_nodes_file,"r"); assert(bdyfile != NULL); 
  
  fscanf(bdyfile,"%d %d", &n_bdy_nodes,&n_bdy_tags); 
  assert(n_bdy_nodes > 0 && n_bdy_tags > 0); 
  
  Realloc(bNds,int,n_bdy_nodes);
  Realloc(bTags,int,n_bdy_nodes);
  Realloc(bXs,double,n_bdy_nodes);
  Realloc(bYs,double,n_bdy_nodes);
  Realloc(bZs,double,n_bdy_nodes);

  bindvar("n_bdy_nodes",&n_bdy_nodes,INT);
  bindvar("n_bdy_tags",&n_bdy_tags,INT);
  bindvar("bNds",bNds,INT);
  bindvar("bTags",bTags,INT);
  bindvar("bXs", bXs,DOUBLE);
  bindvar("bYs", bYs,DOUBLE);
  bindvar("bZs", bZs,DOUBLE);

  for(int ibn = 0; ibn < n_bdy_nodes; ibn++) {
    fscanf(bdyfile,"%d %lf %lf %lf %d", &bNds[ibn], &bXs[ibn], &bYs[ibn], &bZs[ibn], &bTags[ibn]);
    //bNds[ibn] += n_existed_atoms;
  }

  printf("read_bdy_nodes finished, n_bdy_nodes = %d\n",n_bdy_nodes);

  return 0;
}


void XFEMFrame::shift_fem_node_id(int np0)
{
  for(int iele=0;iele<_NELE;iele++) {
    for(int j=0;j<_NNODE_PER_ELEMENT;j++) {
      elements[iele*_NNODE_PER_ELEMENT+j] += np0;
    }
  }

  for(int ibn = 0; ibn < n_bdy_nodes; ibn++) {
    bNds[ibn] += np0;
  }
}




#ifdef _PARALLEL
void XFEMFrame::Broadcast_FEM_Param()
{ 
  int n_int_param, n_double_param, i, j, k, m, iA, iN, ind, start, size,
    *buf_int;
  double *buf_double;
      
  if(myDomain==0) { 
    Master_to_Slave("Broadcast_FEM_Param");    
    //            *bNds+*bTags      *elemenets       
    n_int_param = n_bdy_nodes*2 + _NELE*_NNODE_PER_ELEMENT 
      //*map23  
      +_NDIM
      //_NELE, _NNODE_PER_ELEMENT, _NINT_PER_ELEMENT, _NDIM
      + 4 
      //n_bdy_nodes, n_bdy_tags, n_existed_atoms, nfemNodes
      + 4;
      n_double_param = n_bdy_nodes*3 + _NINT_PER_ELEMENT + _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT;
      
      buf_int = (int *) malloc(n_int_param*sizeof(int));
      buf_double = (double *) malloc(n_double_param*sizeof(double));

      // //int bufs
      start = 0;
      buf_int[start++] = _NELE;
      buf_int[start++] = _NNODE_PER_ELEMENT;
      buf_int[start++] = _NINT_PER_ELEMENT;
      buf_int[start++] = _NDIM;
      buf_int[start++] = n_bdy_nodes;
      buf_int[start++] = n_bdy_tags;
      buf_int[start++] = n_existed_atoms;
      buf_int[start++] = nfemNodes;

      for(i= 0;i<n_bdy_nodes;i++) 
      	buf_int[start++] = bNds[i];

      for(i= 0;i<n_bdy_nodes;i++) 
      	buf_int[start++] = bTags[i];      

      for(i=0;i<_NELE; i++) 
      	for(j = 0;j<_NNODE_PER_ELEMENT;j++) {
      	  ind = i*_NNODE_PER_ELEMENT+j;
      	  buf_int[start++] = elements[ind];
      }

      for(i = 0;i<_NDIM;i++)
      	buf_int[start++] = map23[i];
            
      //INFO_Printf("myDomain = %d, start = %d, n_int_param = %d\n", myDomain, start, n_int_param);
      assert(start == n_int_param);

      // double bufs
      start = 0;
      for(i=0;i<n_bdy_nodes;i++, start+=3) {
      	buf_double[i] = bXs[i];
      	buf_double[2*n_bdy_nodes+i] = bYs[i];
      	buf_double[3*n_bdy_nodes+i] = bZs[i];
      }

      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) 
      	buf_double[start++] = gauss_weight[iA];      

      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) 
      	for(j =0;j<_NDIM;j++)
      	  for(k=0;k<_NDIM;k++)
      	    for(m=0;m<_NDIM;m++)
      	      for(iN=0;iN<_NNODE_PER_ELEMENT;iN++) {
      		ind = (((iA*_NDIM+j)*_NDIM+k)*_NDIM+m)*_NNODE_PER_ELEMENT+iN;
      		buf_double[start++] = dFdu[ind];
      	      }    

      //INFO_Printf("myDomain = %d, start = %d, n_double_param = %d\n", myDomain, start, n_double_param);
      assert(start == n_double_param);
  }

  /* broadcasting */
  MPI_Bcast(&n_int_param, 1,           MPI_INT,   0,MPI_COMM_WORLD);
  MPI_Bcast(&n_double_param, 1,        MPI_INT,   0,MPI_COMM_WORLD);

  if(myDomain!=0) {
    buf_int = (int *) malloc(n_int_param*sizeof(int));
    buf_double = (double *) malloc(n_double_param*sizeof(double));
  }

  MPI_Bcast(buf_int,   n_int_param,    MPI_INT,   0,MPI_COMM_WORLD);
  MPI_Bcast(buf_double,n_double_param, MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  assert(n_int_param>0 && n_double_param>0);

  if(myDomain!=0)
    {
      /* unpacking parameters
       * the following lines can be generated from the above lines by
       * replace regexp:
       *  \(buf_double\[[0-9]+\]\)[ ]+= \([^;]+\);
       *  \2 = \1;
       */
      //int bufs

      start = 0;
      _NELE = buf_int[start++];
      _NNODE_PER_ELEMENT = buf_int[start++];
      _NINT_PER_ELEMENT = buf_int[start++];
      _NDIM = buf_int[start++];
      n_bdy_nodes = buf_int[start++];
      n_bdy_tags =  buf_int[start++];
      n_existed_atoms = buf_int[start++];
      nfemNodes = buf_int[start++];

      //allocate memeories on cpu 1-ncpu
      assert(n_bdy_nodes > 0 && n_bdy_tags > 0 && _NP>0 && allocmultiple !=0);       
      Realloc(bNds,int,n_bdy_nodes);
      Realloc(bTags,int,n_bdy_nodes);
      Realloc(bXs,double,n_bdy_nodes);
      Realloc(bYs,double,n_bdy_nodes);
      Realloc(bZs,double,n_bdy_nodes);
      size = _NELE*_NNODE_PER_ELEMENT;
      Realloc(elements,int,size);
      size = _NINT_PER_ELEMENT;
      Realloc(gauss_weight,double,size);
      size = _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT;
      Realloc(dFdu,double,size);
      Realloc(map23,int,2);

      for(i= 0;i<n_bdy_nodes;i++) {
	bNds[i]  = buf_int[start++];
      }
      for(i= 0;i<n_bdy_nodes;i++) {
	bTags[i] = buf_int[start++] ;
      }

      for(i=0;i<_NELE; i++) 
	for(j = 0;j<_NNODE_PER_ELEMENT;j++) {
	  ind = i*_NNODE_PER_ELEMENT+j;
	  elements[ind] = buf_int[start++] ;
	}

      for(i = 0;i<_NDIM;i++) {
	map23[i] = buf_int[start++];
      }
      
      //INFO_Printf("myDomain = %d, start = %d, n_int_param = %d\n", myDomain, start, n_int_param);      
      assert(start == n_int_param);

      // double bufs
      start = 0;
      for(i=0;i<n_bdy_nodes;i++, start+=3) {
	bXs[i] = buf_double[i];
	bYs[i] = buf_double[2*n_bdy_nodes+i];
	bZs[i] = buf_double[3*n_bdy_nodes+i];
      }

      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
	gauss_weight[iA] = buf_double[start++];
      }

      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) 
	for(j =0;j<_NDIM;j++)
	  for(k=0;k<_NDIM;k++)
	    for(m=0;m<_NDIM;m++)
	      for(iN=0;iN<_NNODE_PER_ELEMENT;iN++) {
		ind = (((iA*_NDIM+j)*_NDIM+k)*_NDIM+m)*_NNODE_PER_ELEMENT+iN;
		dFdu[ind] = buf_double[start++];
	      }    
      //INFO_Printf("About to quit:: myDomain = %d, start = %d, n_double_param = %d\n", myDomain, start, n_double_param);
      assert(start == n_double_param);    
    }

  size = _NP*allocmultiple;
  if (myDomain != 0) {
    Realloc(_SRref, Vector3,size);
    Realloc(_Rref, Vector3,size);
  }
  MPI_Bcast(_SRref,3*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_Rref, 3*size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (myDomain == 1 || myDomain == 0) {
    for (int i = 0;i<_NP;i++) {
      INFO_Printf("broadcast[%d]; NP=%d; size =%d; i= %d\n",myDomain,_NP,size,i);
      INFO_Printf("broadcast[%d];_SRref = [%f,%f], _Rref = [%f,%f]\n",myDomain, _SRref[i][0],_SRref[i][1],_Rref[i][0], _Rref[i][1]);
      INFO_Printf("broadcast[%d]; H[00]=%f; H[1][1]=%f; H[2][2]=%f\n",myDomain,_H[0][0], _H[1][1], _H[2][2]);
    }
  }
  //if add the following lines, then there are Segmentation fault. why????????????
  // if(buf_double!=NULL)
  // free(buf_double);
  // if(buf_int!=NULL)
  // free(buf_int);
  
  return;
}

#endif//_PARALLEL


#ifdef _TEST

/* Main Program Begins */
class XFEMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST













// #if 0
// int XFEMFrame::read_elements(const char* fname)
// {
//     int iE, j, iN;
//     char *buffer; char *pp, *q;
//     int np, nframe;
//     char extfname[MAXFILENAMELEN];
//     double coeff;
//     int readCharCount, ind;

//     strcpy(extfname,fname);
//     LFile::SubHomeDir(extfname,extfname);
//     LFile::LoadToString(extfname,buffer,0);

//     pp=buffer;

//     q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
//     sscanf(q,"%d",&_NELE);

//     Alloc_Elements();
    
//     for(iE=0;iE<_NELE;iE++)
//     {
//        q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
//        //INFO_Printf("%s\n",q);
//        for(j=0;j<_NNODE_PER_ELEMENT;j++)
//        {
//             sscanf(q,"%d%n",&iN,&readCharCount);
//             q += readCharCount;
//             ind=iE*_NNODE_PER_ELEMENT+j;
//             elements[ind] = iN;// +n_existed_atoms;
//             INFO_Printf("%d %d  ",ind,elements[ind]);    
//        }
//     }
//     //INFO_Printf("\n");

//     Free(buffer);
//     return 0;
// }

// #endif //0
