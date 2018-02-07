/* This file contains the main program that
 *  is included in all MD++ programs
 *
 *  Last modified: Yanming Wang, Mar 13, 2015 for phasefield MPI code
 */
#ifdef _STK_MPI
#include <mpi.h>
#endif


#ifdef _USETCL
int Tcl_AppInit(Tcl_Interp *interp)                                  
{                                                                    
    return sim.Tcl_AppInit(interp);                                  
}                                                                    
                                                                     
int MD_Cmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{                                                                    
    FILE *istream, *ostream; int mypipe[2];                          
                                                                     
    if (pipe (mypipe))                                               
    {                                                                
        fprintf (stderr, "Pipe failed.\n");                          
        return EXIT_FAILURE;                                         
    }                                                                
    istream = fdopen(mypipe[0], "r");                                
    ostream = fdopen(mypipe[1], "w");                                
    for(int i=1;i<argc;i++)                                          
        fprintf(ostream,"%s ",(char *)argv[i]);                      
    fprintf(ostream,"\n");                                           
    fclose(ostream);                                                 
    sim.parse_line(istream);                                         
    fclose(istream);                                                 
    return TCL_OK;                                                   
}                                                                    
                                                                     
int MD_SetVar(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{                                                                    
    int i, curn;                                                     
    curn = sim.identify(argv[1]);                                    
    if(curn<0)                                                       
    {                                                                
        WARNING("Unrecognized variable "<<argv[1]);                  
        return TCL_OK;                                               
    }                                                                
    for(i=2;i<argc;i++)                                              
    {                                                                
        switch(sim.vartype[curn])                                    
        {                                                            
        case(INT): sscanf(argv[i],"%d",(int *)sim.varptr[curn]+i-2+sim.shift); break; 
        case(LONG): sscanf(argv[i],"%ld",(long *)sim.varptr[curn]+i-2+sim.shift); break;
        case(DOUBLE): sscanf(argv[i],"%lf",(double *)sim.varptr[curn]+i-2+sim.shift); break;
        case(STRING): strcpy((char *)sim.varptr[curn]+i-2+sim.shift,argv[i]); break;   
        default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); 
        }                                                            
    }                                                                
    return TCL_OK;                                                   
}                                                                    
                                                                     
int MD_GetVar(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{                                                                    
    int i, curn, offset;  char buf[1000];                            
    curn = sim.identify(argv[1]);                                    
    if(curn<0)                                                       
    {                                                                
        WARNING("Unrecognized variable "<<argv[1]);                  
        return TCL_OK;                                               
    }                                                                
    for(i=1;i<argc;i++)                                              
    {                                                                
        if(i==1) offset=0;                                           
        else sscanf(argv[i],"%d",&offset);                           
        if(i>2) sprintf(buf," ");                                    
        switch(sim.vartype[curn])                                    
        {                                                            
        case(INT): sprintf(buf,"%d",*((int *)sim.varptr[curn]+offset+sim.shift)); break; 
        case(LONG): sprintf(buf,"%ld",*((long *)sim.varptr[curn]+offset+sim.shift)); break;
        case(DOUBLE): sprintf(buf,"%.16g",*((double *)sim.varptr[curn]+offset+sim.shift)); break;
        case(STRING): sprintf(buf,"%s",(char *)sim.varptr[curn]+offset+sim.shift); break;
        default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); 
        }                                                            
    }                                                                
    Tcl_SetResult(interp, buf, TCL_VOLATILE);                        
    return TCL_OK;                                                   
}                                                                    


#ifdef _USEOCTAVE
int OCTAVE_Cmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{   
    
#if 0 // parse_and_execute (FILE *) does not work                                                            
    FILE *istream, *ostream; int mypipe[2];                          
    if (pipe (mypipe))                                               
    {                                                                
        fprintf (stderr, "Pipe failed.\n");                          
        return EXIT_FAILURE;                                         
    }                                                                
    istream = fdopen(mypipe[0], "rt");                                
    ostream = fdopen(mypipe[1], "wt");                                
    fprintf(stdout, "OCTAVE:");
    for(int i=1;i<argc;i++)
    {
        fprintf(stdout," %s", argv[i]); 
        fprintf(ostream,"%s ",(char *)argv[i]);                      
    }
    fprintf(stdout,"\n");
    fprintf(ostream,"\n");                                           
    fclose(ostream);                                                 
    parse_and_execute (istream);
    fclose(istream);
#else
    FILE *ostream; char filename[100], command[100]; int pid;                          
    pid = getpid();
    sprintf(filename,"/tmp/octave_scratch_p%d.m",pid);
    ostream = fopen(filename,"wt");
    fprintf(stdout, "OCTAVE:");
    for(int i=1;i<argc;i++)
    {
        fprintf(stdout," %s", argv[i]); 
        fprintf(ostream," %s ",(char *)argv[i]);                      
    }
    fprintf(stdout,"\n");
    fprintf(ostream,"\n");                                           
    fclose(ostream);                                                 
    parse_and_execute (filename);

    // parse_and_execute (FILE *) does not work
    //FILE *istream;
    //istream = fopen(filename,"rt");
    //parse_and_execute (istream);
    //fclose (istream);

    remove(filename);
#endif
    return TCL_OK;                                                   
}                                                                    
int OCTAVE_Run(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{   
    for(int i=1;i<argc;i++)                                          
    {
        fprintf(stdout, "OCTAVE_Run: %s\n", argv[i]);                          
        parse_and_execute (argv[i]);
    }
    return TCL_OK;                                                   
}                                                                    
int OCTAVE_GetVar(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{   
    /* make Octave write the variable into a data file */
    FILE *ostream; char filename[100], command[100]; int pid;                          
    pid = getpid();
    sprintf(filename,"/tmp/octave_scratch_p%d.m",pid);
    ostream = fopen(filename,"wt");
    fprintf(ostream,"fid = fopen('/tmp/octave_variable_p%d.out','w');\n",pid);
    if(argc <= 1)
    {
      fprintf(ostream,"fprintf(fid,'%%.16g\\n',%s);\n",(char *)argv[1]);
    } 
    else if(argc > 2)
    {
      if (strcmp(argv[2],"int")==0)
        fprintf(ostream,"fprintf(fid,'%%d\\n',%s);\n",(char *)argv[1]);
      else if (strcmp(argv[2],"double")==0)
        fprintf(ostream,"fprintf(fid,'%%.16g\\n',%s);\n",(char *)argv[1]);
      else if (strcmp(argv[2],"string")==0)
        fprintf(ostream,"fprintf(fid,'%%s\\n',%s);\n",(char *)argv[1]);
      else
        FATAL("unknown vartype ("<<argv[2]<<")"); 
    }                                                            
    fprintf(ostream,"fprintf(fid,'%%.16g\\n',%s);\n",(char *)argv[1]);
    fprintf(ostream,"fclose(fid);\n");
    fprintf(ostream,"\n");                                           
    fclose(ostream);                                                 
    parse_and_execute (filename);
    remove(filename);

    /* read data file into buf and give it to Tcl */
    FILE *istream; char buf[1000];
    sprintf(filename,"/tmp/octave_variable_p%d.out",pid);
    istream = fopen(filename,"rt");
    fgets(buf,1000,istream); 
    if(buf[strlen(buf)-1]=='\n') buf[strlen(buf)-1]=0;
    fclose(istream);
    remove(filename);

    Tcl_SetResult(interp, buf, TCL_VOLATILE);                        
    return TCL_OK;                                                   
}                                                                    
#endif//_USEOCTAVE

#endif//_USETCL

#ifdef _USEPY

/* execute MD++ commands */
static PyObject*
md_cmd(PyObject *self, PyObject *args)
{
    FILE *istream, *ostream; int mypipe[2];                          
    char *s; int ret_val;
    if(!PyArg_ParseTuple(args, "s:cmd", &s))
        return NULL;
                                                                     
    if (pipe (mypipe))                                               
    {                                                                
        fprintf (stderr, "Pipe failed.\n");                          
        return NULL;                                         
    }                                                                
    istream = fdopen(mypipe[0], "r");                                
    ostream = fdopen(mypipe[1], "w");                                
    fprintf(ostream,"%s ",s);                      
    fprintf(ostream,"\n");                                           
    fclose(ostream);                                                 
    ret_val = sim.parse_line(istream);                                         
    fclose(istream);                                                 
    return Py_BuildValue("i", ret_val);
}

/* get MD++ variable values */
static PyObject*
md_get(PyObject *self, PyObject *args)
{                       
    int i, curn, offset;  char buf[1000];                            
    char *s; int ok;

    offset = 0;
    ok = PyArg_ParseTuple(args, "s|i:get", &s, &offset);

    if (!ok) return NULL;

    curn = sim.identify(s);                                    
    if(curn<0)                                                       
    {                                                                
        WARNING("Unrecognized variable "<<s);                  
        return NULL;
    }                                                                

    switch(sim.vartype[curn])                                    
    {                                                            
    case(INT): return Py_BuildValue("i", *((int *)sim.varptr[curn]+offset+sim.shift)); break; 
    case(LONG): return Py_BuildValue("l", *((long *)sim.varptr[curn]+offset+sim.shift)); break; 
    case(DOUBLE): return Py_BuildValue("d", *((double *)sim.varptr[curn]+offset+sim.shift)); break; 
    case(STRING): return Py_BuildValue("s", ((char *)sim.varptr[curn]+offset+sim.shift)); break; 
    default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); return NULL;
    }                                                            
    return NULL;                                                   
}                                                                    

static PyObject*
md_get_array(PyObject *self, PyObject *args)
{                       
    int curn, istart, iend, iskip, offset;  char buf[1000];                            
    char *s; int ok;
    PyObject *lst, *num;  Py_ssize_t i;

    istart = 0; iend = 10; iskip = 1;
    ok = PyArg_ParseTuple(args, "s|iii:get_array", &s, &istart, &iend, &iskip);

    if (!ok) return NULL;

    curn = sim.identify(s);                                    
    if(curn<0)                                                       
    {                                                                
        WARNING("Unrecognized variable "<<s);                  
        return NULL;
    }                                                                

    lst = PyList_New((iend-istart)/iskip);
    if (!lst) return NULL;

    i = 0;
    for (offset=istart; offset<iend; offset+=iskip)
    {
        switch(sim.vartype[curn])                                    
        {                                                            
        case(INT): num = Py_BuildValue("i", *((int *)sim.varptr[curn]+offset+sim.shift)); break; 
        case(LONG): num = Py_BuildValue("l", *((long *)sim.varptr[curn]+offset+sim.shift)); break; 
        case(DOUBLE): num = Py_BuildValue("d", *((double *)sim.varptr[curn]+offset+sim.shift)); break; 
        case(STRING): num = Py_BuildValue("s", ((char *)sim.varptr[curn]+offset+sim.shift)); break; 
        default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); return NULL;
        }                                                            
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);
        i++;
   }
   return lst;
}                                                                    

static PyObject*
md_set_array(PyObject *self, PyObject *args)
{                       
    int curn, istart, iskip, offset, numLines;  char buf[1000];                            
    char *s; int ok;
    PyObject *lst, *num;  Py_ssize_t i;

    istart = 0; iskip = 1;
    ok = PyArg_ParseTuple(args, "O!s|ii:set_array", &PyList_Type, &lst, 
                            &s, &istart, &iskip);

    if (!ok) return NULL;

    numLines = PyList_Size(lst);
    if (numLines < 0) return NULL;

    curn = sim.identify(s);                                    
    if(curn<0)                                                       
    {                                                                
        WARNING("Unrecognized variable "<<s);                  
        return NULL;
    }                                                                

    offset = istart;
    for (i = 0; i<numLines; i++)
    {
        num = PyList_GetItem(lst, i);
        switch(sim.vartype[curn])                                    
        {                                                            
        case(INT): *((int *)sim.varptr[curn]+offset+sim.shift) = (int) PyInt_AsLong(num); break; 
        case(LONG): *((long *)sim.varptr[curn]+offset+sim.shift) = PyInt_AsLong(num); break; 
        case(DOUBLE): *((double *)sim.varptr[curn]+offset+sim.shift) = PyFloat_AsDouble(num); break; 
        case(STRING): strcpy((char *)sim.varptr[curn]+i-2+sim.shift,PyString_AsString(num)); break;   
        default: FATAL("unknown vartype ("<<sim.vartype[curn]<<")"); return NULL;
        }                                                            
        offset += iskip; 
   }
   return Py_None;
}                                                                    

static PyMethodDef MDMethods[] = {
    {"cmd", md_cmd, METH_VARARGS,
     "Execute string as MD++ commands."},
    {"get", md_get, METH_VARARGS,
     "Get MD++ variable values."},
    {"get_array", md_get_array, METH_VARARGS,
     "Get MD++ variable array values."},
    {"set_array", md_set_array, METH_VARARGS,
     "Set MD++ variable array values."},
    {NULL, NULL, 0, NULL}
};
#endif//_USEPY

#ifdef _PARALLEL
int main_master(int argc, char *argv[]);
int main_slave(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    sim.ParallelInit(&argc,&argv);

    if(sim.myDomain==0)
        main_master(argc,argv);
    else
        main_slave(argc,argv);

    MPI_Finalize();
    return 0;
}
int main_slave(int argc, char *argv[])
{
    sim.initvars();
    strcpy(sim.myname,argv[0]);
    
    sim.WaitForCommand();
    return 0;
}
#else
#define main_master main
#endif //_PARALLEL


/* This is the main program for single processor run */
int main_master(int argc, char *argv[])
{
#ifdef _STK_MPI
    MPI_Init(&argc, &argv);
#endif
    
#ifdef _USETCL
#ifdef _USEOCTAVE
    /* set up a copy of argv for octave */
    int argcp;
    char **argvp, argv0[100];
    strcpy(argv0,argv[0]);
    argvp=(char **)malloc(sizeof(char *)*1);
    argvp[0] = argv0;
    argcp = 1;

    for(int i=0;i<argc;i++)
       fprintf(stdout,"argv[%d]=%s\n",i,argv[i]);

    int embedded;
    octave_main(argcp,argvp,embedded=1);
#endif//_USEOCTAVE
#endif//_USETCL

#ifdef _USEPY
  Py_SetProgramName(argv[0]);  /* optional but recommended */
  Py_Initialize();
  if (argc>1) PySys_SetArgvEx(argc-1,argv+1,0); else PySys_SetArgvEx(0,NULL,0);
  Py_InitModule("mdpp", MDMethods);
#endif//_USEPY

    sim.initvars();
    strcpy(sim.myname,argv[0]);

    for(int i=0;i<argc;i++)
       fprintf(stdout,"argv[%d]=%s\n",i,argv[i]);

#ifndef _USETCL
#ifdef _USEPY
    PyRun_SimpleString("from time import time,ctime\n"
                       "print 'Start time is',ctime(time())\n");
    if (argc>1) 
    {
       PyRun_SimpleFile(SCParser::getfilehandler(argc,argv),argv[1]);
    }
    PyRun_SimpleString("from time import time,ctime\n"
                       "print 'End time is',ctime(time())\n");
    Py_Finalize();
#else
    if (argc>2) 
    {
       sim.set_dirname(argv[2]);
    }
    if (argc>3) sim.set_randseed(argv[3]);
    sim.parse(SCParser::getfilehandler(argc,argv));
#endif//_USEPY
#else
    Tcl_Main_Wrapper(argc, argv);
#endif//_USETCL

#ifdef _STK_MPI
    MPI_Finalize();
#endif
    
    return 0;    
}

