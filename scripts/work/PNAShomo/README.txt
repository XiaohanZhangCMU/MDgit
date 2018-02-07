The following steps will guide you to reproduce all the data concerning the homogeneous
dislocation nucleation in Cu in the PNAS paper.

S. Ryu, K. Kang, W. Cai, "Entropic Effect on Dislocation Nucleation", PNAS, 2011.

Note: All the calculations are performed on cluster "mc2", so that the
executable is eam_mc2.  If you compile MD++ on a different cluster, the
executable name will be different.  You will need to change the names in the
command lines and in the .pbs and .bash scripts accordingly.  You can also
create a symbolic line of your executable to eam_mc2 to avoid changing the
scripts.

Step 1:
   Obtain the zero-K stress-strain curve of Cu under pure shear stress condition.

       make eam build=R SYS=mc2
       mkdir runs/ho-cu-DN
       bin/eam_mc2 scripts/work/PNAShomo/C1_generate_0K_shear_stress_strain_curve.tcl

   This run takes about 5 minutes to complete.
   This run opens a window.  So make sure your X-server is turned on.
   To see the results, do the following.

       cd runs/ho-cu-DN/0K-stressstrain
       octave
         load stress.dat
         load strain.dat
         plot(strain(:,2),-stress(:,4),'o-');


Step 2:
   Prepare zero-K NEB calculation cells.
   To run an interactive job for testing, do the following.

       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.11

   This run will create a folder runs/ho-cu-DN/0K_110eps_112
   where 110eps indicates the shear strain is 0.110 and the slip direction is 112.
   This run takes about 10 seconds to finish.

   To reproduce the data in the paper, run the following jobs.
 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.090 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.092
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.095 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.100 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.105 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.110 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.115 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.120 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.125 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.130 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.135 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.140 
       bin/eam_mc2 scripts/work/PNAShomo/C2_make_0K_pure_shear_stress_cell.tcl 0.145 


Step 3:
   Obtain the equilibrated cell size at finite temperature at a given shear strain. 
   To run an interactive job for testing, do the following.

       bin/eam_mc2 scripts/work/PNAShomo/C3_make_non0K_pure_shear_stress_cell.tcl 300 0.110

   where 300 is temperature (K) and 0.110 is the shear strain.
   This performs one MD simulation (2500 steps) under the NPT ensemble followed by
   25 MD simulations (5000 steps each) under the NVT ensemble with differnt box sizes 
   to equilibrate the stress.  This run takes about 2 hours to finish.

   The stress data during the equilibration run can be ploted by:
       cd runs/ho-cu-DN/300K_110eps_112
       octave
         data=load('stress_trace.dat');
         plot(data(:,1),data(:,2:7));
    (*** need to test ***)

   To reproduce the data in the paper, submit the following job to a parallel cluster.

       mkdir bin1
       cp bin/eam_mc2 bin1

       qsub scripts/work/PNAShomo/PBS/PNAS_Step3a.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step3b.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step3c.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step3d.pbs

   Each job takes about 3 hours to finish.


Step 4:
   Compute minimum energy path (MEP) for all strain and thermal expansions. 
   To run an interactive job for testing, do the following.

       bin/eam_mc2 scripts/work/PNAShomo/C4_MEP_calculation.tcl 200  110  0

   where 200 is the number of relaxation steps (typically we need 2000 for
   convergence), 110 is the shear strain times 1000, 0 is the temperature.
   The energy profile during the MEP relaxation can be monitored by:

       cd runs/ho-cu-DN/0K_120eps_112/NEB
       octave
          load stringeng.out
          plot(stringeng(end:-20:1,3:5:end)','o-')

   To reproduce the data in the paper, submit the following job to a parallel cluster.

       qsub scripts/work/PNAShomo/PBS/PNAS_Step4a.pbs      

   This job performs the MEP search for zero-K cells.  
   So it can run after Step 2 is finished.
   The following jobs can be submitted after Step 3 is finished.

       qsub scripts/work/PNAShomo/PBS/PNAS_Step4b.pbs      
       qsub scripts/work/PNAShomo/PBS/PNAS_Step4c.pbs      
       qsub scripts/work/PNAShomo/PBS/PNAS_Step4d.pbs      
       qsub scripts/work/PNAShomo/PBS/PNAS_Step4e.pbs      

   Each job takes about 5 hours to finish.

Step 5:
   Convert relaxed MEP paths to a series of cn files, to be used for Umbrella Sampling.
   To run an interactive job for testing, do the following.

       bin/eam_mc2 scripts/work/PNAShomo/C5_parse_MEP_path.tcl  110  0

   where 110 is the shear strain times 1000, 0 is the temperature.
   Obviously, this step should only be executed after Step 4 is finished.

   To reproduce the data in the paper, run following shell script.

       scripts/work/PNAShomo/C5_auto.bash

       
Step 6:
   Prepare for Umbrella Sampling.
   To run an interactive job for testing, do the following.

       bin/eam_mc2 scripts/work/PNAShomo/C6_prepare_US.tcl 110 0

   where 110 is the shear strain times 1000, 0 is the temperature.
   This step should only be executed after Step 5 is finished.

   To reproduce the data in the paper, run following shell script.

       scripts/work/PNAShomo/C6_auto.bash

   This step creates many folders such as 

       ho-cu-DN/0K_110eps_112/NEB/run_0
       ho-cu-DN/0K_110eps_112/NEB/run_1
            ...
       ho-cu-DN/0K_110eps_112/NEB/run_11

   These folders contain cn files for initial conditions in Umbrella Sampling simulations. 

   The file ho-cu-DN/0K_110eps_112/NEB/lam_array.dat contains the center
   nucleus size for each Umbrella Sampling window.  Only the beginning lines
   (specified by the maximum # in total#.txt) are important.


Step 7:
   Run Umbrella Sampling.
   This step should only be executed after Step 6 is finished.

   To reproduce the data in the paper, use the following shell scripts to 
   submit jobs to a parallel cluster.

       scripts/work/PNAShomo/C7_auto_qsub_200K.bash
       scripts/work/PNAShomo/C7_auto_qsub_300K.bash
       scripts/work/PNAShomo/C7_auto_qsub_400K.bash
       scripts/work/PNAShomo/C7_auto_qsub_500K.bash

   Each shell script submits Umbrella Sampling jobs for all the strains
   at the given temperature.  At each strain and temperature conditions,
   parallel runs will be performed for up to 24 sampling windows.
   You may not want to run all 4 shell scripts at once because it will
   submit a large number of jobs to your parallel cluster.


Step 8:
   Compute kinetic pre-exponential factor in nucleation rate.

   To run an interactive job for testing, do the following.

       bin/eam_mc2 scripts/PNAShomo/C8_obtain_kinetic_factor.tcl 20 100 120 300

   The value 20 is the window index in US, 100 is the critical nucleus size,
    120 is the shear strain times 1000, 300 is the temperature.
   This will create a folder runs/ho-cu-DN/300K_120eps_112/Kinetic, which
   contains .txt files that will be processed by Matlab.

   This step can be executed 3 hours after Step 7 is started.

   To reproduce the data in the paper, submit the following jobs to a parallel
   cluster.

       qsub scripts/work/PNAShomo/PBS/PNAS_Step8a.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step8b.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step8c.pbs

   After these jobs have finished, run the Matlab file to obtain free energy curves.

       cd runs/ho-cu-DN
       matlab 
        scripts/work/PNAShomo/matlab/USdata_process_RyuCaiOrderParameter.m 

   This matlab file does not work with octave because it uses "fminsearch"
   which is not provided by octave.  You need to transfer the Umbrella Sampling 
   data files to your local computer and run Matlab from there.

      tar -cvf USdata.tar */UMB_*.txt */UMB_*/*.txt */Kinetic*/*.txt */stress_long.dat */strain_long.dat */UMBnobias_*/*.txt

   The other two files, errorbarlogy.m  sclabs.m, need to be in the same folder
   as the USdata_process_RyuCaiOrderParameter.m file.

   Running this Matlab file will produce the free energy profile at the given
   temperature and strain specified at the beginning of this file.

       

Step 9: (optional)
   Compute the stress at each strain and temperature more accurately (than Step 3).

   To reproduce the data in the paper, submit the following jobs to a parallel
   cluster.

       qsub scripts/work/PNAShomo/PBS/PNAS_Step9a.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step9b.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step9c.pbs
       qsub scripts/work/PNAShomo/PBS/PNAS_Step9d.pbs

   Each job takes about 3 hours to finish.

