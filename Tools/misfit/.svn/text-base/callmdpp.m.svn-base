%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File: callmdpp.m
%%
%% Last Modified : Fri Oct  7 22:09:00 2005
%% Wei Cai, caiwei@stanford.edu
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%call MD++ to compute Edata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write MD++ script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp = fopen(sprintf('../../scripts/%s',scriptfile),'w');

fprintf(fp,'# -*-shell-script-*- \n');
fprintf(fp,'# Calculate misfit energy \n');
%
fprintf(fp,'#setnolog \n');
fprintf(fp,'setoverwrite \n');
fprintf(fp,'dirname = ~/Codes/MD++/runs/%s \n',rundir);
fprintf(fp,'#-------------------------------------------- \n');
fprintf(fp,'#Read in potential file \n');
fprintf(fp,'# \n');
if(strcmp(executable,'fs_gpp'))
 fprintf(fp,'potfile = ~/Codes/MD++/potentials/mo_pot readpot \n');
end 
fprintf(fp,'#-------------------------------------------- \n');
fprintf(fp,'#Create Perfect Lattice Configuration \n');
fprintf(fp,'# \n');
fprintf(fp,'crystalstructure = %s latticeconst = %e #(A) \n',crystalstructure,latticeconst);
fprintf(fp,'latticesize = [  %4d %4d %4d %4d  \n',latticesize(1,1),latticesize(1,2),latticesize(1,3),latticesize(1,4));
fprintf(fp,'                 %4d %4d %4d %4d  \n',latticesize(2,1),latticesize(2,2),latticesize(2,3),latticesize(2,4));
fprintf(fp,'                 %4d %4d %4d %4d ]\n',latticesize(3,1),latticesize(3,2),latticesize(3,3),latticesize(3,4));
fprintf(fp,'makecrystal finalcnfile = perf.cn writecn \n');
fprintf(fp,'#------------------------------------------------------------- \n');
%increase maximum number of neighbors per atom
fprintf(fp,'NNM = 100 \n');
fprintf(fp,'#------------------------------------------------------------- \n');
%end
fprintf(fp,'# Calculate misfit energy \n');
fprintf(fp,'input = [ %d  # index of direction normal to misfit plane  1/x 2/y 3/z \n',normaldir);
fprintf(fp,'          %e %e %e   #x0:dx:x1 \n',dxyz(1,1),dxyz(1,2),dxyz(1,3));
fprintf(fp,'          %e %e %e   #y0:dy:y1 \n',dxyz(2,1),dxyz(2,2),dxyz(2,3));
fprintf(fp,'          %e %e %e   #z0:dz:z1 \n',dxyz(3,1),dxyz(3,2),dxyz(3,3));
fprintf(fp,'         ]  # grid number \n');
fprintf(fp,'eval calmisfit \n');
fprintf(fp,'quit \n');
fprintf(fp,'#command = "mv Emisfit.out Gamma110_v2.out " runcommand \n');
fprintf(fp,'#quit \n');
fprintf(fp,'#------------------------------------------------------------- \n');
fprintf(fp,'#Plot Configuration \n');
fprintf(fp,'atomradius = 1.0 bondradius = 0.3 bondlength = 0 \n');
fprintf(fp,'atomcolor = cyan highlightcolor = purple backgroundcolor = gray \n');
fprintf(fp,'bondcolor = red fixatomcolor = yellow \n');
fprintf(fp,'energycolorbar = [ 1 -3.35 -3.32]  highlightcolor = red \n');
fprintf(fp,'plot_atom_info = 3  plotfreq = 10 rotateangles = [ 0 0 0 1 ] \n');
fprintf(fp,'win_width = 600 win_height = 600 \n');
fprintf(fp,'openwin alloccolors rotate saverot refreshnnlist eval plot \n');
fprintf(fp,'sleep quit \n');

fclose(fp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run MD++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system(sprintf('../../bin/%s ../../scripts/%s',executable,scriptfile));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load MD++ output into Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Edata = load(sprintf('../../runs/%s/Emisfit.out',rundir));
Mx=max(Edata(:,1));
My=max(Edata(:,2));
Mz=max(Edata(:,3));

E0=Edata(1,4);
Edata(:,4)=Edata(:,4)-E0;


