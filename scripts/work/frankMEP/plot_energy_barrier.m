
%close all;
clear;
clc;

dirname = '../../../runs/frankMEP/Mar28_0_05/';
epot2_file = strcat(dirname,'/EPOT_2.dat');
epot1_file = strcat(dirname,'/EPOT_1.dat');
EPOT_2=load(epot2_file);
EPOT_1=load(epot1_file);



d_void_epot = 4.63 ; % 2.64; % 4.63; %9.04101350863513;4.63
tmp = load(strcat(dirname,'/nucleus-', num2str(0),'.dat'));
numvoids_0 = size(tmp,1);  

%these are the indices for strained_.cn that have perfect frank-partial dislocations
%Mar28_0_05
valid_disl_index_0 = [ 0 1 2 3 4 5 6 7 8 11 16 18 ]+1;
%Mar28_0_04
%valid_disl_index_0 = [ 0 1 2 3 4 5 6 7 8 11 14 18 ]+1;
%valid_disl_index_0 = [1:length(EPOT_2)];

EPOT_2_0 = EPOT_2(valid_disl_index_0);

for i = 1 : length(EPOT_2_0)
  id = valid_disl_index_0(i);
  tmp = load(strcat(dirname,'/nucleus-', num2str(id-1),'.dat'));
  numvoids = size(tmp,1);  
  NumVoids_0(i) = numvoids;
  EPOT_2_offset_0(i) = EPOT_2_0(i)-(numvoids)*d_void_epot;
  display(EPOT_2_0(i));
  display(numvoids);
end

figure(2);
hold on;
plot(NumVoids_0, EPOT_2_offset_0(1,:)'-EPOT_2_offset_0(1,1),'*-');
%plot(NumVoids_0, EPOT_2_offset_0(1,:)','*-');
