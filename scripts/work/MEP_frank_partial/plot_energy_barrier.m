figure(1);

load '../../../runs/frankMEP/EPOT_2.dat'
load '../../../runs/frankMEP/EPOT_1.dat'

plot(EPOT_2(:,1))

figure(2);

d_void_epot = 4.9448
tmp = load(strcat('../../../runs/frankMEP/nucleus-', num2str(0),'.dat'));
numvoids_0 = size(tmp,1);  

for i = 1 : length(EPOT_2)
  tmp = load(strcat('../../../runs/frankMEP/nucleus-', num2str(i-1),'.dat'));
  numvoids = size(tmp,1);  
  NumVoids(i) = numvoids;
  display(numvoids);
  EPOT_2_offset(i) = EPOT_2(i)-(numvoids-numvoids_0)*d_void_epot;
end

dlmwrite('../../../runs/frankMEP/numvoids.txt', NumVoids')
plot(EPOT_2_offset);
