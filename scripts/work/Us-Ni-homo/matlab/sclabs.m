function out = sclabs(Z3,Z1,Z2,Narray1,Freq1,Prob1,Narray2,Freq2,Prob2,Narray3,Freq3,Prob3)

exp_array1=Freq1./Prob1;
exp_array2=Freq2./Prob2;
exp_array3=Freq3./Prob3;

exp_array1(find(Freq1==0))=0;
exp_array2(find(Freq2==0))=0;
exp_array3(find(Freq3==0))=0;

M1=sum(Freq1);
M2=sum(Freq2);
M3=sum(Freq3);
temp=0;

   
for i=1:length(Narray2)
    ind1=find(Narray1 == Narray2(i));
    [Eind1,zzz]=size(ind1);
    if(Z2==1)
        Eind1=0;
    end
    ind2=i;
    ind3=find(Narray3 == Narray2(i));
    [Eind3,zzz]=size(ind3);
    
%    Narray2(i)
%    temp
%Narray1
%Narray2(i)
%Eind1,Eind3,ind1,Freq1(ind1),Freq2(ind2),Freq3(ind3)
    if (Eind1>0 && Eind3>0 && (Freq1(ind1)+Freq2(ind2)+Freq3(ind3)) >0 )
    temp=temp+exp_array2(ind2)*(Freq1(ind1)+Freq2(ind2)+Freq3(ind3))/...
        ((exp_array1(ind1))*M1/Z1+(exp_array2(ind2))*M2/Z2+...
        (exp_array3(ind3))*M3/Z3);
    elseif (Eind1>0 && Eind3==0 && (Freq1(ind1)+Freq2(ind2)) >0)
    temp=temp+exp_array2(ind2)*(Freq1(ind1)+Freq2(ind2))/...
        ((exp_array1(ind1))*M1/Z1+(exp_array2(ind2))*M2/Z2);
    elseif (Eind1==0 && Eind3>0 && (Freq2(ind2)+Freq3(ind3)) >0 )
    temp=temp+exp_array2(ind2)*(Freq2(ind2)+Freq3(ind3))/...
        ((exp_array2(ind2))*M2/Z2+(exp_array3(ind3))*M3/Z3);
    elseif (Freq2(i)>0)
    temp=temp+exp_array2(ind2)*(Freq2(ind2))/...
        ((exp_array2(ind2))*M2/Z2);
    end
    
end
out = abs(Z2-temp);