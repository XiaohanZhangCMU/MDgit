close all;clear;
kB=8.617343*10^(-5);
Nsite = 832;Ntot=14976;

%%%
nobiasfor0=0;
nobiasfor1=0;

%%% input the field and kBT data
eps=125;
T=300;

%%% cut off for the kinetic factor
fccut=500;

%%% cut off of the US windows side
cut=2;

%%% initial setup
FILE_type='.txt';
eps_val=sprintf('%03deps_112/',eps);
T_val=sprintf('%3dK_',T);
UMB_val=sprintf('UMB_');

% get ncpu
for i=0:100
    CHECK=strcat(T_val,eps_val,UMB_val,num2str(i),'/Freq.txt');
    key1=exist(CHECK);
    if (key1<1)
        ncpu=i;
        break;
    end
end
Z=zeros(1,ncpu);
Z(1)=1;

disp(sprintf('\nMerging Umbrella Sampling (US) data'));
disp('0. We are getting number of windows');
% comput the probability by overlapping window method
data=load(strcat(T_val,eps_val,UMB_val,'arg1',FILE_type));
windowsize=data(2)*2+1-2*cut;
Narray_TOT=zeros(windowsize,ncpu);
Freq_TOT=zeros(windowsize,ncpu);
Prob_TOT=zeros(windowsize,ncpu);
EXP_TOT=zeros(windowsize,ncpu);

disp('1. We are getting Bias Function for Each Windows');
for i=1:ncpu
    %disp(sprintf('%d',i));
    Narraytemp=load(strcat(T_val,eps_val,UMB_val,num2str(i-1),'/Narray',FILE_type));
    Freqtemp=load(strcat(T_val,eps_val,UMB_val,num2str(i-1),'/Freq',FILE_type));
    Probtemp=load(strcat(T_val,eps_val,UMB_val,num2str(i-1),'/Prob',FILE_type));
    
    if (i==1 && nobiasfor0==1)
        Narraytemp=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Narray',FILE_type));
        Freqtemp=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Freq',FILE_type));
        Probtemp=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Prob',FILE_type));
        disp('check0!')
    end
    
    if (i==2 && nobiasfor1==1)
        Narraytemp=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Narray',FILE_type));
        Freqtemp=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Freq',FILE_type));
        Probtemp=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Prob',FILE_type));
        disp('check1!')
    end
    %i
    %Narraytemp
    Narraytemp=Narraytemp(1+cut:end-cut);
    Freqtemp=Freqtemp(1+cut:end-cut);
    Probtemp=Probtemp(1+cut:end-cut);
    
    Narray_TOT(:,i)=Narraytemp;
    Freq_TOT(:,i)=Freqtemp;
    Prob_TOT(:,i)=Probtemp;
    
    exp_array=Freqtemp./Probtemp;
    exp_array(find(Freqtemp==0))=0;
    EXP_TOT(:,i)=exp_array;
end

disp('2. We are getting Z(i)s for Each Windows');
for i=2:ncpu
    %disp(sprintf('%d',i));
    Narray1=load(strcat(T_val,eps_val,UMB_val,num2str(max(i-3,0)),'/Narray',FILE_type));
    Narray2=load(strcat(T_val,eps_val,UMB_val,num2str(i-2),'/Narray',FILE_type));
    Narray3=load(strcat(T_val,eps_val,UMB_val,num2str(i-1),'/Narray',FILE_type));
    Freq1=load(strcat(T_val,eps_val,UMB_val,num2str(max(i-3,0)),'/Freq',FILE_type));
    Freq2=load(strcat(T_val,eps_val,UMB_val,num2str(i-2),'/Freq',FILE_type));
    Freq3=load(strcat(T_val,eps_val,UMB_val,num2str(i-1),'/Freq',FILE_type));
    Prob1=load(strcat(T_val,eps_val,UMB_val,num2str(max(i-3,0)),'/Prob',FILE_type));
    Prob2=load(strcat(T_val,eps_val,UMB_val,num2str(i-2),'/Prob',FILE_type));
    Prob3=load(strcat(T_val,eps_val,UMB_val,num2str(i-1),'/Prob',FILE_type));
    
    if (max(i-3,0)==0 && nobiasfor0==1)
        Narray1=load(strcat(T_val,eps_val,'UMBnobias_',num2str(max(i-3,0)),'/Narray',FILE_type));
        Freq1=load(strcat(T_val,eps_val,'UMBnobias_',num2str(max(i-3,0)),'/Freq',FILE_type));
        Prob1=load(strcat(T_val,eps_val,'UMBnobias_',num2str(max(i-3,0)),'/Prob',FILE_type));
    end
    if (i-2==0 && nobiasfor0==1)
        Narray2=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-2),'/Narray',FILE_type));
        Freq2=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-2),'/Freq',FILE_type));
        Prob2=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-2),'/Prob',FILE_type));
    end
    if (i-1==0 && nobiasfor0==1)
        Narray3=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Narray',FILE_type));
        Freq3=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Freq',FILE_type));
        Prob3=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Prob',FILE_type));
    end
    
    if (max(i-3,0)==1 && nobiasfor1==1)
        Narray1=load(strcat(T_val,eps_val,'UMBnobias_',num2str(max(i-3,0)),'/Narray',FILE_type));
        Freq1=load(strcat(T_val,eps_val,'UMBnobias_',num2str(max(i-3,0)),'/Freq',FILE_type));
        Prob1=load(strcat(T_val,eps_val,'UMBnobias_',num2str(max(i-3,0)),'/Prob',FILE_type));
    end
    if (i-2==1 && nobiasfor1==1)
        Narray2=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-2),'/Narray',FILE_type));
        Freq2=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-2),'/Freq',FILE_type));
        Prob2=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-2),'/Prob',FILE_type));
    end
    if (i-1==1 && nobiasfor1==1)
        Narray3=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Narray',FILE_type));
        Freq3=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Freq',FILE_type));
        Prob3=load(strcat(T_val,eps_val,'UMBnobias_',num2str(i-1),'/Prob',FILE_type));
    end
    
    
    
    Narray1=Narray1(1+cut:end-cut);
    Narray2=Narray2(1+cut:end-cut);
    Narray3=Narray3(1+cut:end-cut);
    Freq1=Freq1(1+cut:end-cut);
    Freq2=Freq2(1+cut:end-cut);
    Freq3=Freq3(1+cut:end-cut);
    Prob1=Prob1(1+cut:end-cut);
    Prob2=Prob2(1+cut:end-cut);
    Prob3=Prob3(1+cut:end-cut);
    Z(i)=Z(i-1);  
   [x,fval] = fminsearch(@sclabs,Z(i),optimset('TolX',Z(i-1)*1e-14),...
       Z(max(1,i-2)),Z(i-1),Narray1,Freq1,Prob1,Narray2,...
            Freq2,Prob2,Narray3,Freq3,Prob3);
    Z(i)=x;
end

Narraymin=load(strcat(T_val,eps_val,UMB_val,num2str(0),'/Narray',FILE_type));
Narraymax=load(strcat(T_val,eps_val,UMB_val,num2str(ncpu-1),'/Narray',FILE_type));
n_min=min(Narraymin(1+cut:end-cut));
n_max=max(Narraymax(1+cut:end-cut));

[aa,bb]=size(Narray_TOT);
p=zeros(1,n_max-n_min+1);
s=0;

disp('3. We are computing probability distribution using Z(i)s');
for l=n_min:n_max
%    for l=20:20
    ind=find (Narray_TOT == l);
    nume=sum(Freq_TOT(ind));
    deno=sum(  EXP_TOT(ind)'.*...
        sum(Freq_TOT(:,ceil(ind/aa)))./Z(ceil(ind/aa)) );
    s=s+1;
    if(nume==0)
       p(s)=0;
    else
       p(s)=nume/deno;
    end
    perr(s)=sqrt( p(s)*1/(deno));
end

n_seq=n_min:n_max;
% Change variable such that n=0 is n(1).
indn0=find(n_seq==0);
nmax=n_max;
n_array=0:1:nmax;
p_array=zeros(1,nmax+1);
perr_array=zeros(1,nmax+1);
for i=1:nmax+1
    p_array(i)=p(indn0+i-1);
    perr_array(i)=perr(indn0+i-1);
end

disp('4. We are drawing raw data');
figure(1)
errorbar(n_array,p_array,perr_array,'.-');
errorbarlogy;
xlabel('n');
ylabel('Probability');
xlim([0 nmax])
disp(sprintf('       Now, sum of raw p=%f with cut=%d',sum(p_array),cut));
perr_r=perr_array./p_array;
perr_r=perr_r(find(perr_r>0));
disp(sprintf('       max(perr/p_array)=%f',max(perr_r)));

disp('5-1. Now processing by "sum(p) = 1 normalization scheme"');
p_type1=p_array/sum(p_array);
perr_type1=perr/sum(p_array);
n_type1nonzero=n_array(find(p_type1>0));
p_type1nonzero=p_type1(find(p_type1>0));
[p_type1c,type1ind]=min(p_type1nonzero);
nc=n_type1nonzero(type1ind);
F_type1c=-kB*T*log(p_type1c/Nsite);
%perr_type1nonzero=perr_type1(find(p_type1>0)); % do this later
%Ferr_type1c=kB*T*Nsite*perr_type1nonzero(type1ind);
ind1=round(type1ind-0.18*length(n_type1nonzero));
ind2=min(round(type1ind+0.18*length(n_type1nonzero)),length(n_type1nonzero));
nblock=n_type1nonzero(ind1:ind2);
Fblock1=-kB*T*log(p_type1nonzero(ind1:ind2)/Nsite);
coef1=polyfit(nblock,Fblock1,2);F1fit=polyval(coef1,nblock);
n1cfit=-coef1(2)/coef1(1)/2;
F1cfit=max(F1fit);
p1cfit=Nsite*exp(-F1cfit/kB/T);
disp(sprintf('        [Here, F1c=-kB T ln(p1c/Nsite), as p1c is for Nsite]'));
disp(sprintf('   a. From simple min\n       nc=%d, p1c=%e Fc=%f (eV)',nc,p_type1c,F_type1c));
disp(sprintf('   b. From fit\n       n1cfit=%f, p1cfit=%e, F1cfit=%f (eV)',n1cfit,p1cfit,F1cfit));

%sum(p_array)
%n_array(1:5)
%break;
disp('5-2. Now processing by "sum_n{n*p(n)} = 1 normalization scheme"');
%p_type2=p_array/( p_array(1)+p_array(1:end)*n_array(1:end)' );
% originally, the above is right expression, however, our MD++ already
% take that into account in C++ source code, so we can simply do
% this simple normalization.
p_type2=p_array/sum(p_array);
p_type2nonzero=p_type2(find(p_type2>0));
n_type2nonzero=n_array(find(p_type2>0));
[p_type2c,type2ind]=min(p_type2nonzero);
F_type2=-kB*T*log(p_type2);
F_type2c=-kB*T*log(p_type2c);
Fblock2=-kB*T*log(p_type2nonzero(ind1:ind2));
coef2=polyfit(nblock,Fblock2,2);F2fit=polyval(coef2,nblock);
n2cfit=-coef2(2)/coef2(1)/2;
F2cfit=max(F2fit);
p2cfit=exp(-F2cfit/kB/T);
disp(sprintf('        [Here F2c=-kB T ln(p2c), as p2c is for a single site]'));
disp(sprintf('   a. From simple min  nc=%d, p2c=%e F2c=%f (eV)',nc,p_type2c,F_type2c));
disp(sprintf('   b. From fit\n       n2cfit=%f, p2c=%e, F2cfit=%f (eV)',n2cfit,p2cfit,F2cfit));
figure(2)
plot(n_array,F_type2,'.-',nblock,F2fit,'ko');
xlabel('n');
ylabel('Free Energy (eV)');
title('Valid only for 5-2 normalization scheme')
xlim([0 nmax])

disp('6. We are getting fc^+ value');
% computation of Zeldovich factor
Zeldo = (-2*coef2(1)/2/pi/kB/T)^0.5;
disp(sprintf('       Zeldovich Factor = %f',Zeldo));
%%% obtain Kinetic factor
n_file=strcat(T_val,eps_val,'Kinetic/n_array.txt');
time_file=strcat(T_val,eps_val,'Kinetic/time_array.txt');
key2=exist(n_file);
if (key2==0)
    disp('Kinetic factor data is missing, run the simulation!');
    break;
end

n_arrayK=load(n_file);
indn=find(n_arrayK(:,2)<50);
n_arrayK=n_arrayK(indn,:);
time_array=load(time_file)*10^(-12);
dnsqr=mean(n_arrayK(1:fccut,:).^2);
errn=std(n_arrayK(1:fccut,:).^2);

k=fit_line(time_array(1:4),dnsqr(1:4),[10^14]);
fcp=k/2;

disp(sprintf('       fc+ = %e',fcp));
disp(sprintf('       fc+ * zeldo = %e',fcp*Zeldo));
figure(5)
subplot(121)
plot(time_array(1:10),dnsqr(1:10),'o-',time_array(1:4),k*time_array(1:4),'k-')
xlabel('time(s)');
ylabel('<dN^2>');
subplot(122)
plot(time_array(1:10),n_arrayK(1:fccut,1:10)')
xlabel('time(s)');
ylabel('dN');

J_new=Ntot*fcp*Zeldo*exp(-F1cfit/(kB*T));
J_old=Ntot*fcp*Zeldo*exp(-F2cfit/(kB*T));
disp(sprintf('J_1=%e, J_2=%e',J_new,J_old));

disp(sprintf('\n==========================')); 
disp(sprintf('Publish Following Data Set')); 
disp(sprintf('==========================')); 
disp(sprintf('       Zeldovich Factor = %f',Zeldo));
disp(sprintf('       fc+ = %e',fcp));
disp(sprintf('       fc+ * Zeldo = %e',fcp*Zeldo));
if (  (nobiasfor0+nobiasfor1) > 0 )
   disp(sprintf('    nc=%f, pc=%e, Fc=%f (eV)',n2cfit,p2cfit,F2cfit));   
   disp(sprintf('J = Ntot*fcp*Zeldo*exp(-Fc/(kB*T)) = %e where Ntot=%d',J_old,Ntot));
else
   disp(sprintf('    nc=%f, pc=%e, Fc=%f (eV)',n1cfit,p1cfit,F1cfit)); 
   disp(sprintf('J = Ntot*fcp*Zeldo*exp(-Fc/(kB*T)) = %e where Ntot=%d',J_new,Ntot));
end

indMAX=(find(n_array==max(Narray_TOT(:,3))));
if 1
    disp(sprintf('==========================')); 
    disp(sprintf('Dubug of Window 1,2,3')); 
    disp(sprintf('==========================')); 
    disp(sprintf('open up DUBUG window at figure(1000)'));
    figure(1000)
    semilogy(Narray_TOT(:,1),Prob_TOT(:,1),'ro-',...
        Narray_TOT(:,2),Prob_TOT(:,2),'bo-',...
        Narray_TOT(:,3),Prob_TOT(:,3),'mo-',...
        n_array(1:indMAX),p_array(1:indMAX),'k.-');
    Pabs=Freq_TOT(:,1)/(sum(Freq_TOT(:,1))/Nsite)/Nsite;
    disp(sprintf('%13d| ',Narray_TOT(:,1)));
    disp(sprintf('%13d| ',Freq_TOT(:,1)));
    disp(sprintf('%10.6e| ',Pabs*Nsite));
    disp(sprintf('%10.6e| ',Pabs));
    disp(sprintf('%10.6e| ',p_type1(1:7)));
    
    %disp(sprintf('\n'));
    %disp(sprintf('%10.6e| ',Prob_TOT(:,1)./sum(Freq_TOT(:,1))/Nsite));
end
