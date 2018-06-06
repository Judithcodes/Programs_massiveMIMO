

function [f_k_arr]=get_f_k(BS_density,UE_density,K_Total)

%hold on; 
%   BS_density=BS_density*1e-6;
%    UE_density=UE_density*1e-6;
 q=30;
 %plot(x,nbinpdf(x,q,20/(20+q*4)),'p-');
 %grid on;
 f_k_arr=zeros(1,K_Total);
 sum_total=0;
 sum=nbinpdf(0,q,UE_density/(UE_density+q*BS_density));
 for x=K_Total:1:1000
     sum_total=sum_total+nbinpdf(x,q,UE_density/(UE_density+q*BS_density));
 end
 for idx_k=1:K_Total
     k=idx_k;
%      sum=nbinpdf(0,q,UE_density/(UE_density+q*BS_density));
     if k<=K_Total
          f_k=nbinpdf(k,q,UE_density/(UE_density+q*BS_density))/(1-sum);
     end
     
  
%      if k>K_Total
%          fprintf('k value wrong');
%      end
     f_k_arr(idx_k)=f_k;
 end
 
 f_k_SIM=importdata('f_k_100_300_08_32_.mat');
 
 x=1:32;
 a=x.*f_k_arr;
 k_ave=0;
 aa=0;
  for i=1:32
      k_ave=k_ave+a(i);
      aa=aa+f_k_arr(i);
  end

  f_k_arr(K_Total)=1-aa
  plot(1:32,f_k_arr,'x',1:32,f_k_SIM,'o');
    k_ave
end
