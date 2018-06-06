function []=main_DL_mMIMO(BS_Density,UE_Density)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: M_mimo_anal_LANL_limUE_UL.m
% Authors: Xuefeng Yao
% Version: 8.0
% Date: 15/6/2017
% Description: Analytical Results for the Massive MIMO UP LINK
% Condersidering Los and NLos with Limited UE in each BS cell
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modifications in this version from V1: using the distance between typical and typical Bs instead of R_E 
BS_Density_str=num2str(BS_Density);
UE_Density_str=num2str(BS_Density);
E_betaLL_name=['E_beta_',BS_Density_str,'_',UE_Density_str,'_LoS_LoS.mat'];
E_betaLNL_name=['E_beta_',BS_Density_str,'_',UE_Density_str,'_LoS_NLoS.mat'];
E_betaNLL_name=['E_beta_',BS_Density_str,'_',UE_Density_str,'_NL0S_LoS.mat'];
E_betaNLNL_name=['E_beta_',BS_Density_str,'_',UE_Density_str,'_NL0S_NLoS.mat'];
f_k_name=['ANAF_K_',BS_Density_str,'_',UE_Density_str,'.mat'];
ccdf_name=['ANAcov_pro_0',BS_Density_str,'_',UE_Density_str,'.mat'];
ccdf_threshold=['ANAccdf',BS_Density_str,'_',UE_Density_str,'.mat'];
alpha_LoS = 2.42;
alpha_NLoS = 4.28;
A_LoS_dB_Macro=30.8;
A_NLoS_dB_Macro=2.7;
A_LoS = 10^(-1*A_LoS_dB_Macro/10);
A_NLoS = 10^(-1*A_NLoS_dB_Macro/10);
d0=300;

BS_density=round(BS_Density)/1000000;
BS_density_km=BS_density*1000000;
len_BS_density=length(BS_density);
UE_density=round(UE_Density)/1000000;
UE_density_km=UE_density*1000000;
e=0.8;
K_total=32;
M=64;
Nt_noise=zeros(M,1);
%%%%%%%%%%%HOW ABOUT Q WHEN NO BS OFF %%%%%%%%%%%%%
q_value=3.5;%Q value related to the distribution of cell size
Nt_factor=0;
Pt_dBm=23;
Pt_W=200;% in mW
Pt_BS=39811;
du=0.1;
%Power_Signal_Los=0;
Power_Signal_NLos=0;
Power_error=0;
Power_interference=0;
SINR_LoS=0;
min_r=1;
max_r=3000;
dr=1;
r00=min_r:dr:max_r;
length_r00=length(r00);
BS_Average_num=BS_density*4*max_r^2;
noise_rou2=0;
r_ll_1=min_r:dr:max_r;
r_00_1=min_r:dr:max_r;
length_r00_1=length(r_00_1);
length_r_ll=length(r_ll_1);
intermediate_P_Error=repmat({NaN(1,length_r_ll)},1,10);

  [BS_Density_DynOn,f_k_array,f_LoS]=getDynOnLamda(BS_density*1000000,UE_density*1000000,M,K_total);
save(f_k_name,'f_k_array');
% save('f_LoS_10_300','f_LoS');
% f_k_array=importdata('f_K_10_100_32_DL.mat');
% f_LoS=importdata('f_LoS_10_100.mat');
% f_LoS=[0.6,0.4];
%BS_Density_DynOn=BS_density;        % this code only work when UE density is much larger than BS density
%f_k_array=importdata('f_k_array_5_100.mat');        %The distribution of per-cell UE number get from simulation
%BS_Density_DynOn=BS_Density_DynOn/1000000;
BS_Density_DynOn=BS_density*(1-(1/(1+UE_density/(q_value*BS_density))^q_value));
k=1:1:K_total;
sum_k=0;
for kk=1:1:K_total
    sum_k=kk*f_k_array(kk)+sum_k;
end
k_average=sum_k-3;
% k_average=round(k_average)-1;
BS_Density_DynOn_ServeKth=(k_average/K_total)*BS_Density_DynOn;
% BS_density=BS_Density_DynOn;% test if the impact of BS density is big
% f_R_00_k_LoS=importdata('f_r00_L.mat');
% f_R_00_k_NLoS=importdata('f_r00_NL.mat');
% BS_Density_DynOn=BS_Density_DynOn_ServeKth;
%%
E_beta_ll_LoS=NaN(1,length_r_ll);
E_beta_ll_NLoS=NaN(1,length_r_ll);
%     Power_Signal_LoS=NaN(1,length_r00_1);
%     Power_Signal_NLoS=NaN(1,length_r00_1);
f_LoS_u_term_DynOn=@(u) 2*pi*BS_Density_DynOn.*u.*Pr_L(u);
f_NLoS_u_term_DynOn=@(u) 2*pi*BS_Density_DynOn.*u.*(1-Pr_L(u));
f_LoS_u_term_AllBS=@(u) 2*pi*BS_density.*u.*Pr_L(u);
f_NLoS_u_term_AllBS=@(u) 2*pi*BS_density.*u.*(1-Pr_L(u));
f_LoS_Scheduled_UE=@(u) A_LoS*2*pi*BS_Density_DynOn_ServeKth.*(u.^(1-alpha_LoS)).*Pr_L(u);
f_NLoS_Scheduled_UE=@(u) A_NLoS*2*pi*BS_Density_DynOn_ServeKth.*(u.^(1-alpha_NLoS)).*(1-Pr_L(u));
f_LoS_Scheduled_UE_IN2=@(u) ((A_LoS)^2)*2*pi*BS_Density_DynOn.*(u.^(1-2*alpha_LoS)).*Pr_L(u);
f_NLoS_Scheduled_UE_IN2=@(u) ((A_NLoS)^2)*2*pi*BS_Density_DynOn.*(u.^(1-2*alpha_NLoS)).*(1-Pr_L(u));
f_LoS_Scheduled_UE_IN2_P1=@(u) ((A_LoS)^1)*2*pi*BS_Density_DynOn.*(u.^(1-1*alpha_LoS)).*Pr_L(u);
f_NLoS_Scheduled_UE_IN2_P1=@(u) ((A_NLoS)^1)*2*pi*BS_Density_DynOn.*(u.^(1-1*alpha_NLoS)).*(1-Pr_L(u));
%%
R_e=sqrt(1/(BS_density*pi));%Exclusion ball model R_e also need modification (Larger R_e)
R_e_have_lthUE=sqrt(1/(BS_Density_DynOn_ServeKth*pi));
R_e_1=min(length_r00,R_e^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)));
R_e_2=min(length_r00,R_e^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)));
R_e_LoS=R_e;
R_e_NloS=R_e;
R_e_LoS_NLoS=318^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
R_e_NLoS_LoS=148^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
R_e_LoS_have_lthUE=R_e_have_lthUE;
R_e_NloS_have_lthUE=R_e_have_lthUE;
R_e_LoS_NLoS_have_lthUE=R_e_LoS_have_lthUE^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
R_e_NLoS_LoS_have_lthUE=R_e_LoS_NLoS_have_lthUE^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
%%

%     f_LoS=importdata('f_LoS.mat');
%     f_LoS(1)=f_LoS(1);
%      f_LoS(2)=f_LoS(2);

%     Power_Error=NaN(1,length_r00_1);
%     Power_Inter_Part2_intra=NaN(1,length_r00_1);
%     Power_Inter_Part2_inter=NaN(1,length_r00_1);
%     Power_Thermal_Noise=NaN(1,length_r00_1);
rou_noise_2=0;
N0 = 10^((-104 + 9)/10)*1e-3;
% K_total=10;
%%
%%%%%%%%%%Computing the power not related to r00_1%%%%%%%%%%%%%%%%%%%%%%%%%
%computing the interference power part1 intra-cell interference
r_00_k=min_r:dr:max_r;
length_r00k=length(r_00_k);
Beta_00_k_LoS=@(u) dr*A_LoS^(1-e)*2*pi*BS_density.*(u.^(alpha_LoS*e-alpha_LoS+1)).*Pr_L(u).*exp(-1*integral(f_LoS_u_term_AllBS,0,u)+...
    (-1*integral(f_NLoS_u_term_AllBS,0, u^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));%%BS density in this equation should be the dynamic On BS density
Beta_00_k_NLoS=@(u) dr*A_NLoS^(1-e)*2*pi*BS_density.*(u.^(alpha_NLoS*e-alpha_NLoS+1)).*(1-Pr_L(u)).*exp(-1*integral(f_NLoS_u_term_AllBS,0,u)+...
    (-1*integral(f_LoS_u_term_AllBS,0, u^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));%%Same as above one
f_LoS_r_00_1=@(u) 2*pi*BS_density*Pr_L(u)*u*exp(-1*integral(f_LoS_u_term_AllBS,0,u)+(-1*integral(f_NLoS_u_term_AllBS,0, u.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
f_NLoS_r_00_1=@(u) 2*pi*BS_density*(1-Pr_L(u))*u*exp(-1*integral(f_NLoS_u_term_AllBS,0,u)+(-1*integral(f_LoS_u_term_AllBS,0, u.^(alpha_NLoS/alpha_LoS)*...
    ((A_LoS/A_NLoS)^(1/alpha_LoS)))));

%      f_k_array=zeros(1,K_total);
k_array_nonzero=1:1:10;
k_array_zero=0:1:9;
%      for idx_k=1:10;
%          f_k_array(idx_k)=get_f_k(idx_k,BS_density_km,UE_density_km);
%      end

%      Power_Inter_P1_intra=sum(Power_Inter_P1_intra_single.*k_array_zero.*f_k_array);
%%%%%%%%%computing the interference power part1 inter-cell %%%%%%%%%%%%%%%
Beta_ll_k_LoS=@(u) dr*A_LoS^(-e)*2*pi*BS_density.*(u.^(alpha_LoS*e+1)).*Pr_L(u).*exp(-1*integral(f_LoS_u_term_AllBS,0,u)+...
    (-1*integral(f_NLoS_u_term_AllBS,0, u^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));%%BS density in this equation should be Dynamic ON BS density
Beta_ll_k_NLoS=@(u) dr*A_NLoS^(-e)*2*pi*BS_density.*(u.^(alpha_NLoS*e+1)).*(1-Pr_L(u)).*exp(-1*integral(f_NLoS_u_term_AllBS,0,u)+...
    (-1*integral(f_LoS_u_term_AllBS,0, u^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));%%BS Density in this equation should be Dynamic ON BS density
r_ll_k=min_r:dr:max_r;
length_rllk=length(r_ll_k);
% sum_betallk=0;
% for idx_r_ll_k=1:length_rllk
%     fprintf('The %d-th r_ll point of totally %d SINR points...111\n', idx_r_ll_k, length_rllk);
%     sum_betallk=sum_betallk+Beta_ll_k_LoS(r_ll_k(idx_r_ll_k))+Beta_ll_k_NLoS(r_ll_k(idx_r_ll_k));
% end
%%
E_beta_Campbell_ol_k_LoS=integral(f_LoS_Scheduled_UE_IN2,1*R_e_LoS,inf)+...
    integral(f_NLoS_Scheduled_UE_IN2,1*R_e_LoS_NLoS,inf); % as numerator components of IN2 when 1-th UE of l-th BS and l-th BS are LoS transmission
E_beta_Campbell_ol_k_NLoS=integral(f_LoS_Scheduled_UE_IN2,R_e_NLoS_LoS,inf)+...
    integral(f_NLoS_Scheduled_UE_IN2,1*R_e_NloS,inf); % as nuerator components of IN2 when 1-th UE of l-th BS and l-th BS are NLoS transmission
E_beta_Campbell_ol_k_LoS_all=integral(f_LoS_Scheduled_UE_IN2_P1,1*R_e_LoS,inf)+...
    integral(f_NLoS_Scheduled_UE_IN2,1*R_e_LoS_NLoS,inf); % as numerator components of IN2 when 1-th UE of l-th BS and l-th BS are LoS transmission
E_beta_Campbell_ol_k_NLoS_all=integral(f_LoS_Scheduled_UE_IN2_P1,R_e_NLoS_LoS,inf)+...
    integral(f_NLoS_Scheduled_UE_IN2,1*R_e_NloS,inf); % as nuerator components of IN2 when 1-th UE of l-th BS and l-th BS are NLoS transmission
% E_beta_Campbell_ol_k_LoS_all=integral(f_LoS_Scheduled_UE_IN2,1*R_e_LoS,inf)+...
%     integral(f_NLoS_Scheduled_UE_IN2,1*R_e_LoS_NLoS,inf); % as numerator components of IN2 when 1-th UE of l-th BS and l-th BS are LoS transmission
% E_beta_Campbell_ol_k_NLoS_all=integral(f_LoS_Scheduled_UE_IN2,R_e_NLoS_LoS,inf)+...
%     integral(f_NLoS_Scheduled_UE_IN2,1*R_e_NloS,inf); % as nuerator components of IN2 when 1-th UE of l-th BS and l-th BS are NLoS transmission

% k_LoS=k_average*(1);
% k_NLoS=k_average*f_LoS(2);
 delta_2_k=0;%sum_betallk*E_beta_Campbell_ol_k;
%%

E_beta_00k_1_e_NLoS=0;
E_beta_00k_e_LoS=0;
E_beta_00k_e_NLoS=0;
E_P5=@(k) 0;
E_P6=@(k) 0;
E_P7= 0;
E_P8= 0;
Power_Error_LoS_term=0;
E_beta_llk_e_sum_beta_0lk_LoS=zeros(1,length_r00);
E_beta_llk_e_sum_beta_0lk_NLoS=zeros(1,length_r00);
E_beta_llk_e_sum_beta_0lk_LoS_NLoS=zeros(1,length_r00);
E_beta_llk_e_sum_beta_0lk_NLoS_NLoS=zeros(1,length_r00);
% E_beta_llk_x_NLoS=importdata('E_beta_llk_x_NLoS.mat');
% E_beta_llk_x_LoS=importdata('E_beta_llk_x_LoS.mat');
E_beta_00k_1_e_LoS=0;
E_beta_00k_1_e_NLoS=0;
%%
flag_exist=0;

%%
if flag_exist==1;
    E_beta_llk_e_sum_beta_0lk_LoS=importdata('./E_beta/E_beta_64_10_100_32_08_LoS_LoS_1.mat');
    E_beta_llk_e_sum_beta_0lk_NLoS=importdata('./E_beta/E_beta_64_10_100_32_08_LoS_NLoS_2.mat');
    E_beta_llk_e_sum_beta_0lk_LoS_NLoS=importdata('./E_beta/E_beta_64_10_100_32_08_NLoS_LoS_3.mat');
    E_beta_llk_e_sum_beta_0lk_NLoS_NLoS=importdata('./E_beta/E_beta_64_10_100_32_08_NLoS_NLoS_4.mat');
else
    
    for idx_r=1:length_r00
        r_LoS_T_NLoS= r_00_k(idx_r).^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
        r_NLoS_T_LoS=r_00_k(idx_r).^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
        fprintf('The %d-th r_ll point of totally %d SINR points...111\n', idx_r, length_rllk);
        for idx_r_ll_k=1:length_rllk
            
            
            %%
            %many differences with UL case need to be modified cafefully
            E_beta_llk_e_sum_beta_0lk_LoS(idx_r)=E_beta_llk_e_sum_beta_0lk_LoS(idx_r)+dr*(A_LoS.*r_00_k(idx_r_ll_k).^(-alpha_LoS)).^(-e)*...
                (integral(f_LoS_Scheduled_UE,r_00_1(idx_r),inf)+...
                integral(f_NLoS_Scheduled_UE,r_LoS_T_NLoS,inf))*...%Delta_l for LoS case
                f_LoS_r_00_1(r_00_k(idx_r_ll_k));% As the components of computing the interference caused by other-cell scheduled UEs together with UL power control
            E_beta_llk_e_sum_beta_0lk_NLoS(idx_r)=E_beta_llk_e_sum_beta_0lk_NLoS(idx_r)+dr*(A_NLoS.*r_00_k(idx_r_ll_k).^(-alpha_NLoS)).^(-e)*...
                (integral(f_NLoS_Scheduled_UE,r_00_1(idx_r),inf)+...% Delta_l for NLoS case
                integral(f_LoS_Scheduled_UE,r_NLoS_T_LoS,inf))*...%same as above
                f_NLoS_r_00_1(r_00_k(idx_r_ll_k));
            E_beta_llk_e_sum_beta_0lk_LoS_NLoS(idx_r)=E_beta_llk_e_sum_beta_0lk_LoS_NLoS(idx_r)+dr*(A_LoS.*r_00_k(idx_r_ll_k).^(-alpha_LoS)).^(-e)*...
                (integral(f_LoS_Scheduled_UE,r_00_1(idx_r),inf)+...
                integral(f_NLoS_Scheduled_UE,r_LoS_T_NLoS,inf))*...%Delta_l for LoS case
                f_LoS_r_00_1(r_00_k(idx_r_ll_k));%Delta_l for LoS case when typical UE and typical BS are NLoS transmission
            E_beta_llk_e_sum_beta_0lk_NLoS_NLoS(idx_r)=E_beta_llk_e_sum_beta_0lk_NLoS_NLoS(idx_r)+dr*(A_NLoS.*r_00_k(idx_r_ll_k).^(-alpha_NLoS)).^(-e)*...
                (integral(f_NLoS_Scheduled_UE,r_00_1(idx_r),inf)+...% Delta_l for NLoS case
                integral(f_LoS_Scheduled_UE,r_NLoS_T_LoS,inf))*...%same as above
                f_NLoS_r_00_1(r_00_k(idx_r_ll_k));%Delta_l for LoS other cell interference calculated by campbell theorem when typical UE and typical BS are NLoS transmission
            if idx_r_ll_k==length_rllk
                save(E_betaLL_name,'E_beta_llk_e_sum_beta_0lk_LoS');
                save(E_betaLNL_name,'E_beta_llk_e_sum_beta_0lk_NLoS');
                save(E_betaNLL_name,'E_beta_llk_e_sum_beta_0lk_LoS_NLoS');
                save(E_betaNLNL_name,'E_beta_llk_e_sum_beta_0lk_NLoS_NLoS');
            end
            %%
            %         E_beta_00k_1_e_LoS=E_beta_00k_1_e_LoS+dr*(A_LoS.*r_00_k(idx_r_ll_k).^(-alpha_LoS)).^(1-e)*...
            %             f_LoS_r_00_1(r_00_k(idx_r_ll_k));% components for the 1-st part of averageble IN2 (denominotor) for LoS case
            %         E_beta_00k_1_e_NLoS=E_beta_00k_1_e_NLoS+dr*(A_NLoS.*r_00_k(idx_r_ll_k).^(-alpha_NLoS)).^(1-e)*...
            %             f_NLoS_r_00_1(r_00_k(idx_r_ll_k));%same as above for NLoS case
            %     E_beta_00k_e_LoS=E_beta_00k_e_LoS+dr*(A_LoS.*r_00_k(idx_r_ll_k).^(-alpha_LoS)).^(-e)*f_LoS_r_00_1(r_00_k(idx_r_ll_k));
            %     E_beta_00k_e_NLoS=E_beta_00k_e_NLoS+dr*(A_NLoS.*r_00_k(idx_r_ll_k).^(-alpha_NLoS)).^(-e)*f_NLoS_r_00_1(r_00_k(idx_r_ll_k));
            %     E_P5=@(k) E_P5(k)+ dr*Power_Inter5_00kL_ll1(r_00_k(idx_r_ll_k),k)*f_LoS_r_00_1(r_00_k(idx_r_ll_k));
            %     E_P6=@(k) E_P6(k)+dr*Power_Inter6_00kNL_ll1(r_00_k(idx_r_ll_k),k)*f_NLoS_r_00_1(r_00_k(idx_r_ll_k));
            %     E_P7= E_P7+dr*Power_Inter7_llkL_ll1(r_00_k(idx_r_ll_k))*f_LoS_r_00_1(r_00_k(idx_r_ll_k));
            %     E_P8= E_P8+dr*Power_Inter8_llkNL_ll1(r_00_k(idx_r_ll_k))*f_NLoS_r_00_1(r_00_k(idx_r_ll_k));
        end
    end
end
 
 for idx_r_ll_k=1:length_rllk
     E_beta_00k_1_e_LoS=E_beta_00k_1_e_LoS+dr*(A_LoS.*r_00_k(idx_r_ll_k).^(-alpha_LoS)).^(1-e)*...
         f_LoS_r_00_1(r_00_k(idx_r_ll_k));% components for the 1-st part of averageble IN2 (denominotor) for LoS case
     E_beta_00k_1_e_NLoS=E_beta_00k_1_e_NLoS+dr*(A_NLoS.*r_00_k(idx_r_ll_k).^(-alpha_NLoS)).^(1-e)*...
         f_NLoS_r_00_1(r_00_k(idx_r_ll_k));%same as above for NLoS case
 end
    
    Power_Signal_LoS=@(x) ...   % consider as a function of the distance between typical UE and original BS(r00) when ty[ical BS and UE are LoS transmission
        (M+1)*Pt_BS./k_average.*(A_LoS.*x.^(-1*alpha_LoS)).^(2-e);
    Power_Signal_NLoS=@(x) ...   % a function of r_001 when the typical BS and UE are NLoS transmission
        (M+1)*Pt_BS./k_average.*(A_NLoS.*x.^(-1*alpha_NLoS)).^(2-e);% consider as a function of the distance between typical UE and original BS(r00)
    Power_Inter_P1_LoS= @(x) (((A_LoS.*x.^(-1*alpha_LoS)).^(1-e)+rou_noise_2+ E_beta_llk_e_sum_beta_0lk_LoS(x)+E_beta_llk_e_sum_beta_0lk_NLoS(x))).*...
        (Pt_BS*E_beta_Campbell_ol_k_LoS_all+Pt_BS.*(A_LoS.*x.^(-1*alpha_LoS)));
%     +(Pt_BS/k_average*M*E_beta_Campbell_ol_k_LoS/(...
%         f_LoS(1)/f_LoS(1)*(E_beta_00k_1_e_LoS+E_beta_llk_e_sum_beta_0lk_LoS+E_beta_llk_e_sum_beta_0lk_NLoS)+...
%         f_LoS(2)/f_LoS(2)*(E_beta_00k_1_e_NLoS+E_beta_llk_e_sum_beta_0lk_LoS_NLoS+E_beta_llk_e_sum_beta_0lk_NLoS_NLoS))));
    Power_Inter_P1_NLoS= @(x) (((A_NLoS.*x.^(-1*alpha_NLoS)).^(1-e)+rou_noise_2+ E_beta_llk_e_sum_beta_0lk_LoS_NLoS(x)+E_beta_llk_e_sum_beta_0lk_NLoS_NLoS(x))).*...
        (Pt_BS*E_beta_Campbell_ol_k_NLoS_all+Pt_BS.*(A_NLoS.*x.^(-1*alpha_NLoS)));
%         +(Pt_BS/k_average*M*E_beta_Campbell_ol_k_NLoS/(...
%         f_LoS(1)/f_LoS(1)*(E_beta_00k_1_e_LoS+E_beta_llk_e_sum_beta_0lk_LoS+E_beta_llk_e_sum_beta_0lk_NLoS)+...
%         f_LoS(2)/f_LoS(2)*(E_beta_00k_1_e_NLoS+E_beta_llk_e_sum_beta_0lk_LoS_NLoS+E_beta_llk_e_sum_beta_0lk_NLoS_NLoS)))...
%         );
    Power_Inter_P2_LoS= @(x) Pt_BS/k_average.*(A_LoS.*x.^(-1*alpha_LoS)).^(2-e);
    Power_Inter_P2_NLoS= @(x) Pt_BS/k_average.*(A_NLoS.*x.^(-1*alpha_NLoS)).^(2-e);
    SINR_LoS_function=@(r00) Power_Signal_LoS(r00)./(Power_Inter_P1_LoS(r00)-Power_Inter_P2_LoS(r00));
    SINR_NLoS_function=@(r00) Power_Signal_NLoS(r00)./(Power_Inter_P1_NLoS(r00)-Power_Inter_P2_NLoS(r00));
    
    
    %%
  N_gamma=13;
sigma=(N_gamma)*(factorial(N_gamma)^(-1/(N_gamma)));
intermediate_SINR= @(x,T) 0;
CCDE=@(x,T) 0;

for exp_N=1:N_gamma
    %            f_LoS_r_00_1=@(u) 2*pi*BS_density*Pr_L(u)*u*exp(-1*integral(f_LoS_u_term_AllBS,0,u)+(-1*integral(f_NLoS_u_term_AllBS,0, u.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
    %            f_NLoS_r_00_1=@(u) 2*pi*BS_density*(1-Pr_L(u))*u*exp(-1*integral(f_NLoS_u_term_AllBS,0,u)+(-1*integral(f_LoS_u_term_AllBS,0, u.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    scale=prod(exp_N+1:N_gamma)/(prod(1:(N_gamma-exp_N)))*(-1)^(exp_N+1);
    if exp_N==1
        intermediate_SINR_1=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_1_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==2
        intermediate_SINR_2=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_2_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==3
        intermediate_SINR_3=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_3_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==4
        intermediate_SINR_4=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_4_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==5
        intermediate_SINR_5=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_5_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==6
        intermediate_SINR_6=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_6_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==7
        intermediate_SINR_7=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_7_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==8
        intermediate_SINR_8=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_8_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==9
        intermediate_SINR_9=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_9_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==10
        intermediate_SINR_10=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_10_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==11
        intermediate_SINR_11=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_11_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==12
        intermediate_SINR_12=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_12_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    if exp_N==13
        intermediate_SINR_13=@(x,T) 2*pi*BS_density*Pr_L(x)*x*scale.*exp((-1*exp_N*T*sigma./SINR_LoS_function(x))+...
            (-1*integral(f_LoS_u_term_AllBS,0,x))+(-1*integral(f_NLoS_u_term_AllBS,0, x.^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS)))));
        intermediate_SINR_13_NLoS=@(x,T) 2*pi*BS_density*(1-Pr_L(x))*x*scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x)+...
            (-1*integral(f_NLoS_u_term_AllBS,0,x))+(-1*integral(f_LoS_u_term_AllBS,0, x.^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS)))));
    end
    %            if exp_N==14
    %                intermediate_SINR_14=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_14_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    %            if exp_N==15
    %                intermediate_SINR_15=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_15_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    %             if exp_N==16
    %                intermediate_SINR_16=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_16_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    %            if exp_N==17
    %                intermediate_SINR_17=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_17_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    %            if exp_N==18
    %                intermediate_SINR_18=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_18_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    %            if exp_N==19
    %                intermediate_SINR_19=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_19_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    %            if exp_N==20
    %                intermediate_SINR_20=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_LoS_function(x));
    %                intermediate_SINR_20_NLoS=@(x,T) scale.*exp(-1*exp_N*T*sigma./SINR_NLoS_function(x));
    %            end
    
    % scale.*intermediate_SINR(x,T).^(exp_N)
    
end
intermediate_SINR=@(x,T) intermediate_SINR_1(x,T)+intermediate_SINR_2(x,T)+intermediate_SINR_3(x,T)+intermediate_SINR_4(x,T)+intermediate_SINR_5(x,T)+...
    intermediate_SINR_6(x,T)+intermediate_SINR_7(x,T)+intermediate_SINR_8(x,T)+intermediate_SINR_9(x,T)+intermediate_SINR_10(x,T)+...
    intermediate_SINR_11(x,T)+intermediate_SINR_12(x,T)+intermediate_SINR_13(x,T);%+...
%            intermediate_SINR_16(x,T)+intermediate_SINR_17(x,T)+intermediate_SINR_18(x,T)+intermediate_SINR_19(x,T)+intermediate_SINR_20(x,T);
intermediate_SINR_NLoS=@(x,T) intermediate_SINR_1_NLoS(x,T)+intermediate_SINR_2_NLoS(x,T)+intermediate_SINR_3_NLoS(x,T)+intermediate_SINR_4_NLoS(x,T)+intermediate_SINR_5_NLoS(x,T)+...
    intermediate_SINR_6_NLoS(x,T)+intermediate_SINR_7_NLoS(x,T)+intermediate_SINR_8_NLoS(x,T)+intermediate_SINR_9_NLoS(x,T)+intermediate_SINR_10_NLoS(x,T)+...
    intermediate_SINR_11_NLoS(x,T)+intermediate_SINR_12_NLoS(x,T)+intermediate_SINR_13_NLoS(x,T);%+...
%            intermediate_SINR_16_NLoS(x,T)+intermediate_SINR_17_NLoS(x,T)+intermediate_SINR_18_NLoS(x,T)+intermediate_SINR_19_NLoS(x,T)+intermediate_SINR_20_NLoS(x,T);

intermediate_SINR_Ek=@(x,T) 0;
intermediate_SINR_NLoS_Ek=@(x,T) 0;
Threshold_dB=-10:2:30;
Threshold=10.^(Threshold_dB./10);
length_Threshold=length(Threshold);
sum_Pr=zeros(1,length_Threshold);
sum_Pr_NLoS=zeros(1,length_Threshold);
sum_Pr_final=zeros(1,length_Threshold);
%        for k=1:K_total
%            intermediate_SINR_Ek=@(x,T) intermediate_SINR(x,T,k)*f_k_array(k)+intermediate_SINR_Ek(x,T);
%            intermediate_SINR_NLoS_Ek=@(x,T) intermediate_SINR_NLoS(x,T,k)*f_k_array(k)+intermediate_SINR_NLoS_Ek(x,T);
%        end
r_final=0.5:0.5:500;
for idx_T=1:length_Threshold
    %temp=@(x) intermediate_SINR(x,Threshold_SINR(idx_T)).*f_LoS_r_00_1(x);
    %temp_in=integral(temp,0,1);
    for idx_r_ll_k=1:length_r00_1
        fprintf('The %d-th r_ll point of totally %d SINR points...333\n', idx_r_ll_k, length_rllk);
        sum_Pr(idx_T)=sum_Pr(idx_T)+dr*intermediate_SINR(r_00_1(idx_r_ll_k),Threshold(idx_T));
        sum_Pr_NLoS(idx_T)=sum_Pr_NLoS(idx_T)+dr*intermediate_SINR_NLoS(r_00_1(idx_r_ll_k),Threshold(idx_T));
    end
    if Threshold(idx_T)==1
        cov_pro=sum_Pr(idx_T)+sum_Pr_NLoS(idx_T);
        save(ccdf_name,'cov_pro');
    end
    
end

sum_Pr_final=sum_Pr+sum_Pr_NLoS;
save(ccdf_threshold,'sum_Pr_final');

hold on;
xlabel('SINR Threshold in dB');
ylabel('CCDF_analytical');
figure(666)
hold on ;
plot(Threshold_dB,sum_Pr_final,'o');
grid on;
box on;