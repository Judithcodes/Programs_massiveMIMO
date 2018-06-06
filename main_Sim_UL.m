function []=main_Sim_UL(total_antennas,BS_density,UE_density_array,kk,ee)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: main_simmulation uplink massive mimo limited users
% Authors: Xuefeng Yao
% Version: 3.0
% Date: 15/3/2017
% Description: Main function for massive mimo uplink HPPP(PPP) model considering
% limited users ans and Los and Nlos 
% 
% 
% 
% 
% Copyright(c): For personal study only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%This version works but not perfect %%%%%%%%%%%%%%%%%%%%%%%%%%%
% In tihs version the UEs are generated total randomly. TRY TO GENERATE UEs
% IN EACH BS USING NEGATIVE BINAMIAL DISTRIBUTION MIGHT BE BETTER!!!!
total_antennas_str=num2str(total_antennas);
BS_density_str=num2str(BS_density);
UE_density_str_first=num2str(UE_density_array(1));
UE_density_str_last=num2str(UE_density_array(end));
ee_str=num2str(ee*10);
kk_str=num2str(kk);
ase_name=['CCDFSIM',total_antennas_str,'_',BS_density_str,'_',UE_density_str_first,'_',UE_density_str_last,'_',ee_str,'_',kk_str];
ase_name_x=['CCDFSIM',total_antennas_str,'_',BS_density_str,'_',UE_density_str_first,'_',UE_density_str_last,'_',ee_str,'_',kk_str,'_x'];




%%
% set some global variables 
global d0;
global Alpha_LoS;
global Alpha_NLoS;
global A_LoS;
global A_NLoS;
global N0;
global SINR_thresold;
global M
global global_Test;
% variables settings 
M=64;
%Pt_dBm=23;
Pt_W=200;
Alpha_LoS=2.09;
Alpha_LoS_Pico=2.09;
Alpha_LoS_Macro=2.09;
Alpha_NLoS=3.75;
Alpha_NLoS_Pico=3.75;
Alpha_NLoS_Macro=3.75;
A_LoS_dB = 41.1;
A_LoS_dB_Pico=41.1;
A_LoS_dB_Macro=41.1;%30.8;
A_NLoS_dB = 32.9;
A_NLoS_dB_Pico=32.9;
A_NLoS_dB_Macro=32.9;%;
A_LoS = 10^(-1*A_LoS_dB/10);
A_NLoS = 10^(-1*A_NLoS_dB/10);
A_LoS_Macro = 10^(-1*A_LoS_dB_Macro/10);
A_NLoS_Macro = 10^(-1*A_NLoS_dB_Macro/10);
A_LoS_Pico = 10^(-1*A_LoS_dB_Pico/10);
A_NLoS_Pico= 10^(-1*A_NLoS_dB_Pico/10);
%BS_density=4;
%len_BS_density=length(BS_density);
max_r=2e3;
enlarged_factor=1.5;
e=ee;

exp_num=1000;
 Nt_noise=zeros(M,1);
 
  UE_density_test=UE_density_array;
 length_Diff_UE_density=length(UE_density_test);
%  BS_density_test=[10,100];
%  BS_num_antennas=1000./BS_density_test;
%  length_BS_test=length(BS_density_test);

%  BS_HEight=30;
 UE_noise_figure = 9;
 PICO_noise_figure = 13; 
 N0 = 10^((-104 + UE_noise_figure)/10)*1e-3;     % 10MHz!!!
%  N0_UL = 10^((-104 + PICO_noise_figure)/10)*1e-3;     % 10MHz!!!
SINR_Threshold_array_dB=-80:0.5:40;
Threshold=10.^(SINR_Threshold_array_dB/10);
SINR_CCDF_dB=NaN(1,length(SINR_Threshold_array_dB)-1);
SINR_CCDF=NaN(length_Diff_UE_density,length(SINR_Threshold_array_dB)-1);

area_spec_efficiency=zeros(1,length_Diff_UE_density);
x_array = 1/2 * (Threshold(1:end-1) + Threshold(2:end));
usable_x_idx = find(x_array >= 1, 1);
dx_array=Threshold(2:end)-Threshold(1:end-1);
length_FixedThreshold=length(x_array);
x_array_dB=10*log10(x_array);
BS_density_Dynamic=zeros(1,exp_num);
ASE_array_exp=zeros(1,exp_num);
BS_density_DynOnOff=zeros(1,length_Diff_UE_density);
%%
      if BS_density<11
            Alpha_LoS=Alpha_LoS_Macro;
            Alpha_NLoS=Alpha_NLoS_Macro;
            A_LoS=A_LoS_Macro;
            A_NLoS=A_NLoS_Macro;
            Pt_W=200;
            Pt_W_UE=200;
            M=total_antennas;
            BS_HEight=0;
            K_total=kk;
            PL_Model='36.828Macro';
      else
            Alpha_LoS=Alpha_LoS_Pico;
            Alpha_NLoS=Alpha_NLoS_Pico;
            A_LoS=A_LoS_Pico;
            A_NLoS=A_NLoS_Pico;
            Pt_W=0.2;
            Pt_W_UE=0.2;
            M=total_antennas;
            BS_HEight=0;
            K_total=kk;
            PL_Model='36.828Macro';
      end

     %repmat({NaN(len_BS_density_array, exp_num)}, 1, otherPowers_num);
%      SINR_dB=zeros(length_Diff_UE_density,exp_num);
     SINR_dB=repmat({zeros(1,K_total)},length_Diff_UE_density,exp_num);
     SINR=repmat({zeros(1,K_total)},length_Diff_UE_density,exp_num);
     SINR_NLoS=zeros(1,exp_num);
     SINR_LoS=zeros(1,exp_num);
     %SINR=zeros(length_Diff_UE_density,exp_num);
     SINR_Threshold_array_dB=-20:0.01:40;
     Threshold=10.^(SINR_Threshold_array_dB./10);
     SINR_CCDF_dB=NaN(1,length(SINR_Threshold_array_dB));
     SINR_CCDF=NaN(length_Diff_UE_density,length(SINR_Threshold_array_dB));
     %%
     test_ave_SINR=0;
     test_n=0;
     flag_NLoS=0;
     flag_LoS=0;
for idx_UE_density_test=1:length_Diff_UE_density

    for exp_idx=1:exp_num
        flag_out=0;
    % Report the simulation progress
            if mod(exp_idx, exp_num/100) == 0
                fprintf('\nSimulation progress: %3.1f%% completed......', (exp_idx/exp_num*100));
            end
            enlarged_considered_area=(enlarged_factor*max_r/1000*2)^2;
            %HPPP for BS
            total_BS_num=0;
            while total_BS_num==0
                total_BS_num=poissrnd(BS_density*enlarged_considered_area);
            end
            flag_Ful_BS=zeros(1,total_BS_num);
            %%%%%%%%%%set BS related matrix used fro caculation%%%%%%%%%%%%% 
            UE_servedBy_BS=zeros(total_BS_num,K_total);
            UE_typical_BS=zeros(total_BS_num,K_total);
            flag_UE_servedBy_BS=zeros(total_BS_num,K_total);
            distance_UE_BS=zeros(1,total_BS_num);
            distance_UE_typical=zeros(1,total_BS_num);
            PathLoss_UE_BS=zeros(1,total_BS_num);
            check_BS_full=zeros(1,K_total);
            BS_Avaliable_sequence=zeros(total_BS_num,K_total);
            
             %%%%%%%%%%%%set BS related matrix for caculation%%%%%%%%%%%%%%%
            BS_coordinates = [0, unifrnd(-1*enlarged_factor*max_r,enlarged_factor*max_r, 1, total_BS_num-1) ...
                + 1i * unifrnd(-1*enlarged_factor*max_r,enlarged_factor*max_r, 1, total_BS_num-1)]; %BS coordinates with typical BS at original
            total_UE_num=0;
            while total_UE_num==0
                total_UE_num=poissrnd(UE_density_test(idx_UE_density_test)*enlarged_considered_area);
            end
            UE_coordinates=unifrnd(-1*enlarged_factor*max_r,enlarged_factor*max_r,1,total_UE_num)+1i*unifrnd(-1*enlarged_factor*max_r,enlarged_factor*max_r,1,total_UE_num);
            %%%%%%%%%%%%set UE related matrix for caculation%%%%%%%%%%%%%%
            Power_compensation_UE=zeros(1,total_UE_num);
            
            Served_BS_idx=zeros(1,total_UE_num);
            PathLoss_UE_typical=zeros(1,total_UE_num);
            flag_BS_UE=zeros(total_BS_num,K_total);
            sequence_k=1:1:K_total;
            flag_BS_avaliable_sequence=cell(1,total_BS_num);
            
          
            %%%%%%%%%%%%%set UE related matrix for caculation%%%%%%%%%%%%% 
           %%%%%%%%%%%%%%|||||||||||||||||||||||||||||||||||||||||||||||%%%%%%%%%%%
           %%%%%%%%%%%%%set the avaliable status matrix of the polit sequences in BS %%%%%%%%%%%%%%%%
            for BS_temp_idx=1:total_BS_num
                BS_Avaliable_sequence(BS_temp_idx,:)=1:1:K_total;    
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %traverse all UE in enlarged area
            for UE_idx=1:total_UE_num
%                      if mod(UE_idx, total_UE_num/100) == 0
%                     fprintf('\nUE traverse progress: %3.1f%% completed......', (UE_idx/total_UE_num*100));
%                      end
                    for BS_idx=1:total_BS_num
                        distance_UE_BS(BS_idx)=abs(UE_coordinates(UE_idx)-BS_coordinates(BS_idx));
                        distance_UE_BS(BS_idx)=sqrt(distance_UE_BS(BS_idx)^2+BS_HEight^2);
                        
                        if rand(1)<(get_Pr_LoS(distance_UE_BS(BS_idx),PL_Model))
                            PathLoss_UE_BS(BS_idx)=A_LoS*distance_UE_BS(BS_idx)^(-1*Alpha_LoS);
                        else
                            PathLoss_UE_BS(BS_idx)=A_NLoS*distance_UE_BS(BS_idx)^(-1*Alpha_NLoS);
                        end

                    end
                    if UE_idx==1&&(abs(real(BS_coordinates(1)))>max_r||abs(imag(BS_coordinates(1)>max_r)))
                        flag_out=1;
                        
                    end
                     [max_pathloss_UE, BS_maxpathloss_idx]=max(PathLoss_UE_BS);
                    PL_currentUE_TypicalBS=PathLoss_UE_BS(1);
                    if  flag_Ful_BS(BS_maxpathloss_idx)==1;
                        %BS_maxpathloss_idx
                        continue;
                    end
                    
                    
                    distance_UE_typical(UE_idx)=abs(UE_coordinates(UE_idx));
                    temp_avaliable_sequence=BS_Avaliable_sequence(BS_maxpathloss_idx,:);
                    temp_avaliable_idx=find(temp_avaliable_sequence>0);
                    if isempty(temp_avaliable_idx)
                        flag_Ful_BS(1,BS_maxpathloss_idx)=1;
                        %BS_maxpathloss_idx
                        continue;
                    end
                    temp_k=ceil(rand(1)*length(temp_avaliable_idx));
                    temp_k_real=temp_avaliable_idx(temp_k);
                    if (BS_maxpathloss_idx)==1&&(temp_avaliable_idx(1)==1)
                        temp_k_real=1;
                        r_001=abs(UE_coordinates(UE_idx)-BS_coordinates(BS_maxpathloss_idx));
                        if max_pathloss_UE==A_NLoS*r_001^(-Alpha_NLoS)
                            flag_NLoS=1;
                        else
                            flag_LoS=1;

                        end
                    end
                    BS_Avaliable_sequence(BS_maxpathloss_idx,temp_k_real)=0;
                    UE_servedBy_BS(BS_maxpathloss_idx,temp_k_real)=max_pathloss_UE;
                    UE_typical_BS(BS_maxpathloss_idx,temp_k_real)=PL_currentUE_TypicalBS;                                   
            end
            %%
            if  flag_out==1
                continue;
            end
            flag=0;
            for idx_BS_l=1:total_BS_num
                if abs(real(BS_coordinates(idx_BS_l)))>max_r||abs(imag(BS_coordinates(idx_BS_l)))>max_r
                    UE_servedBy_BS(idx_BS_l-flag,:)=[];
                    UE_typical_BS(idx_BS_l-flag,:)=[];
                    flag=flag+1;
                end
            end
            if UE_servedBy_BS(1,1)==0
                continue;
            end
            BS_density_Dynamic(exp_idx)=length(find(sum(UE_servedBy_BS,2)~=0))/(((max_r/1000)*2)^2);
            UE_arrayIntypical=UE_servedBy_BS(1,:);
            Idx_arrayInTypical=find(UE_arrayIntypical);
            num_UEInTypical=length(UE_arrayIntypical(UE_arrayIntypical>0));
            %%
            sum_rate=0;
            for idx_InTypicalBS=1:num_UEInTypical
                idx_current=Idx_arrayInTypical(idx_InTypicalBS);
                total_BS_num_thining=length(UE_servedBy_BS(:,idx_current));
                U00_Matrix=zeros(M,total_BS_num_thining);
                %%%%%%%%%%Uplink Trainig Stage%%%%%%%%%%%%%%%%%
                %UE_In_Typical_BS=UE_servedBy_BS(1,:);
                Idx_Tagged_UE=idx_current;
                %Idx_Tagged_UE
                UE_Tagged=UE_servedBy_BS(1,Idx_Tagged_UE);
                if isempty(UE_Tagged)
                    continue;
                end
                
                %Note that the zero elememts at UE_BS matrix
                %should be removed before caculating the power
                %compensation element
                %                                             if Idx_Tagged_UE~=1
                %                                                     temp_a=UE_servedBy_BS(1,Idx_Tagged_UE);
                %                                                     UE_servedBy_BS(1,Idx_Tagged_UE)=UE_servedBy_BS(1,1);
                %                                                     UE_servedBy_BS(1,1)=temp_a;
                %                                                     temp_b=UE_typical_BS(1,Idx_Tagged_UE);
                %                                                     UE_typical_BS(1,Idx_Tagged_UE)=UE_typical_BS(1,1);
                %                                                     UE_typical_BS(1,1)=temp_b;
                %                                                     Idx_Tagged_UE=1;
                %                                             end
                if isempty(UE_servedBy_BS)
                    continue;
                end
                
                UE_Same_Sequence_Served=UE_servedBy_BS(:,Idx_Tagged_UE);
                idx_BSHasUESameAsTypicalUE=find(UE_Same_Sequence_Served>0);
                UE_Same_Sequence_Served=UE_Same_Sequence_Served(UE_Same_Sequence_Served>0);
                UE_Same_Sequence_Typical=UE_typical_BS(:,Idx_Tagged_UE);
                PL_UE_Same_sequence_with_typicalUE=UE_Same_Sequence_Typical(UE_Same_Sequence_Typical>0);
                %%%%%%%%%%%fix till here                       PL_UE_Same_sequence_with_typicalUE=PL_UE_Same_sequence_with_typicalUE(PL_UE_Same_sequence_with_typicalUE>0);
                total_Nonzero=length(UE_Same_Sequence_Served);
                U00_observation=zeros(M,1);
                % Channel_rayleigh=randn(M,1)+1i*randn(M,1);
                Channel_training=zeros(M,total_Nonzero);
                
                
                %%%%%%%%Caculating the Observation of U00 %%%%%%%%%%%%%%%%%%
                Pathloss_factor=sqrt(Pt_W_UE*(UE_Same_Sequence_Served.^-e).*(PL_UE_Same_sequence_with_typicalUE));
                
                for idx_same=1:total_Nonzero
                    Channel_training(:,idx_same)=sqrt(1/2)*(randn(M,1)+1i*randn(M,1));
                    U00_Matrix(:,idx_same)=Pathloss_factor(idx_same,1)*sqrt(eye(M))*(Channel_training(:,idx_same));
                    if idx_same==1
                        Channel_Original=U00_Matrix(:,1)./sqrt(Pt_W_UE*UE_Same_Sequence_Served(1,1).^-e);
                    end
                end
                %
                %
                %
                U00_observation=sum(U00_Matrix,2);
                %%%%%%%%Caculating the Observation of U00 %%%%%%%%%%%%%%%%%%
                
                Nt_noise=0;%sqrt(1/2)*(randn(M,1)+1i*randn(M,1));
                U00_observation=U00_observation+Nt_noise;
                g_00=U00_observation';
                %%%%Estimate the channel Using MMSE Estimator
                Sum_Pathloss_factor=Pt_W_UE*((UE_Same_Sequence_Served.^-e).*PL_UE_Same_sequence_with_typicalUE);
                %UE_Same_Sequence_Served.^-e
                %Sum_Pathloss_factor=Pt_W.*UE_servedBy_BS(:,Idx_Tagged_UE).^(-e).*UE_typical_BS(:,Idx_Tagged_UE)+Nt_factor;
                Sum_Pathloss_factor=sum(Sum_Pathloss_factor);
                Channel_Estimated=((sqrt(Pt_W_UE)*sqrt(UE_Same_Sequence_Served(1,1).^-e)*PL_UE_Same_sequence_with_typicalUE(1,1))/Sum_Pathloss_factor)*U00_observation;
                Channel_Error=Channel_Original-Channel_Estimated;%1*M vector


                %Channel_Estimated=(sqrt(1/2)*(randn(M,1)+1i*randn(M,1)));
                %Channel_Estimated=U00_observation;
                %g_00=Channel_Original;
                
                
                
                
                %%%%%%%%%%Uplink Training Stage%%%%%%%%%%%%%%%%%
                %%
                Power_Interference=0;
                %%%%%%%%%%%%%%Caculating The SINR%%%%%%%%%%%%%%%
                temp_power_sig=abs(g_00*Channel_Estimated);
                temp_power_sig=temp_power_sig^2;
                Power_Signal=Pt_W*(UE_Same_Sequence_Served(1,1).^-e)*temp_power_sig;
                UEs_InTypicalBS=UE_servedBy_BS(1,:);
                
                %Power_Signal=Pt_W*UE_servedBy_BS(1,Idx_Tagged_UE).^(-e)*sum(abs(Channel_Estimated.*Channel_Estimated).^2);
                temp_power_error=(abs(g_00*Channel_Error))^2;
                Power_Estimate_Error=Pt_W*(UE_Same_Sequence_Served(1,1).^-e)*temp_power_error;
                

                %PL_Inter_BS=zeros(1,total_BS_num_thining);
                
                %                                 for idx_Inter_BS=1:total_BS_num_thining
                %
                % %                                      if Distance_BS_TypicalUE(idx_Inter_BS)<300
                % %                                                 if rand(1)<(1-Distance_BS_TypicalUE(idx_Inter_BS)/300)
                % %                                                    PL_Inter_BS(idx_Inter_BS)=A_LoS*Distance_BS_TypicalUE(idx_Inter_BS)^(-Alpha_LoS);
                % %                                                 else
                % %                                                     PL_Inter_BS(idx_Inter_BS)=A_NLoS*Distance_BS_TypicalUE(idx_Inter_BS)^(-Alpha_NLoS);
                % %                                                 end
                % %                                             else
                % %                                                 PL_Inter_BS(idx_Inter_BS)=A_NLoS*Distance_BS_TypicalUE(idx_Inter_BS)^(-Alpha_NLoS);
                % %                                       end
                %                                 end
                
                for idx_l=1:total_BS_num_thining
                    
                    temp_Served=UE_servedBy_BS(idx_l,:);
                    idx_InCurrentBS=find(temp_Served>0);
                    temp_Served_nonZero=temp_Served(temp_Served>0);
                    idx_InCurrentBS=find(temp_Served);
                    temp_Typical=UE_typical_BS(idx_l,:);
                    
                    temp_Typical_nonZero=temp_Typical(temp_Typical>0);
                    K_InCurrentBS=length(temp_Served_nonZero);
                    
                    for idx_k=1:K_InCurrentBS
                        idx_k_real=idx_InCurrentBS(idx_k);
                        if idx_k_real==idx_InTypicalBS&&idx_l==1
                            continue;
                        end
                        
                        if idx_k_real==10000000000
                            temp_Typical_nonZero
                            idx_k
                            idx_k_real
                            Channel_training
                            %temp_power_inter=(g_00.*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*ones(M,1)))'*(g_00.*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*ones(M,1)));
                            temp_power_inter=(abs(g_00*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*Channel_training(:,idx_k_real))))^2;
                            %temp_power_inter=(g_00.*(sqrt(temp_Typical_nonZero(idx_k))))'*(g_00.*(sqrt(temp_Typical_nonZero(idx_k))));
                        else
                            channel_Inter=sqrt(1/2)*(randn(M,1)+1i*randn(M,1));
                            %temp_power_inter=(g_00.*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*ones(M,1)))'*(g_00.*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*ones(M,1)));
                            temp_power_inter=(abs(g_00*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*channel_Inter))^2);
                            %temp_power_inter=(g_00.*(sqrt(temp_Typical_nonZero(idx_k))))'*(g_00.*(sqrt(temp_Typical_nonZero(idx_k))));
                            
                        end
                        %temp_power_inter=(g_00.*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*ones(M,1)))'*(g_00.*(sqrt(temp_Typical_nonZero(idx_k))*sqrt(eye(M))*ones(M,1)));
                        if idx_l==1
                            Channel_BS_TypicalUE=Channel_Original;
                        end
                        
                        
                        
                        %((sqrt(Pt_W)*sqrt(UE_Same_Sequence_Served(1,1).^-e)*UE_Same_Sequence_Typical(1,1))/Sum_Pathloss_factor)
                        Power_Interference=Power_Interference+Pt_W*(temp_Served_nonZero(idx_k).^-e)*temp_power_inter;
                        % temp_Served_nonZero(idx_k).^-e
                    end
                end
                
                Power_Noise=N0;
                SINR_dB{idx_UE_density_test,exp_idx}(idx_InTypicalBS)=10*log10(Power_Signal/(Power_Estimate_Error+Power_Interference+Power_Noise));
                SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS)=Power_Signal/(Power_Estimate_Error+Power_Interference+Power_Noise);
                if flag_NLoS==1
                   SINR_NLoS(exp_idx)=Power_Signal/(Power_Estimate_Error+Power_Interference+Power_Noise); 
                   flag_NLoS=0; 
                end
                if flag_LoS==1
                   SINR_LoS(exp_idx)=Power_Signal/(Power_Estimate_Error+Power_Interference+Power_Noise); 
                   flag_LoS=0; 
                end
                SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS);
                test_ave_SINR=test_ave_SINR+SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS);
                test_n=test_n+1;
                test_ave_SINR/test_n;
                %                                             SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS)
%                 if SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS)<1
%                     SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS)=0;
%                 end
                sum_rate=log2(1+ SINR{idx_UE_density_test,exp_idx}(idx_InTypicalBS))+sum_rate;
                sum_rate;
                %10*log10(SINR)
                %%%%%%%%%%%%%%Caculating The SINR%%%%%%%%%%%%%%%
                %SINR(idx_UE_density_test,exp_idx)
            end
                           
                    
                    
    end    
    for idx=1:length(SINR_Threshold_array_dB)
        num_consideredUE_NLoS=0;
        num_qualifiedUE_NLoS=0;
        for idx_exp=1:exp_num
            num_consideredUE_NLoS=num_consideredUE_NLoS+length(SINR{idx_UE_density_test,idx_exp}(SINR{idx_UE_density_test,idx_exp}>0));
            num_qualifiedUE_NLoS=num_qualifiedUE_NLoS+length(SINR{idx_UE_density_test,idx_exp}(SINR{idx_UE_density_test,idx_exp}>Threshold(idx)));
        end
        
        SINR_CCDF(idx_UE_density_test,idx)=num_qualifiedUE_NLoS/num_consideredUE_NLoS;
%         save('Threshold_x_e08_U100_B10_K10.mat','SINR_Threshold_array_dB');
%         save('CCDF_e08_U100_B10_K10.mat','SINR_CCDF');
%         SINR_CCDF_dB(idx_test,idx)=length(find(SINR_dB(idx_test,:)>SINR_Threshold_array_dB(idx)))/exp_num;
    end
    hold on;
    plot(SINR_Threshold_array_dB,SINR_CCDF);
    grid on;
    box on;
    save(ase_name,'SINR_CCDF');
end

% save(ase_name,'SINR_CCDF');

%         for idx_test=1:length_Diff_UE_density
%             for idx=1:length(SINR_Threshold_array_dB)
%                 SINR_CCDF(idx_test,idx)=length(find(SINR(idx_test,:)>Threshold(idx)))/exp_num;
%                 SINR_CCDF_dB(idx_test,idx)=length(find(SINR_dB(idx_test,:)>SINR_Threshold_array_dB(idx)))/exp_num;
%             end
%         end
% %         figure(1);
%        
%         xlabel('SINR Threshold in dB');
%         ylabel('CCDF');
%         plot(SINR_Threshold_array_dB,SINR_CCDF(1,:),'r',SINR_Threshold_array_dB,SINR_CCDF(2,:),'g',SINR_Threshold_array_dB,SINR_CCDF(3,:),'b',SINR_Threshold_array_dB,SINR_CCDF(4,:),'k');
%         %set(gca,'YLim',[0 1]);
%         set(gca,'YTick',0:0.1:1);
%         set(gca,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'}); 
%          grid on;
%         box on;
%         hold on;
% SINR_NLoS=SINR_NLoS(SINR_NLoS>0);
ccdf_NLoS=zeros(1,length(Threshold));
ccdf_LoS=zeros(1,length(Threshold));
for idx=1:length(SINR_Threshold_array_dB)
    num_consideredUE_NLoS=length(SINR_NLoS(SINR_NLoS>0));
    num_consideredUE_LoS=length(SINR_LoS(SINR_LoS>0));
    num_qualifiedUE_NLoS=0;
    num_qualifiedUE_LoS=0;
    for idx_exp=1:exp_num
        num_qualifiedUE_LoS=num_qualifiedUE_LoS+(SINR_LoS(idx_exp)>Threshold(idx));
        num_qualifiedUE_NLoS=num_qualifiedUE_NLoS+(SINR_NLoS(idx_exp)>Threshold(idx));
    end
    ccdf_NLoS(idx)=num_qualifiedUE_NLoS/num_consideredUE_NLoS;
    ccdf_LoS(idx)=num_qualifiedUE_LoS/num_consideredUE_LoS;
end
hold on;
plot(SINR_Threshold_array_dB,ccdf_NLoS,'k',SINR_Threshold_array_dB,ccdf_LoS,'b');
box on ;
grid on; 
end
        