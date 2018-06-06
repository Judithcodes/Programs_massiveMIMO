% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: main_sim.m
% Authors: Xuefeng Yao
% Version: 4.0
% Date: 2017-06-01
% Description: Main function of the simulation for the
% modified HPPP model considering LoS and NLoS transmissions (UAS1)
% with the assumption of dynamic cell on-off
% Copyright(c): For personal study only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Add multiple random seeds in the following array to generate multiple MAT files in order not to cause memory overflow
rand_seed_bank = [601];

for rand_seed_idx = 1:length(rand_seed_bank)
    % close all;
    clearvars -global;
    clearvars -except rand_seed*;
    dbstop if error;
    
    % Load the random seed
    rand_seed = rand_seed_bank(rand_seed_idx);
    rand('state',rand_seed);
    randn('state',rand_seed);
    
    % Global variables
    
    %     global R1;                  % Constant R_1 in the LoS probability function in [m]
    %     global R2;                  % Constant R_2 in the LoS probability function in [m]
    global d0;                  % Constant d_0 in the LoS probability function in [m]
    %     global a;                   % Constant a in the LoS probability function in [m]
    global alpha_LoS;           % The path loss exponent for the LoS transmission
    global alpha_NLoS;          % The path loss exponent for the NLoS transmission
    global A_LoS;               % The path loss at r=1 for the LoS transmission
    global A_NLoS;              % The path loss at r=1 for the NLoS transmission
    % I changed the TX power to a BS density specific value on June 1st, 2015!!!
    % global P_tx;                % The transmission power of a BS in [W]
    global N0;                  % The AWGN noise power in [W]
    global UE_noise_figure;     % The UE's noise figure in [dB]
    global SINR_threshold;      % The SINR threshold for the definition of coverage
    
    
    % Set the simulation scheme and the number of experiments
    SCH = 'Prop';                       % [Prop]: The proposed scheme. [Jeff]: The existing scheme (J. G. Andrews, TCOM2011)
    sub_scheme = 'smallestPL';          	% 'nearestBS' or the subscheme of Jeff's work, 'lower' or 'upper'
    LoS_Pr_model = 'Linear';            % 'Linear' or '36.828' or '36.814' or 'Stepfn'
    small_cell_type = 'Xsparse';         % the type of the considered small cell network
    variable_Tx_power_enabler = 0;      % enabler of variable Tx power with respect to BS density
    multi_path_fading_type = 'Rayleigh';  % multi-path fading type: 'Rayleigh' or 'Rician'
    exp_num = 1e4;                     	% The number of experiments
    epslong = 1e-6;                     % A very small constant number
    figure_plot_enabler = 0;            % enabler of plotting figures
    
    % Set the BS density in BSs/m^2
    if strcmpi(small_cell_type, 'debug')
        BS_density_step = 0.5;
    else
        BS_density_step = 0.1;
    end
    
    if strcmpi(small_cell_type, 'XXXsparse')
        BS_density_array = 10.^[-8:BS_density_step:-7].';
    elseif strcmpi(small_cell_type, 'XXsparse')
        BS_density_array = 10.^[-7:BS_density_step:-6].';
    elseif strcmpi(small_cell_type, 'Xsparse')
        BS_density_array = 10.^[-6:BS_density_step:-3.1].';
    elseif strcmpi(small_cell_type, 'sparse')
        BS_density_array = 10.^[-3.1:BS_density_step:-2.4].';
    elseif strcmpi(small_cell_type, 'medium')
        BS_density_array = 10.^[-2.4:BS_density_step:-2.1].';
    elseif strcmpi(small_cell_type, 'dense')
        BS_density_array = 10.^[-2.1:BS_density_step:-2].';
    elseif strcmpi(small_cell_type, 'Xdense')
        BS_density_array = 10.^[-2:BS_density_step:-1].';
    elseif strcmpi(small_cell_type, 'XXdense')
        BS_density_array = 10.^[-1:BS_density_step:1].';
    elseif strcmpi(small_cell_type, 'debug')
        BS_density_array = 10.^[-6:BS_density_step:-2].';
    else
        error('undefined "small_cell_type"');
    end
    len_BS_density_array = length(BS_density_array);
    
    % Set the UE density in UEs/m^2
    UE_density = Inf;
    
    % Get the probability of turning on a small cell
    b_value = 3.5;
    q_value = 3.5;
    Pr_cell_off = ...
        (b_value*BS_density_array).^q_value ./ ...
        (b_value*BS_density_array+UE_density).^q_value;      	% The probability of cell turning off
    activated_BS_density_array = ...
        BS_density_array .* (1-Pr_cell_off);                 	% The average number of active small cells per km^2
    
    
    
    % ============================
    % Basic Parameter Setting
    % ============================
    
    alpha_LoS = 2.09;
    alpha_NLoS = 3.75;
    if strcmpi(SCH, 'Jeff') == 1 && strcmpi(sub_scheme, 'lower') == 1
        alpha_NLoS = 2.09
    end
    A_LoS = 10^(-4.11);
    A_NLoS = 10^(-3.29);
    %P_tx = 10^(24/10)*1e-3;
    UE_noise_figure = 9;
    N0 = 10^((-104 + UE_noise_figure)/10)*1e-3;
    SINR_threshold = 1;
    
    if variable_Tx_power_enabler == 1
        % Get the average size of a vorinoi cell [m^2] and the TX powers
        targe_SNR_dB_for_edge_UE = 15; % targe SNR for cell-edge UE in dB
        P_tx = NaN(1, len_BS_density_array);
        P_tx_dBm = NaN(1, len_BS_density_array);    % in dBm
        P_tx_max = 10^(46/10)*1/1e3;                % in Watts
        for BS_density_idx = 1:len_BS_density_array
            BS_density = BS_density_array(BS_density_idx);
            Voronoi_cell_ave_size = q_value/b_value/BS_density;
            Voronoi_cell_equi_radius = sqrt(Voronoi_cell_ave_size/pi);
            % calculate the average path gain at the cell edge
            %         average_path_gain = ...
            %             get_LoS_Pr_v3(SCH, LoS_Pr_model, Voronoi_cell_equi_radius) * ...
            %             A_LoS * Voronoi_cell_equi_radius^(-1*alpha_LoS) + ...
            %             (1-get_LoS_Pr_v3(SCH, LoS_Pr_model, Voronoi_cell_equi_radius)) * ...
            %             A_NLoS * Voronoi_cell_equi_radius^(-1*alpha_NLoS);
            
            average_path_gain = ...
                A_NLoS * Voronoi_cell_equi_radius^(-1*alpha_NLoS);
            
            P_tx(BS_density_idx) = 10^(targe_SNR_dB_for_edge_UE/10) * N0 / average_path_gain;
            if P_tx(BS_density_idx) > P_tx_max
                P_tx(BS_density_idx) = P_tx_max;
            end
            P_tx_dBm(BS_density_idx) = 10*log10(P_tx(BS_density_idx)*1000);
        end
    else
        for BS_density_idx = 1:len_BS_density_array
            P_tx_dBm(BS_density_idx) = 24;    % in dBm
            P_tx(BS_density_idx) = 10^(P_tx_dBm(BS_density_idx)/10)*1/1e3;    % in Watts
        end
    end
    
    dx_dB = 0.25;
    if strcmpi(small_cell_type, 'debug')
        SINR_threshold_array_dB = [0]; % for debug purposes
    else
        SINR_threshold_array_dB = [-60:dx_dB:60];
    end
    SINR_threshold_array = 10.^(SINR_threshold_array_dB/10);
    
    len_SINR_threshold_array = length(SINR_threshold_array);
    dx_array = SINR_threshold_array(2:end)-SINR_threshold_array(1:end-1);
    
    % numerical range of r
    min_r = 1e-3;
    enlarged_max_r_factor = 2.5;
    
    if strcmpi(small_cell_type, 'XXdense')
        max_r = 1e3
        dr = 0.05
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'Xdense')
        max_r = 1e3
        dr = 0.1
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'dense')
        max_r = 2.5e3
        dr = 0.25
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'medium')
        max_r = 2.5e3
        dr = 0.25
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'sparse')
        max_r = 2.5e3
        dr = 0.25
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'Xsparse')
        max_r = 5e3
        dr = 0.5
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'XXsparse')
        max_r = 1e4
        dr = 1
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'XXXsparse')
        max_r = 2.5e4
        dr = 2.5
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    elseif strcmpi(small_cell_type, 'debug')
        max_r = 5e4
        dr = 5
        min_average_BS_num = BS_density_array(1) * 4 * max_r^2
        max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    else
        error('undefined "small_cell_type"');
    end
    
    r_array = [min_r:dr:max_r];
    len_r_array = length(r_array);
    
    % quick check of the LoS probability function
    LoS_Pr_r_array = get_LoS_Pr_v3(SCH, LoS_Pr_model, r_array);
    NLoS_Pr_r_array = 1 - LoS_Pr_r_array;
    if figure_plot_enabler == 1
        figure(100)
        plot(r_array, LoS_Pr_r_array);
        axis([0, 1000, 0, 1]);
    end
    
    if variable_Tx_power_enabler == 1
        str_Tx_power = 'varPower';
    else
        str_Tx_power = ['fixP' num2str(P_tx_dBm(1)) 'dBm'];
    end
    
    % Naming the MAT files
    MATname_scenario = ['[Scenario]' 'd0' num2str(d0)];
    % MATname_data = ['[SCH_' SCH ']' 'd0_' num2str(round(d0)) '_LoSMod_' LoS_Pr_model([1:2 end-2:end]) '_RS' num2str(rand_seed)]
    MATname_data = ['[Sim_SCH_' SCH '_' sub_scheme ']' '_LoSMod_' LoS_Pr_model([1:2 end-2:end]) ...
        '_' str_Tx_power ...
        '_' multi_path_fading_type ...
        '_' num2str(len_SINR_threshold_array) 'SINRs' ...
        '_BSdensity_' small_cell_type '_UEdensity_' num2str(UE_density*1e6) 'UEperSqkm_RS' num2str(rand_seed)]
    MATname_brief_data = ['[Brief]' MATname_data]
    
    
    % initialize the registers for results
    %p_cov = NaN(len_BS_density_array, len_SINR_threshold_array);
    area_spec_efficiency = NaN(len_BS_density_array, 1);
    Count_in_cov = zeros(len_BS_density_array, len_SINR_threshold_array);
    SINR_full_result = NaN(len_BS_density_array, exp_num);
    Dyn_on_BSs_per_square_km = NaN(len_BS_density_array, exp_num);
    ave_dyn_on_BSs_per_square_km = NaN(len_BS_density_array, 1);
    effective_BS_num = NaN(len_BS_density_array, exp_num);
    
    tstart=tic;
    for BS_density_idx = 1:len_BS_density_array
        BS_density = BS_density_array(BS_density_idx)
        fprintf('\nThe %d-th BS density of the total %d BS densities......', BS_density_idx, len_BS_density_array);
        
        for exp_idx = 1:exp_num
            
            % Report the simulation progress
            if mod(exp_idx, exp_num/100) == 0
                fprintf('\nSimulation progress: %3.1f%% completed......', (exp_idx/exp_num*100));
            end
            
            considered_area_size = (2*max_r)^2;                             % The size of considered area in m^2
            enlarged_considered_area_size = (enlarged_max_r_factor*max_r)^2;              	% The size of considered area in m^2
            
            total_BS_num = 0;
            while total_BS_num == 0,
                total_BS_num = poissrnd(BS_density * enlarged_considered_area_size);   	% The number of total BSs in a "2*max_r" by "2*max_r" square
            end
            BS_cpos = unifrnd(-enlarged_max_r_factor/2*max_r,enlarged_max_r_factor/2*max_r, 1, total_BS_num) ...
                + 1i * unifrnd(-enlarged_max_r_factor/2*max_r,enlarged_max_r_factor/2*max_r, 1, total_BS_num);            % Complex coordinates of all BSs
            
            if UE_density < 1e6
                
                total_UE_num = 0;
                while total_UE_num < 2,
                    total_UE_num = poissrnd(UE_density * enlarged_considered_area_size);   	% The number of total BSs in a "2*max_r" by "2*max_r" square
                end
                % be careful! the UE at the origin could activate another BS!
                UE_cpos = [0 ...
                    unifrnd(-enlarged_max_r_factor/2*max_r,enlarged_max_r_factor/2*max_r, 1, total_UE_num-1)+ ...
                    1i * unifrnd(-enlarged_max_r_factor/2*max_r,enlarged_max_r_factor/2*max_r, 1, total_UE_num-1)];            % Complex coordinates of all UEs
                
                % Check which BSs in the enlarged considered area should be powered on
                BS_on_indicator_array = zeros(1, total_BS_num);
                for UE_idx = 1:total_UE_num
                    UE_to_BS_distance_temp = abs(UE_cpos(UE_idx) - BS_cpos);
                    LoS_probability_array = get_LoS_Pr_v3(SCH, LoS_Pr_model, UE_to_BS_distance_temp);
                    if UE_idx == 1
                        LoS_probability_array_for_UE1 = LoS_probability_array;
                    end
                    
                    long_term_RSRP = NaN(1, total_BS_num);
                    for BS_idx = 1:total_BS_num
                        if rand(1) < LoS_probability_array(BS_idx)
                            long_term_RSRP(BS_idx) = ...
                                P_tx(BS_density_idx) * A_LoS * UE_to_BS_distance_temp(BS_idx)^(-1*alpha_LoS);
                        else
                            long_term_RSRP(BS_idx) = ...
                                P_tx(BS_density_idx) * A_NLoS * UE_to_BS_distance_temp(BS_idx)^(-1*alpha_NLoS);
                        end
                    end
                    if UE_idx == 1
                        long_term_RSRP_for_UE1 = long_term_RSRP;
                    end
                    
                    [dummy_max_RSRP, association_BS_idx] = max(long_term_RSRP);
                    BS_on_indicator_array(association_BS_idx) = 1;
                    
                    if sum(BS_on_indicator_array) == total_BS_num
                        break;
                    end
                end
                
                % Check which BSs in the considered area have been powered on
                effective_BS_on_indicator_array = zeros(1, total_BS_num);
                effective_BS_indices = [];                                      % Indices of effective BSs
                for BS_idx = 1:total_BS_num
                    if BS_on_indicator_array(BS_idx) == 1 ...
                            && abs(real(BS_cpos(BS_idx))) <= max_r ...
                            && abs(imag(BS_cpos(BS_idx))) <= max_r
                        effective_BS_on_indicator_array(BS_idx) = 1;
                        effective_BS_indices = [effective_BS_indices BS_idx];
                    end
                end
                
                % Record the number of dynamically powered-on BSs per square km
                Dyn_on_BSs_per_square_km(BS_density_idx, exp_idx) = ...
                    sum(effective_BS_on_indicator_array)/(considered_area_size/1e6);
            else % infinite UEs, thus all BSs are powered on
                effective_BS_indices = 1:total_BS_num;
            end
            
            effective_BS_num(BS_density_idx, exp_idx) = length(effective_BS_indices);
            effective_r = abs(BS_cpos(effective_BS_indices));   % Distances r of all BSs
            
            
            % Generate the signal strengths from all BSs to the
            % interested UE
            % [sorted_effective_r, sorted_BS_indices] = sort(effective_r, 'ascend');
            
            if UE_density < 1e6
                LoS_probability_array = LoS_probability_array_for_UE1(effective_BS_indices);
                %LoS_probability_array = sorted_LoS_probability_array_for_UE1(sorted_BS_indices);
                
                sub_sample_long_term_RSRP_for_UE1 = long_term_RSRP_for_UE1(effective_BS_indices);
                %sorted_long_term_RSRP_for_UE1 = sorted_long_term_RSRP_for_UE1(sorted_BS_indices);
            else
                LoS_probability_array = get_LoS_Pr_v3(SCH, LoS_Pr_model, effective_r);
            end
            
            %h = zeros(1,effective_BS_num(BS_density_idx, exp_idx));
            %h = exprnd(1,1,effective_BS_num(BS_density_idx, exp_idx));
            long_term_RSRP = NaN(1, effective_BS_num(BS_density_idx, exp_idx));
            short_term_RSRP = NaN(1, effective_BS_num(BS_density_idx, exp_idx));
            
            for effective_BS_idx = 1:effective_BS_num(BS_density_idx, exp_idx)
                if rand(1) < LoS_probability_array(effective_BS_idx) % LoS path
                    long_term_RSRP(effective_BS_idx) = ...
                        P_tx(BS_density_idx) * A_LoS * effective_r(effective_BS_idx)^(-1*alpha_LoS);
                    if strcmpi(multi_path_fading_type, 'Rayleigh') == 1
                        h = exprnd(1,1,1);
                    elseif strcmpi(multi_path_fading_type, 'Rician') == 1
                        % K_factor = 32;
                        K_factor_dB = 13 - 0.03 * effective_r(effective_BS_idx); % 3GPP & 3GPP2 SCM model, distance-dependent K
                        K_factor = 10^(K_factor_dB/10);
                        channel_coeff = sqrt(K_factor/(K_factor+1)) + ...
                            sqrt(1/(K_factor+1))*(1/sqrt(2))*(randn(1) + 1i*randn(1));
                        h = abs(channel_coeff)^2;
                    else
                        error('unknown multi-path fading!');
                    end
                else % NLoS path
                    long_term_RSRP(effective_BS_idx) = ...
                        P_tx(BS_density_idx) * A_NLoS * effective_r(effective_BS_idx)^(-1*alpha_NLoS);
                    h = exprnd(1,1,1);
                end
                short_term_RSRP(effective_BS_idx) = long_term_RSRP(effective_BS_idx) * h;
            end
            if UE_density < 1e6
                if sum(abs(long_term_RSRP - sub_sample_long_term_RSRP_for_UE1)) > 1e-3
                    error('unknown error with the long-term RSRP');
                end
            end
            [sorted_long_term_RSRP, sorted_BS_indices_based_on_long_term_RSRP] = sort(long_term_RSRP, 'descend');
            
            sorted_short_term_RSRP = short_term_RSRP(sorted_BS_indices_based_on_long_term_RSRP);
            
            SINR_full_result(BS_density_idx, exp_idx) = ...
                sorted_short_term_RSRP(1)/(sum(sorted_short_term_RSRP)-sorted_short_term_RSRP(1)+N0);
            
            % Count the times of in-coverage for various SINR thresholds
            for SINR_threshold_idx = 1:len_SINR_threshold_array
                SINR_threshold = SINR_threshold_array(SINR_threshold_idx);
                if SINR_full_result(BS_density_idx, exp_idx) > SINR_threshold
                    Count_in_cov(BS_density_idx, SINR_threshold_idx) = ...
                        Count_in_cov(BS_density_idx, SINR_threshold_idx) + 1;
                end
            end % end of SINR loop
            
        end % end of experiment loop
        ave_dyn_on_BSs_per_square_km(BS_density_idx) = mean(Dyn_on_BSs_per_square_km(BS_density_idx, :));
        % calculate the coverage probability
        p_cov = Count_in_cov/exp_num;
        save(MATname_data);
    end % end of BS density loop
    telapsed_hours = toc(tstart)/3600
    
    
    BS_cpos = [];
    UE_cpos = [];
    LoS_probability_array = [];
    effective_BS_indices = [];
    effective_r = [];
    long_term_RSRP = [];
    short_term_RSRP = [];
    sorted_BS_indices_based_on_long_term_RSRP = [];
    sorted_effective_r = [];
    sorted_long_term_RSRP = [];
    sorted_short_term_RSRP = [];
    SINR_full_result = [];
    save(MATname_brief_data);
    
    % plot the figures
    
    QoS_SINR_th = 1;
    % coverage probability vs lambda
    if figure_plot_enabler == 1
        figure(1)
        box on;
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            semilogx(BS_density_array*1e6, p_cov(:, find(SINR_threshold_array >= QoS_SINR_th,1)), 'rx-.');
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            semilogx(BS_density_array*1e6, p_cov(:, find(SINR_threshold_array >= QoS_SINR_th,1)), 'bx-.');
            % semilogx(BS_density_array, p_cov(:, find(SINR_threshold_array_dB==0,1)), 'b');
        else
            error('undefined scheme');
        end
        hold on;
        xlabel('BS density [1/km^2] \it{\lambda}');
        ylabel('Probability of \it{SINR}\rm{>}\it{\gamma}');
        axis([min(BS_density_array)*1e6,max(BS_density_array)*1e6,0,1]);
        %axis([min(BS_density_array)*1e6,1e4,0,1]);
        
        % legend('Conventional HPPP', 'Modified HPPP considering LoS/NLoS');
    end
    
    % CCDF of SINR
    if figure_plot_enabler == 1
        figure(2)
        box on;
        xlabel('\it{SINR} \rm{[dB]}');
        ylabel('CCDF');
        hold on;
        for BS_density_idx = 1:len_BS_density_array
            if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
                plot(SINR_threshold_array_dB, p_cov(BS_density_idx, :), 'r');
            elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
                plot(SINR_threshold_array_dB, p_cov(BS_density_idx, :), 'b');
            else
                error('undefined scheme');
            end
        end
        
        % axis([min(SINR_threshold_array_dB),max(SINR_threshold_array_dB),0,1]);
        
        % legend('Conventional HPPP', 'Modified HPPP considering LoS/NLoS');
    end
    
    
    %     x_array = 1/2 * (SINR_threshold_array(1:end-1) + SINR_threshold_array(2:end));
    %     usable_x_idx = find(x_array >= QoS_SINR_th, 1);
    %     for BS_density_idx = 3:len_BS_density_array
    %         numerical_SINR_PDF = (1-p_cov(BS_density_idx,2:end)) - (1-p_cov(BS_density_idx,1:end-1));
    %
    %         %     figure(99)
    %         %     % plot(10*log10(x_array), numerical_SINR_PDF);
    %         %     plot(x_array, numerical_SINR_PDF);
    %
    %         area_spec_efficiency_exact(BS_density_idx) = ...
    %             sum(log2(1+x_array(usable_x_idx:end)) .* numerical_SINR_PDF(usable_x_idx:end) .* dx_array(usable_x_idx:end)) * BS_density_array(BS_density_idx);
    %         area_spec_efficiency_LB(BS_density_idx) = ...
    %             sum(log2(1+QoS_SINR_th) .* numerical_SINR_PDF(usable_x_idx:end) .* dx_array(usable_x_idx:end)) * BS_density_array(BS_density_idx);
    %         area_spec_efficiency_UB(BS_density_idx) = ...
    %             log2(1 + sum(x_array(usable_x_idx:end) .* numerical_SINR_PDF(usable_x_idx:end) .* dx_array(usable_x_idx:end))) * BS_density_array(BS_density_idx);
    %     end
    %
    %
    %     % coverage probability vs lambda
    %     figure(33)
    %     box on;
    %     if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
    %         semilogx(BS_density_array, area_spec_efficiency_exact, 'r');
    %         hold on;
    %         semilogx(BS_density_array, area_spec_efficiency_LB, 'k');
    %     elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
    %         % plot(BS_density_array, area_spec_efficiency_exact, 'b');
    %         semilogx(BS_density_array, area_spec_efficiency_exact, 'b');
    %         hold on;
    %         semilogx(BS_density_array, area_spec_efficiency_LB, 'k');
    %     else
    %         error('undefined scheme');
    %     end
    %
    %     xlabel('BS density [1/m^2] \it{\lambda}');
    %     ylabel('Area spectral efficiency [bps/Hz/m^2]');
    %
    %     legend('Exact value', 'Lower bound');
    %
    %
    %
    %     deri_gx = 1/2 * (1./(1+SINR_threshold_array(1:end-1)) + 1./(1+SINR_threshold_array(2:end)));
    %     for BS_density_idx = 1:len_BS_density_array
    %         F_bar_x = 1/2 * (p_cov(BS_density_idx,1:end-1) + p_cov(BS_density_idx,2:end));
    %         area_spec_efficiency(BS_density_idx) = ...
    %             sum(F_bar_x .* deri_gx .* dx_array) * BS_density_array(BS_density_idx);
    %     end
    %
    %     % coverage probability vs lambda
    %     figure(3)
    %     box on;
    %     if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
    %         plot(BS_density_array, area_spec_efficiency, 'r');
    %     elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
    %         plot(BS_density_array, area_spec_efficiency, 'b');
    %         % semilogx(BS_density_array, area_spec_efficiency, 'b');
    %     else
    %         error('undefined scheme');
    %     end
    %     hold on;
    %     xlabel('BS density \it{\lambda}');
    %     ylabel('Area spectral efficiency [bps/Hz/km^2]');
    %     % axis([min(BS_density_array),max(BS_density_array),0,1]);
    %
    %     % legend('Conventional HPPP', 'Modified HPPP considering LoS/NLoS');
    %
    %
    %
    %
    %     deri_gx = 1/2 * (1./(1+SINR_threshold_array(1:end-1)) + 1./(1+SINR_threshold_array(2:end)));
    %     for BS_density_idx = 1:len_BS_density_array
    %         F_bar_x = 1/2 * (p_cov(BS_density_idx,1:end-1) + p_cov(BS_density_idx,2:end));
    %         area_spec_efficiency(BS_density_idx) = ...
    %             sum(F_bar_x .* deri_gx .* dx_array);
    %     end
    %
    %     % coverage probability vs lambda
    %     figure(4)
    %     box on;
    %     if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
    %         semilogx(BS_density_array, area_spec_efficiency, 'r');
    %     elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
    %         semilogx(BS_density_array, area_spec_efficiency, 'b');
    %         % semilogx(BS_density_array, area_spec_efficiency, 'b');
    %     else
    %         error('undefined scheme');
    %     end
    %     hold on;
    %     xlabel('BS density \it{\lambda}');
    %     ylabel('Area spectral efficiency [bps/Hz/km^2]');
    %     % axis([min(BS_density_array),max(BS_density_array),0,1]);
    %
    %     % legend('Conventional HPPP', 'Modified HPPP considering LoS/NLoS');
    
end