% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: main.m
% Authors: Xuefeng Yao
% Version: 13.0
% Date: 2017-06-25
% Description: Main function of the study on network analysis based on
% modified HPPP model considering LoS and NLoS transmissions
% [UAS 1: the UE associated with the BS with the smallest path loss]
% Copyright(c): For personal study only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% Add multiple random seeds in the following array to generate multiple MAT files in order not to cause memory overflow
rand_seed_bank = [1300];

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
    
    global R1;                  % Constant R_1 in the LoS probability function in [m]
    global R2;                  % Constant R_2 in the LoS probability function in [m]
    global d0;                  % Constant d_0 in the LoS probability function in [m]
    global a;                   % Constant a in the LoS probability function in [m]
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
    sub_scheme = 'smallestPL';        	% 'smallestPL' or the subscheme of Jeff's work, 'lower' or 'upper'
    LoS_Pr_model = 'Linear';            % 'Linear' or '36.828' or '36.814' or 'Stepfn' or '3psLin' or '2psRvS'
    small_cell_type = 'debug';          % the type of the considered small cell network: 'debug' or 'seeking_best_lambda' or 'xxx'
    variable_Tx_power_enabler = 0;      % enabler of variable Tx power with respect to BS density
    plot_enabler = 0;                   % the enabler of figure plotting in the program
    epslong = 1e-6;                     % A very small constant number
    
    % Set the BS density in BSs/m^2
    if strcmpi(small_cell_type, 'debug')
        BS_density_step = 0.1;
    elseif strcmpi(small_cell_type, 'seeking_best_lambda')
        BS_density_step = 0.01;
    else
        BS_density_step = 0.1;
    end
    
    if strcmpi(small_cell_type, 'XXsparse')
        BS_density_array = 10.^[-10:BS_density_step:-7].';
    elseif strcmpi(small_cell_type, 'Xsparse')
        BS_density_array = 10.^[-7:BS_density_step:-4].';
    elseif strcmpi(small_cell_type, 'sparse')
        BS_density_array = 10.^[-4:BS_density_step:-2.5].';
    elseif strcmpi(small_cell_type, 'medium')
        BS_density_array = 10.^[-2.5:BS_density_step:-1.8].';
    elseif strcmpi(small_cell_type, 'dense')
        BS_density_array = 10.^[-1.8:BS_density_step:-1.3].';
    elseif strcmpi(small_cell_type, 'Xdense')
        BS_density_array = 10.^[-1.3:BS_density_step:-1].';
    elseif strcmpi(small_cell_type, 'XXdense')
        BS_density_array = 10.^[-1:BS_density_step:1].';
    elseif strcmpi(small_cell_type, 'debug')
        %BS_density_array = 10.^[-11:BS_density_step:-1].';
        BS_density_array = 10.^[-6:BS_density_step:-2].';
    elseif strcmpi(small_cell_type, 'seeking_best_lambda')
        if strcmpi(LoS_Pr_model, 'Linear') == 1
            BS_density_array = [10:BS_density_step:25]./1e6;
        elseif strcmpi(LoS_Pr_model, '36.828') == 1
            BS_density_array = [70:BS_density_step:110]./1e6;
        else
            error('unknown path loss model!');
        end
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
    A_LoS = 10^(-1*(10.38-2.09*3));
    A_NLoS = 10^(-1*(14.54-3.75*3));
    
    if strcmpi(LoS_Pr_model, '3psLin') == 1
        d0 = 300;
        d1 = 900;
        alpha_LoS = 2.09;
        alpha_NLoS = 3.75;
        alpha_NLoS_add1 = 4.75;
        %A_NLoS*d0^(-1*alpha_NLoS)
        %A_NLoS_add1*d0^(-1*alpha_NLoS_add1)
        A_NLoS_add1 = A_NLoS*d0^(alpha_NLoS_add1 - alpha_NLoS);
        alpha_NLoS_add2 = 5.75;
        %A_NLoS_add1*d1^(-1*alpha_NLoS_add1)
        %A_NLoS_add2*d1^(-1*alpha_NLoS_add2)
        A_NLoS_add2 = A_NLoS_add1*d1^(alpha_NLoS_add2 - alpha_NLoS_add1);
        
        test_r_array = 1:1:5e3;
        pl_nlos_linear = zeros(1, length(test_r_array));
        pl_los_linear = zeros(1, length(test_r_array));
        for test_r_idx = 1:length(test_r_array)
            test_r = test_r_array(test_r_idx);
            if test_r <= d0
                pl_nlos_linear(test_r_idx) = A_NLoS*test_r^(-1*alpha_NLoS);
            elseif test_r > d0 && test_r <= d1
                pl_nlos_linear(test_r_idx) = A_NLoS_add1*test_r^(-1*alpha_NLoS_add1);
            elseif test_r > d1
                pl_nlos_linear(test_r_idx) = A_NLoS_add2*test_r^(-1*alpha_NLoS_add2);
            end
            pl_los_linear(test_r_idx) = A_LoS*test_r^(-1*alpha_LoS);
        end
    else
        test_r_array = 1:1:5e3;
        pl_nlos_linear = A_NLoS*test_r_array.^(-1*alpha_NLoS);
        pl_los_linear = A_LoS*test_r_array.^(-1*alpha_LoS);
    end
    if plot_enabler
        figure(10)
        semilogx(test_r_array, -10*log10(pl_nlos_linear));
        hold on;
        semilogx(test_r_array, -10*log10(pl_los_linear), 'r');
        xlabel('Distance [m]');
        ylabel('Pathloss [dB]');
        legend('NLoS', 'LoS');
    end
    
    %P_tx = 10^(24/10)*1e-3;
    UE_noise_figure = 9; % in dB
    N0 = 10^((-104 + UE_noise_figure)/10)*1e-3;% 10 MHz!
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
            %             get_LoS_Pr_v2(SCH, LoS_Pr_model, Voronoi_cell_equi_radius) * ...
            %             A_LoS * Voronoi_cell_equi_radius^(-1*alpha_LoS) + ...
            %             (1-get_LoS_Pr_v2(SCH, LoS_Pr_model, Voronoi_cell_equi_radius)) * ...
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
    
    % plot the power vs lambda
    if plot_enabler
        figure(20)
        semilogx(BS_density_array*1e6, P_tx_dBm, 'k');
        box on;
        axis([10^-4, 10^4, -20, 60]);
        xlabel('BS density \it{\lambda}  \rm{[BSs/km^2]}');
        ylabel('BS Tx power [dBm]');
    end
    
    dx_dB = 0.25;
    if strcmpi(small_cell_type, 'debug')
        SINR_threshold_array_dB = [0]; % for debug purposes
    elseif strcmpi(small_cell_type, 'seeking_best_lambda')
        SINR_threshold_array_dB = [-3 0 3 6 10]; % for 'seeking_best_lambda'
    else
        SINR_threshold_array_dB = [-60:dx_dB:60];
    end
    SINR_threshold_array = 10.^(SINR_threshold_array_dB/10);
    
    len_SINR_threshold_array = length(SINR_threshold_array);
    dx_array = SINR_threshold_array(2:end)-SINR_threshold_array(1:end-1);
    
    % numerical range of r
    min_r = 1e-3;
    
    if strcmpi(small_cell_type, 'XXsparse')
        max_r = 1e5
        dr = 10
    elseif strcmpi(small_cell_type, 'Xsparse')
        max_r = 1e4
        dr = 1
    elseif strcmpi(small_cell_type, 'sparse')
        max_r = 5e3
        dr = 0.5
    elseif strcmpi(small_cell_type, 'medium')
        max_r = 2.5e3
        dr = 0.25
    elseif strcmpi(small_cell_type, 'dense')
        max_r = 1e3
        dr = 0.1
    elseif strcmpi(small_cell_type, 'Xdense')
        max_r = 1e3
        dr = 0.1
    elseif strcmpi(small_cell_type, 'XXdense')
        max_r = 1e3
        dr = 0.05
    elseif strcmpi(small_cell_type, 'debug')
        max_r = 1e4
        dr = 1
    elseif strcmpi(small_cell_type, 'seeking_best_lambda')
        max_r = 1e4
        dr = 1
    else
        error('undefined "small_cell_type"');
    end
    min_average_BS_num = BS_density_array(1) * 4 * max_r^2
    max_average_BS_num = BS_density_array(end) * 4 * max_r^2
    
    r_array = [min_r:dr:max_r];
    len_r_array = length(r_array);
    
    % quick check of the LoS probability function
    LoS_Pr_r_array = get_LoS_Pr_v3(SCH, LoS_Pr_model, r_array);
    NLoS_Pr_r_array = 1 - LoS_Pr_r_array;
    if plot_enabler
        if strcmpi(LoS_Pr_model, '2psRvS') == 1
            LoS_Pr_r_array_reverseS = get_LoS_Pr_v3(SCH, '36.828', r_array);
            figure(100)
            xlabel('Distance \it{r} \rm{[km]}');
            ylabel('\rm{Pr}^{\rm{ L}}\rm{(}\it{r}\rm{)}');
            hold on;
            box on;
            plot(r_array/1e3, LoS_Pr_r_array_reverseS, 'k-');
            plot(r_array/1e3, LoS_Pr_r_array, 'k-.');
            axis([0, 300/1e3, 0, 1.1]);
            legend('\rm{Pr}^{\rm{L}}\rm{(}\it{r}\rm{)} of 3GPP Case 2', ...
                'Approximation using piece-wise linear functions');

        else
            figure(100)
            hold on;
            box on;
            plot(r_array/1e3, LoS_Pr_r_array, 'k');
            axis([0, 500/1e3, 0, 1]);
        end
    end
    
    if variable_Tx_power_enabler == 1
        str_Tx_power = 'varPower';
    else
        str_Tx_power = ['fixP' num2str(P_tx_dBm(1)) 'dBm'];
    end
    
    % Naming the MAT files
    MATname_scenario = ['[Scenario]' 'd0' num2str(d0)];
    % MATname_data = ['[SCH_' SCH ']' 'd0_' num2str(round(d0)) '_LoSMod_' LoS_Pr_model([1:2 end-2:end]) '_RS' num2str(rand_seed)]
    MATname_data = ['[SCH_' SCH '_' sub_scheme ']' '_LoSMod_' LoS_Pr_model([1:2 end-2:end]) ...
        '_' str_Tx_power ...
        '_' num2str(len_SINR_threshold_array) 'SINRs' ...
        '_BSdensity_' small_cell_type '_UEdensity_' num2str(UE_density*1e6) 'UEperSqkm_RS' num2str(rand_seed)]
    MATname_brief_data = ['[Brief]' MATname_data]
    
    
    x1 = d0^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS))
    
    if d0 < r_array(1)
        r_lower_indices = [];
        r_middle_indices = [];
        r_upper_indices = find(r_array >= d0);
    else
        r_lower_indices = find(r_array < x1);
        r_middle_indices = find(r_array([r_lower_indices(end)+1:end]) < d0) + r_lower_indices(end);
        r_upper_indices = find(r_array >= d0);
    end
    
    %     if strcmpi(LoS_Pr_model, '3psLin') == 1
    %         r_lower_indices = find(r_array < x1);
    %         r_middle_indices = find(r_array([r_lower_indices(end)+1:end]) < d0) + r_lower_indices(end);
    %         r_upper_indices = find(r_array([r_middle_indices(end)+1:end]) < d1) + r_middle_indices(end);
    %         r_Xupper_indices = find(r_array >= d1);
    %     end
    
    % initialize the registers for results
    p_cov = NaN(len_BS_density_array, len_SINR_threshold_array);
    area_spec_efficiency = NaN(len_BS_density_array, 1);
    
    intermediate_results = repmat({NaN(len_SINR_threshold_array, len_r_array)}, 1, 12);
    
    %     for BS_density_idx = 1:len_BS_density_array
    %         BS_density = BS_density_array(BS_density_idx)
    geom_density_term_LoS = NaN(len_BS_density_array, len_r_array);
    geom_density_term_NLoS = NaN(len_BS_density_array, len_r_array);
    
    du = dr;
    u_array = r_array;
    d0_idx = find(u_array >= d0, 1);
    
    for SINR_threshold_idx = 1:len_SINR_threshold_array
        SINR_threshold = SINR_threshold_array(SINR_threshold_idx)
        fprintf('The %d-th SINR point of totally %d SINR points...\n', SINR_threshold_idx, len_SINR_threshold_array);
        
        Laplace_term_LoS = NaN(len_BS_density_array, len_r_array);
        Laplace_term_NLoS = NaN(len_BS_density_array, len_r_array);
        T_cov = NaN(len_BS_density_array,6);
        
        for r_idx = 1:len_r_array
            
            % Report the simulation progress
            if mod(r_idx, len_r_array/100) == 0
                fprintf('\nSimulation progress: %3d%% completed......', floor(r_idx/len_r_array*100));
            end
            
            r = r_array(r_idx);
            
            if r < x1
                if strcmpi(LoS_Pr_model, '36.828') == 1 ...
                        || strcmpi(LoS_Pr_model, 'Stepfn') == 1 ...
                        || strcmpi(LoS_Pr_model, '2psRvS') == 1
                    %                     u_array_case1_lw = [r:du:d0-epslong];
                    %                     u_array_case1_up = [d0+epslong:du:max_r];
                    u_array_case1_lw = u_array([r_idx:d0_idx]);
                    u_array_case1_up = u_array([d0_idx+1:end]);
                    
                    r1 = r^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
                    r2 = r^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
                    r1_idx = find(u_array >= r1, 1);
                    r2_idx = find(u_array >= r2, 1);
                    
                    if SINR_threshold_idx == 1
                        geom_density_term_LoS(:, r_idx) = ...
                            -2*pi*BS_density_array.* sum(LoS_Pr_r_array([1:r_idx]) .* ...
                            u_array([1:r_idx]) .* du)-2*pi*BS_density_array .* ...
                            sum(NLoS_Pr_r_array([1:r1_idx]) .* u_array([1:r1_idx]) .* du);
                        
                        geom_density_term_NLoS(:, r_idx) = ...
                            -2*pi*BS_density_array.* sum(LoS_Pr_r_array([1:r2_idx]) .* ...
                            u_array([1:r2_idx]) .* du) - 2*pi*BS_density_array .* ...
                            sum(NLoS_Pr_r_array([1:r_idx]) .* u_array([1:r_idx]) .* du);
                    end
                    
                    intermediate_results{1}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + (SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_lw.^alpha_LoS) * du);
                    intermediate_results{2}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_lw.^alpha_NLoS) * du);
                    intermediate_results{3}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + (SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_up.^alpha_LoS) * du);
                    intermediate_results{4}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_up.^alpha_NLoS) * du);
                    
                    Laplace_term_LoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    
                    intermediate_results{5}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_lw.^alpha_LoS) * du);
                    intermediate_results{6}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + (SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_lw.^alpha_NLoS) * du);
                    intermediate_results{7}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_up.^alpha_LoS) * du);
                    intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + (SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_up.^alpha_NLoS) * du);
                    
                    Laplace_term_NLoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                elseif strcmpi(LoS_Pr_model, 'Linear') == 1 || strcmpi(LoS_Pr_model, '3psLin') == 1
                    
                    r1 = r^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
                    r2 = r^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
                    
                    if SINR_threshold_idx == 1
                        geom_density_term_LoS(:, r_idx) = ...
                            -1*pi*BS_density_array.*r.^2 + 2*pi*BS_density_array./(3*d0)*(r.^3 - r1.^3);
                        
                        geom_density_term_NLoS(:, r_idx) = ...
                            -1*pi*BS_density_array.*r2.^2 + 2*pi*BS_density_array./(3*d0)*(r2.^3-r.^3);
                    end
                    
                    intermediate_results{1}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_LoS, 1, (SINR_threshold*r^alpha_LoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_LoS, 1, (SINR_threshold*r^alpha_LoS)^(-1), r));
                    intermediate_results{2}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_LoS, 2, (SINR_threshold*r^alpha_LoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_LoS, 2, (SINR_threshold*r^alpha_LoS)^(-1), r));
                    intermediate_results{3}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_NLoS, 2, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_NLoS, 2, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), r1));
                    if strcmpi(LoS_Pr_model, '3psLin') == 1
                        intermediate_results{4}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(1, alpha_NLoS_add1, 1, A_LoS/A_NLoS_add1*(SINR_threshold*r^alpha_LoS)^(-1), d1) - ...
                            get_rou_value(1, alpha_NLoS_add1, 1, A_LoS/A_NLoS_add1*(SINR_threshold*r^alpha_LoS)^(-1), d0)) + ...
                            (get_rou_value(2, alpha_NLoS_add2, 1, A_LoS/A_NLoS_add2*(SINR_threshold*r^alpha_LoS)^(-1), d1));
                    else
                        intermediate_results{4}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(2, alpha_NLoS, 1, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), d0));
                    end
                    
                    Laplace_term_LoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                        2*pi*activated_BS_density_array/d0 * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array/d0 * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                    %                             2*pi*activated_BS_density_array/d0 * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array/d0 * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    
                    %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(1, alpha_NLoS, 1, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), max_r) - ...
                    %                         get_rou_value(1, alpha_NLoS, 1, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), d0)));
                    
                    
                    
                    
                    intermediate_results{5}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_LoS, 1, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_LoS, 1, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), r2));
                    intermediate_results{6}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_LoS, 2, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_LoS, 2, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), r2));
                    intermediate_results{7}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_NLoS, 2, (SINR_threshold*r^alpha_NLoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_NLoS, 2, (SINR_threshold*r^alpha_NLoS)^(-1), r));
                    if strcmpi(LoS_Pr_model, '3psLin') == 1
                        intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(1, alpha_NLoS_add1, 1, A_NLoS/A_NLoS_add1*(SINR_threshold*r^alpha_NLoS)^(-1), d1) - ...
                            get_rou_value(1, alpha_NLoS_add1, 1, A_NLoS/A_NLoS_add1*(SINR_threshold*r^alpha_NLoS)^(-1), d0)) + ...
                            (get_rou_value(2, alpha_NLoS_add2, 1, A_NLoS/A_NLoS_add2*(SINR_threshold*r^alpha_NLoS)^(-1), d1));
                    else
                        intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(2, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), d0));
                    end
                    
                    Laplace_term_NLoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                        2*pi*activated_BS_density_array/d0 * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array/d0 * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                    %                             2*pi*activated_BS_density_array/d0 * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array/d0 * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                else
                    error('unknown path loss model!');
                end
                %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), max_r) - ...
                %                         get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), d0)));
                
            elseif r >= x1 && r < d0
                if strcmpi(LoS_Pr_model, '36.828') == 1 ...
                        || strcmpi(LoS_Pr_model, 'Stepfn') == 1 ...
                        || strcmpi(LoS_Pr_model, '2psRvS') == 1
                    %                     u_array_case1_lw = [r:du:d0-epslong];
                    %                     u_array_case1_up = [d0+epslong:du:max_r];
                    u_array_case1_lw = u_array([r_idx:d0_idx]);
                    u_array_case1_up = u_array([d0_idx+1:end]);
                    
                    r1 = r^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
                    r2 = r^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
                    r1_idx = find(u_array >= r1, 1);
                    r2_idx = find(u_array >= r2, 1);
                    
                    if SINR_threshold_idx == 1
                        geom_density_term_LoS(:, r_idx) = ...
                            -2*pi*BS_density_array.* sum(LoS_Pr_r_array([1:r_idx]) .* ...
                            u_array([1:r_idx]) .* du)-2*pi*BS_density_array .* ...
                            sum(NLoS_Pr_r_array([1:r1_idx]) .* u_array([1:r1_idx]) .* du);
                        
                        geom_density_term_NLoS(:, r_idx) = ...
                            -2*pi*BS_density_array.* sum(LoS_Pr_r_array([1:r2_idx]) .* ...
                            u_array([1:r2_idx]) .* du) - 2*pi*BS_density_array .* ...
                            sum(NLoS_Pr_r_array([1:r_idx]) .* u_array([1:r_idx]) .* du);
                    end
                    
                    intermediate_results{1}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + (SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_lw.^alpha_LoS) * du);
                    intermediate_results{2}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_lw.^alpha_NLoS) * du);
                    intermediate_results{3}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + (SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_up.^alpha_LoS) * du);
                    intermediate_results{4}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1) * u_array_case1_up.^alpha_NLoS) * du);
                    
                    Laplace_term_LoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    
                    
                    intermediate_results{5}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_lw.^alpha_LoS) * du);
                    intermediate_results{6}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([r_idx:d0_idx]) .* ...
                        u_array_case1_lw ./ (1 + (SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_lw.^alpha_NLoS) * du);
                    intermediate_results{7}(SINR_threshold_idx, r_idx) = ...
                        sum(LoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_up.^alpha_LoS) * du);
                    intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                        sum(NLoS_Pr_r_array([d0_idx+1:end]) .* ...
                        u_array_case1_up ./ (1 + (SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case1_up.^alpha_NLoS) * du);
                    
                    Laplace_term_NLoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                elseif strcmpi(LoS_Pr_model, 'Linear') == 1 || strcmpi(LoS_Pr_model, '3psLin') == 1
                    
                    r1 = r^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
                    r2 = r^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
                    
                    if SINR_threshold_idx == 1
                        geom_density_term_LoS(:, r_idx) = ...
                            -1*pi*BS_density_array.*r.^2 + 2*pi*BS_density_array./(3*d0)*(r.^3 - r1.^3);
                        geom_density_term_NLoS(:, r_idx) = ...
                            -1*pi*BS_density_array.*d0.^2 + 2*pi*BS_density_array./(3*d0)*(d0.^3-r.^3);
                    end
                    
                    intermediate_results{1}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_LoS, 1, (SINR_threshold*r^alpha_LoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_LoS, 1, (SINR_threshold*r^alpha_LoS)^(-1), r));
                    intermediate_results{2}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_LoS, 2, (SINR_threshold*r^alpha_LoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_LoS, 2, (SINR_threshold*r^alpha_LoS)^(-1), r));
                    intermediate_results{3}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_NLoS, 2, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_NLoS, 2, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), r1));
                    if strcmpi(LoS_Pr_model, '3psLin') == 1
                        intermediate_results{4}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(1, alpha_NLoS_add1, 1, A_LoS/A_NLoS_add1*(SINR_threshold*r^alpha_LoS)^(-1), d1) - ...
                            get_rou_value(1, alpha_NLoS_add1, 1, A_LoS/A_NLoS_add1*(SINR_threshold*r^alpha_LoS)^(-1), d0)) + ...
                            (get_rou_value(2, alpha_NLoS_add2, 1, A_LoS/A_NLoS_add2*(SINR_threshold*r^alpha_LoS)^(-1), d1));
                    else
                        intermediate_results{4}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(2, alpha_NLoS, 1, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), d0));
                    end
                    
                    
                    
                    Laplace_term_LoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                        2*pi*activated_BS_density_array/d0 * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array/d0 * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                        -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{1}(SINR_threshold_idx, r_idx) + ...
                    %                             2*pi*activated_BS_density_array/d0 * intermediate_results{2}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array/d0 * intermediate_results{3}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{4}(SINR_threshold_idx, r_idx));
                    
                    %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(1, alpha_NLoS, 1, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), max_r) - ...
                    %                         get_rou_value(1, alpha_NLoS, 1, A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1), d0)));
                    
                    %                     intermediate_results{5}(SINR_threshold_idx, r_idx) = ...
                    %                         (get_rou_value(1, alpha_LoS, 1, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), d0) - ...
                    %                         get_rou_value(1, alpha_LoS, 1, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), r));
                    %                     intermediate_results{6}(SINR_threshold_idx, r_idx) = ...
                    %                         (get_rou_value(1, alpha_LoS, 2, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), d0) - ...
                    %                         get_rou_value(1, alpha_LoS, 2, A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1), r));
                    intermediate_results{7}(SINR_threshold_idx, r_idx) = ...
                        (get_rou_value(1, alpha_NLoS, 2, (SINR_threshold*r^alpha_NLoS)^(-1), d0) - ...
                        get_rou_value(1, alpha_NLoS, 2, (SINR_threshold*r^alpha_NLoS)^(-1), r));
                    if strcmpi(LoS_Pr_model, '3psLin') == 1
                        intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(1, alpha_NLoS_add1, 1, A_NLoS/A_NLoS_add1*(SINR_threshold*r^alpha_NLoS)^(-1), d1) - ...
                            get_rou_value(1, alpha_NLoS_add1, 1, A_NLoS/A_NLoS_add1*(SINR_threshold*r^alpha_NLoS)^(-1), d0)) + ...
                            (get_rou_value(2, alpha_NLoS_add2, 1, A_NLoS/A_NLoS_add2*(SINR_threshold*r^alpha_NLoS)^(-1), d1));
                    else
                        intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(2, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), d0));
                    end
                    
                    Laplace_term_NLoS(:, r_idx) = ...
                        (-2*pi*activated_BS_density_array/d0 * intermediate_results{7}(SINR_threshold_idx, r_idx) - ...
                        2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                    % %                         (-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                    % %                         2*pi*activated_BS_density_array/d0 * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                    
                    %                             exp(-2*pi*activated_BS_density_array * intermediate_results{5}(SINR_threshold_idx, r_idx) + ...
                    %                             2*pi*activated_BS_density_array/d0 * intermediate_results{6}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array/d0 * intermediate_results{7}(SINR_threshold_idx, r_idx) + ...
                    %                             -2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    
                else
                    error('unknown path loss model!');
                end
                %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), max_r) - ...
                %                         get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), d0)));
                
            elseif r >= d0
                
                if strcmpi(LoS_Pr_model, '3psLin') == 1
                    
                    if SINR_threshold_idx == 1
                        geom_density_term_LoS(:, r_idx) = deal(-1);
                        geom_density_term_NLoS(:, r_idx) = -1*pi*BS_density_array.*r.^2;
                    end
                    
                    Laplace_term_LoS(:, r_idx) = deal(-1);
                    
                    if r < d1
                        intermediate_results{8}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(1, alpha_NLoS_add1, 1, A_NLoS_add1/A_NLoS_add1*(SINR_threshold*r^alpha_NLoS_add1)^(-1), d1) - ...
                            get_rou_value(1, alpha_NLoS_add1, 1, A_NLoS_add1/A_NLoS_add1*(SINR_threshold*r^alpha_NLoS_add1)^(-1), r)) + ...
                            (get_rou_value(2, alpha_NLoS_add2, 1, A_NLoS_add1/A_NLoS_add2*(SINR_threshold*r^alpha_NLoS_add1)^(-1), d1));
                        
                        Laplace_term_NLoS(:, r_idx) = ...
                            (-2*pi*activated_BS_density_array * intermediate_results{8}(SINR_threshold_idx, r_idx));
                    elseif r >= d1
                        
                        %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(2, alpha_LoS, 1, (SINR_threshold*r^alpha_LoS)^(-1), r))).* ...
                        %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(2, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), r)));
                        
                        intermediate_results{12}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(2, alpha_NLoS_add2, 1, (SINR_threshold*r^alpha_NLoS_add2)^(-1), r));
                        
                        Laplace_term_NLoS(:, r_idx) = ...
                            (-2*pi*activated_BS_density_array * intermediate_results{12}(SINR_threshold_idx, r_idx));
                        %                             exp(-2*pi*activated_BS_density_array * intermediate_results{12}(SINR_threshold_idx, r_idx));
                        
                        %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), max_r) - ...
                        %                         get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), r)));
                    else
                        error('undefined value of "r"');
                    end
                else
                    if strcmpi(LoS_Pr_model, '36.828') == 1 ...
                            || strcmpi(LoS_Pr_model, 'Stepfn') == 1 ...
                            || strcmpi(LoS_Pr_model, '2psRvS') == 1
                        %                     u_array_case2 = [r:du:max_r];
                        u_array_case2 = u_array([r_idx:end]);
                        
                        r1 = r^(alpha_LoS/alpha_NLoS)*((A_NLoS/A_LoS)^(1/alpha_NLoS));
                        r2 = r^(alpha_NLoS/alpha_LoS)*((A_LoS/A_NLoS)^(1/alpha_LoS));
                        r1_idx = find(u_array >= r1, 1);
                        r2_idx = find(u_array >= r2, 1);
                        
                        if SINR_threshold_idx == 1
                            geom_density_term_LoS(:, r_idx) = ...
                                -2*pi*BS_density_array.* sum(LoS_Pr_r_array([1:r_idx]) .* ...
                                u_array([1:r_idx]) .* du)-2*pi*BS_density_array .* ...
                                sum(NLoS_Pr_r_array([1:r1_idx]) .* u_array([1:r1_idx]) .* du);
                            
                            geom_density_term_NLoS(:, r_idx) = ...
                                -2*pi*BS_density_array.* sum(LoS_Pr_r_array([1:r2_idx]) .* ...
                                u_array([1:r2_idx]) .* du) - 2*pi*BS_density_array .* ...
                                sum(NLoS_Pr_r_array([1:r_idx]) .* u_array([1:r_idx]) .* du);
                        end
                        
                        intermediate_results{9}(SINR_threshold_idx, r_idx) = ...
                            sum(LoS_Pr_r_array([r_idx:end]) .* ...
                            u_array_case2 ./ (1 + (SINR_threshold*r^alpha_LoS)^(-1) * u_array_case2.^alpha_LoS) * du);
                        intermediate_results{10}(SINR_threshold_idx, r_idx) = ...
                            sum(NLoS_Pr_r_array([r_idx:end]) .* ...
                            u_array_case2 ./ (1 + A_LoS/A_NLoS*(SINR_threshold*r^alpha_LoS)^(-1) * u_array_case2.^alpha_NLoS) * du);
                        
                        Laplace_term_LoS(:, r_idx) = ...
                            (-2*pi*activated_BS_density_array * intermediate_results{9}(SINR_threshold_idx, r_idx) + ...
                            -2*pi*activated_BS_density_array * intermediate_results{10}(SINR_threshold_idx, r_idx));
                        %                             exp(-2*pi*activated_BS_density_array * intermediate_results{9}(SINR_threshold_idx, r_idx) + ...
                        %                             -2*pi*activated_BS_density_array * intermediate_results{10}(SINR_threshold_idx, r_idx));
                        
                        intermediate_results{11}(SINR_threshold_idx, r_idx) = ...
                            sum(LoS_Pr_r_array([r_idx:end]) .* ...
                            u_array_case2 ./ (1 + A_NLoS/A_LoS*(SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case2.^alpha_LoS) * du);
                        intermediate_results{12}(SINR_threshold_idx, r_idx) = ...
                            sum(NLoS_Pr_r_array([r_idx:end]) .* ...
                            u_array_case2 ./ (1 + (SINR_threshold*r^alpha_NLoS)^(-1) * u_array_case2.^alpha_NLoS) * du);
                        
                        Laplace_term_NLoS(:, r_idx) = ...
                            (-2*pi*activated_BS_density_array * intermediate_results{11}(SINR_threshold_idx, r_idx) + ...
                            -2*pi*activated_BS_density_array * intermediate_results{12}(SINR_threshold_idx, r_idx));
                        
                        %                             exp(-2*pi*activated_BS_density_array * intermediate_results{11}(SINR_threshold_idx, r_idx) + ...
                        %                             -2*pi*activated_BS_density_array * intermediate_results{12}(SINR_threshold_idx, r_idx));
                        
                    elseif strcmpi(LoS_Pr_model, 'Linear') == 1
                        
                        if SINR_threshold_idx == 1
                            geom_density_term_LoS(:, r_idx) = deal(-1);
                            geom_density_term_NLoS(:, r_idx) = -1*pi*BS_density_array.*r.^2;
                        end
                        
                        Laplace_term_LoS(:, r_idx) = deal(-1);
                        %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(2, alpha_LoS, 1, (SINR_threshold*r^alpha_LoS)^(-1), r))).* ...
                        %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(2, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), r)));
                        
                        intermediate_results{12}(SINR_threshold_idx, r_idx) = ...
                            (get_rou_value(2, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), r));
                        
                        Laplace_term_NLoS(:, r_idx) = ...
                            (-2*pi*activated_BS_density_array * intermediate_results{12}(SINR_threshold_idx, r_idx));
                        %                             exp(-2*pi*activated_BS_density_array * intermediate_results{12}(SINR_threshold_idx, r_idx));
                        
                        %                         exp(-2*pi*activated_BS_density_array * (get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), max_r) - ...
                        %                         get_rou_value(1, alpha_NLoS, 1, (SINR_threshold*r^alpha_NLoS)^(-1), r)));
                        
                    else
                        error('unknown path loss model!');
                    end
                end
            else
                error('undefined value of "r"');
            end
        end % end of r loop to get the second-layer integral result
        
        
        for BS_density_idx = 1:len_BS_density_array
            BS_density = BS_density_array(BS_density_idx);
            
            % LoS part
            if ~isempty(r_lower_indices)
                T_cov(BS_density_idx,1) = sum(...
                    exp(-1*SINR_threshold*r_array(1, r_lower_indices).^alpha_LoS*N0/(P_tx(BS_density_idx)*A_LoS) + ...
                    geom_density_term_LoS(BS_density_idx, r_lower_indices) + ...
                    Laplace_term_LoS(BS_density_idx, r_lower_indices)) .* ...
                    2*pi*BS_density .* r_array(1, r_lower_indices) .* ...
                    LoS_Pr_r_array(r_lower_indices) * dr);
            else
                T_cov(BS_density_idx,1) = 0;
            end
            
            if ~isempty(r_middle_indices)
                T_cov(BS_density_idx,2) = sum(...
                    exp(-1*SINR_threshold*r_array(1, r_middle_indices).^alpha_LoS*N0/(P_tx(BS_density_idx)*A_LoS) + ...
                    geom_density_term_LoS(BS_density_idx, r_middle_indices) + ...
                    Laplace_term_LoS(BS_density_idx, r_middle_indices)) .* ...
                    2*pi*BS_density .* r_array(1, r_middle_indices) .* ...
                    LoS_Pr_r_array(r_middle_indices) * dr);
            else
                T_cov(BS_density_idx,2) = 0;
            end
            
            if ~isempty(r_upper_indices)
                T_cov(BS_density_idx,3) = sum(...
                    exp(-1*SINR_threshold*r_array(1, r_upper_indices).^alpha_LoS*N0/(P_tx(BS_density_idx)*A_LoS) + ...
                    geom_density_term_LoS(BS_density_idx, r_upper_indices) + ...
                    Laplace_term_LoS(BS_density_idx, r_upper_indices)) .* ...
                    2*pi*BS_density .* r_array(1, r_upper_indices) .* ...
                    LoS_Pr_r_array(r_upper_indices) * dr);
            else
                T_cov(BS_density_idx,3) = 0;
            end
            
            % NLoS part
            if ~isempty(r_lower_indices)
                T_cov(BS_density_idx,4) = sum(...
                    exp(-1*SINR_threshold*r_array(1, r_lower_indices).^alpha_NLoS*N0/(P_tx(BS_density_idx)*A_NLoS) + ...
                    geom_density_term_NLoS(BS_density_idx, r_lower_indices) + ...
                    Laplace_term_NLoS(BS_density_idx, r_lower_indices)) .* ...
                    2*pi*BS_density .* r_array(1, r_lower_indices) .* ...
                    NLoS_Pr_r_array(r_lower_indices) * dr);
            else
                T_cov(BS_density_idx,4) = 0;
            end
            
            if ~isempty(r_middle_indices)
                T_cov(BS_density_idx,5) = sum(...
                    exp(-1*SINR_threshold*r_array(1, r_middle_indices).^alpha_NLoS*N0/(P_tx(BS_density_idx)*A_NLoS) + ...
                    geom_density_term_NLoS(BS_density_idx, r_middle_indices) + ...
                    Laplace_term_NLoS(BS_density_idx, r_middle_indices)) .* ...
                    2*pi*BS_density .* r_array(1, r_middle_indices) .* ...
                    NLoS_Pr_r_array(r_middle_indices) * dr);
            else
                T_cov(BS_density_idx,5) = 0;
            end
            
            if ~isempty(r_upper_indices)
                T_cov(BS_density_idx,6) = sum(...
                    exp(-1*SINR_threshold*r_array(1, r_upper_indices).^alpha_NLoS*N0/(P_tx(BS_density_idx)*A_NLoS) + ...
                    geom_density_term_NLoS(BS_density_idx, r_upper_indices) + ...
                    Laplace_term_NLoS(BS_density_idx, r_upper_indices)) .* ...
                    2*pi*BS_density .* r_array(1, r_upper_indices) .* ...
                    NLoS_Pr_r_array(r_upper_indices) * dr);
            else
                T_cov(BS_density_idx,6) = 0;
            end
            
            p_cov(BS_density_idx, SINR_threshold_idx) = sum(T_cov(BS_density_idx,:));
        end % end of BS density loop
    end % end of SINR threshold loop (for CCDF)
    
    save(MATname_data);
    
    intermediate_results = [];
    Laplace_term_LoS = [];
    Laplace_term_NLoS = [];
    geom_density_term_LoS = [];
    geom_density_term_NLoS = [];
    
    save(MATname_brief_data);
    
    % plot the figures
    
    QoS_SINR_th = 1;
    % coverage probability vs lambda
    if plot_enabler
        figure(1)
        box on;
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            semilogx(BS_density_array*1e6, p_cov(:, find(SINR_threshold_array >= QoS_SINR_th,1)), 'r');
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            semilogx(BS_density_array*1e6, p_cov(:, find(SINR_threshold_array >= QoS_SINR_th,1)), 'b');
            % semilogx(BS_density_array, p_cov(:, find(SINR_threshold_array_dB==0,1)), 'b');
        else
            error('undefined scheme');
        end
        hold on;
        xlabel('BS density [1/km^2] \it{\lambda}');
        ylabel('Probability of \it{SINR}\rm{>}\it{\gamma}');
        axis([min(BS_density_array*1e6),max(BS_density_array*1e6),0,1]);
        
        % legend('Conventional HPPP', 'Modified HPPP considering LoS/NLoS');
    end
    
    % CCDF of SINR
    if plot_enabler
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
    
    if len_SINR_threshold_array > 1
        x_array = 1/2 * (SINR_threshold_array(1:end-1) + SINR_threshold_array(2:end));
        usable_x_idx = find(x_array >= QoS_SINR_th, 1);
        for BS_density_idx = 1:len_BS_density_array
            %         numerical_SINR_PDF = (1-p_cov(BS_density_idx,2:end)) - (1-p_cov(BS_density_idx,1:end-1));
            %
            %         area_spec_efficiency_exact(BS_density_idx) = ...
            %             sum(log2(1+x_array(usable_x_idx:end)) .* numerical_SINR_PDF(usable_x_idx:end) .* dx_array(usable_x_idx:end)) * BS_density_array(BS_density_idx);
            %         area_spec_efficiency_LB(BS_density_idx) = ...
            %             sum(log2(1+QoS_SINR_th) .* numerical_SINR_PDF(usable_x_idx:end) .* dx_array(usable_x_idx:end)) * BS_density_array(BS_density_idx);
            %         area_spec_efficiency_UB(BS_density_idx) = ...
            %             log2(1 + sum(x_array(usable_x_idx:end) .* numerical_SINR_PDF(usable_x_idx:end) .* dx_array(usable_x_idx:end))) * BS_density_array(BS_density_idx);
            
            CCDF_array = (p_cov(BS_density_idx,1:end-1) + p_cov(BS_density_idx,2:end))/2;
            area_spec_efficiency_exact(BS_density_idx) = BS_density_array(BS_density_idx)*1e6 * (...
                (1/log(2)) * sum(CCDF_array(usable_x_idx:end) ./ (1+x_array(usable_x_idx:end)) .* dx_array(usable_x_idx:end)) + ...
                log2(1+x_array(usable_x_idx))*CCDF_array(usable_x_idx));
        end
        
        
        % coverage probability vs lambda
        if plot_enabler
            figure(3)
            box on;
            if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
                semilogx(BS_density_array*1e6, area_spec_efficiency_exact, 'r');
                hold on;
                % semilogx(BS_density_array*1e6, area_spec_efficiency_LB, 'k');
            elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
                % plot(BS_density_array, area_spec_efficiency_exact, 'b');
                semilogx(BS_density_array*1e6, area_spec_efficiency_exact, 'b');
                hold on;
                % semilogx(BS_density_array, area_spec_efficiency_LB, 'k');
            else
                error('undefined scheme');
            end
            
            xlabel('BS density [1/km^2] \it{\lambda}');
            ylabel('Area spectral efficiency [bps/Hz/km^2]');
            
            legend('Exact value');
        end
    end
    
    
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