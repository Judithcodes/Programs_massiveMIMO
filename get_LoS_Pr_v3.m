% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: get_LoS_Pr_v3.m
% Authors: Ming Ding
% Version: 3.0
% Date: 2015-06-25
% Description: Get the probability of LoS transmission (3 piece included)
% Copyright(c): For personal study only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LoS_Pr] = get_LoS_Pr_v3(SCH, LoS_Pr_model, r_array)

global R1;                  % Constant R_1 in the LoS probability function in [m]
global R2;                  % Constant R_2 in the LoS probability function in [m]
global d0;                  % Constant d_0 in the LoS probability function in [m]
global a;                   % Constant a in the LoS probability function in [m]


LoS_Pr = NaN(1, length(r_array));

for r_idx = 1:length(r_array)
    r = r_array(r_idx);
    
    if strcmpi(LoS_Pr_model, '36.828') == 1
        R1 = 156;
        R2 = R1/(log(10))^2;
        d0 = R1/log(10);
        a = 5;
        
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            if r <= d0
                LoS_Pr(r_idx) = 1 - a * exp(-1*R1./r);
            elseif r > d0
                LoS_Pr(r_idx) = a * exp(-1*r./R2);
            else
                error('undefined value for "LoS_Pr"');
            end
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            LoS_Pr(r_idx) = zeros(1, length(r));
        else
            error('undefined scheme');
        end
        
    elseif strcmpi(LoS_Pr_model, '2psRvS') == 1
        d0 = 20;
        d1 = 120;
        
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            if r <= d0
                LoS_Pr(r_idx) = 1;
            elseif r > d0 && r <= d1
                LoS_Pr(r_idx) = 1 - (r - d0) * (1/(d1 - d0));
            elseif r > d1
                LoS_Pr(r_idx) = 0;
            else
                error('undefined value for "LoS_Pr"');
            end
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            LoS_Pr(r_idx) = zeros(1, length(r));
        else
            error('undefined scheme');
        end
        
    elseif strcmpi(LoS_Pr_model, '36.814') == 1
        d0 = 18;
        R1 = 18;
        R2 = 36;
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            if r <= d0
                LoS_Pr(r_idx) = 1;
            elseif r > d0
                LoS_Pr(r_idx) = 1 - (1-R1/r) .* (1-exp(-1*r/R2));
            else
                error('undefined value for "LoS_Pr"');
            end
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            LoS_Pr(r_idx) = zeros(1, length(r));
        else
            error('undefined scheme');
        end
        
    elseif strcmpi(LoS_Pr_model, 'Linear') == 1
        
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            d0 = 300;
            if r <= d0
                LoS_Pr(r_idx) = 1 - r/d0;
            elseif r > d0
                LoS_Pr(r_idx) = 0;
            else
                error('undefined value for "LoS_Pr"');
            end
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            d0 = 1e-6;
            LoS_Pr(r_idx) = zeros(1, length(r));
        else
            error('undefined scheme');
        end
        
    elseif strcmpi(LoS_Pr_model, 'Stepfn') == 1
        
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            d0 = 300/2;
            if r <= d0
                LoS_Pr(r_idx) = 1;
            elseif r > d0
                LoS_Pr(r_idx) = 0;
            else
                error('undefined value for "LoS_Pr"');
            end
        elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
            d0 = 1e-6;
            LoS_Pr(r_idx) = zeros(1, length(r));
        else
            error('undefined scheme');
        end
        
    elseif strcmpi(LoS_Pr_model, '3psLin') == 1    
        if strcmpi(SCH, 'Prop') == 1                % The proposed scheme
            d0 = 300;
            if r <= d0
                LoS_Pr(r_idx) = 1 - r/d0;
            elseif r > d0
                LoS_Pr(r_idx) = 0;
            else
                error('undefined value for "LoS_Pr"');
            end
%         elseif strcmpi(SCH, 'Jeff') == 1            % The conventional HPPP model
%             d0 = 1e-6;
%             LoS_Pr(r_idx) = zeros(1, length(r));
        else
            error('undefined scheme');
        end
        
    else
        error('undefined LoS probability model');
    end
end