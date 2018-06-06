
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: get_rou_value.m
% Authors: Ming Ding
% Version: 1.0
% Date: 2015-01-05
% Description: Get the integral result for the function \rou
% Copyright(c): For personal study only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rou_value] = get_rou_value(rou_switch, alpha, beta, t, d)

% alpha = 3.75;
% beta = 2;
% t = 13.2;
% d = 0.0001;
% rou_switch = 2;

% rou_switch
% alpha
% beta
% t=0.000016;
% d

computation_method_code = 3;

switch rou_switch
    case 1,
        if computation_method_code == 1 ||computation_method_code == 10
            dx = 0.01;
            x_array = [0:dx:d];
            summation_Q = sum(x_array.^beta ./ (1+t*x_array.^alpha) * dx);
            rou_value = summation_Q;
        end
        
        if computation_method_code == 2 ||computation_method_code == 10
            % Integrate f(x) = x.^beta ./ (1+t*x.^alpha) from 0 to d:
            f = @(x) x.^beta ./ (1+t*x.^alpha);
            numerical_Q = integral(f,0,d);
            rou_value = numerical_Q;
        end
        
        if computation_method_code == 3 ||computation_method_code == 10
            analytical_Q = d^(beta+1) / (beta+1) * ...
                hypergeom([1, (beta+1)/alpha], 1+(beta+1)/alpha, -1*(t*d^alpha));
            rou_value = analytical_Q;
        end
    case 2,
        if computation_method_code == 1 ||computation_method_code == 10
            dx = 0.01;
            x_array = [d:dx:1e5];
            summation_Q = sum(x_array.^beta ./ (1+t*x_array.^alpha) * dx);
            rou_value = summation_Q;
        end
        
        if computation_method_code == 2 ||computation_method_code == 10
            % Integrate f(x) = x.^beta ./ (1+t*x.^alpha) from d to infinity:
            f = @(x) x.^beta ./ (1+t*x.^alpha);
            numerical_Q = integral(f,d,Inf);
            rou_value = numerical_Q;
        end
        
        if computation_method_code == 3 ||computation_method_code == 10
            analytical_Q = d^(beta+1-alpha) / (t*(alpha-beta-1)) * ...
                hypergeom([1, 1-(beta+1)/alpha], 2-(beta+1)/alpha, -1/(t*d^alpha));
            rou_value = analytical_Q;
        end
    otherwise,
        error('undefined rou_switch');
%         rou_value = analytical_Q
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the validity of the hypergeometric function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u = 10;
% v = 2;
% mu = 1.95;
% beta = 3;
% 
% % Integrate f(x) = x.^(mu-1) ./ (1+beta*x).^v from u to infinity:
% f = @(x) x.^(mu-1) ./ (1+beta*x).^v
% 
% dx = 0.01;
% x_array = [u:dx:1e5];
% 
% summation_Q = sum(x_array.^(mu-1) ./ (1+beta*x_array).^v * dx)
% 
% numerical_Q = integral(f,u,Inf)
% 
% analytical_Q = u^(mu-v) / (beta^v * (v-mu)) * hypergeom([v, v-mu], v-mu+1, -1/(beta*mu))
% 





% u = 10;
% v = 1;
% mu = 2.05;
% beta = 3;
% 
% % Integrate f(x) = x.^(mu-1) ./ (1+beta*x).^v from 0 to u:
% f = @(x) x.^(mu-1) ./ (1+beta*x).^v
% 
% dx = 0.01;
% x_array = [0:dx:u];
% 
% summation_Q = sum(x_array.^(mu-1) ./ (1+beta*x_array).^v * dx)
% 
% numerical_Q = integral(f,0,u)
% 
% analytical_Q = u^(mu) / mu * hypergeom([v, mu], 1+mu, -1*beta*u)
% 


