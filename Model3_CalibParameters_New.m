cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
clc

%%% Calibrate to two points (1993 and 2022)
%% Declare parameters values (determined outside of the model)
% Z_l_0, g_z_l, P_0, g_P, Mu, SD

Z_l_0 = 0.2030;
g_z_l = 0.8285; %(growth over the 30 years period)

Mu = 0.5;

SD = 1;



%% Exogenous variables 
% Land series
T_series = [6.8763, 5.9502];

% Price series
P_series = [0.9082,  1.4156];

% Z_l series
Z_l_series = [Z_l_0, Z_l_0*(1+ g_z_l)];


%% Set-up the grid of parameters for searching
% Z_m_0, H_bar_0, H_bar_0
% g_z_m, g_z_t, g_H_bar

Z_m_0_list = logspace(log10(0.1), log10(4), 75);
g_z_m_list = linspace(0.2,1.2,51);

Z_t_0_list = logspace(log10(0.01), log10(3), 60);
g_z_t_list = linspace(0.2,1.2,51);

H_bar_0_list = logspace(log10(0.1), log10(100), 200);
g_H_bar_list = linspace(0.2,1.2,51);

%% Declare target moments

L_a0 = 0.5288;
L_a1 = 0.3021;

W_a_0_to_W_m_rw_0 = 0.5;
W_a_1_to_W_m_rw_1 = 0.4178;

Y_a_0_to_Y_0 = 0.0878;
Y_a_1_to_Y_1 = 0.0633;

Actual_MM = [L_a0, W_a_0_to_W_m_rw_0, Y_a_0_to_Y_0; L_a1, W_a_1_to_W_m_rw_1, Y_a_1_to_Y_1];

%% Weight Matrix. Give less weight to relative wage since it is less reliable and not the main focus.
Weight = [0.175, 0.15, 0.175; 0.175, 0.15, 0.175];

%% Do random serch.
% Draw random combination of parameters. Do it 500,000 times. 
% There are 119,385,900,000 possible combination of parameter values. 
n = 10000; %(10,000 for now)

% set seed
rng(2024);

% Matrix to store results.
Search_results_Mat = zeros(n,7);

i = 1;
while i <= n
    % Randomly pick parameter values
    Z_m_0 = randsample(Z_m_0_list, 1);
    g_z_m = randsample(g_z_m_list, 1);
    Z_t_0 = randsample(Z_t_0_list, 1);
    g_z_t = randsample(g_z_t_list, 1);
    H_bar_0 = randsample(H_bar_0_list, 1);
    g_H_bar = randsample(g_H_bar_list, 1);
    
    % Record Parameters
    Param_combi = [Z_m_0; g_z_m; Z_t_0; g_z_t; H_bar_0; g_H_bar];
    Search_results_Mat(i,1:6) = Param_combi;

    % Create exogenous variables (Z_mt, Z_tt, H_bar)
    Z_m_series = [Z_m_0, Z_m_0 * (1+g_z_m)];
    Z_t_series = [Z_t_0, Z_t_0 * (1+g_z_t)];
    H_bar_series = [H_bar_0, H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (2*3)
    soln_mat = zeros(2,3);

    % Solve the model for the two period.
    for t = 1:2
        Z_lt = Z_l_series(t);
        Z_tt = Z_t_series(t);
        Z_mt = Z_m_series(t);

        Hbar_t = H_bar_series(t);
        P_t = P_series(t);
        T_t = T_series(t);
        
        % Choose a reasonable starting point
        x0 = [0.5; Z_lt; Z_mt];

        % Find the equilibrium
        x_star = fsolve(@(x)Final_Model_Function(x, Z_lt, Z_tt, Z_mt, Mu, Hbar_t, SD, P_t, T_t), x0);

        % Calculate W_m * avg h in sector m
        h_in_m = (Hbar_t +  ...
        SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)));

        Real_world_W_m = x_star(3)*h_in_m;

        % Calculate Y_a_t, Y_m_t
        Y_a_t = ((Z_lt*x_star(1))^((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu))...
        ^(Mu/(Mu - 1));
        Y_m_t = Z_mt * h_in_m;

        % Calculate the simulated moments
        L_a_sim = x_star(1);
        W_a_to_W_m_rw_sim = x_star(2)/Real_world_W_m;
        Y_a_to_Y_sim = Y_a_t/(Y_a_t + Y_m_t);

        % Store in the matrix.
        soln_mat(t,:) = [L_a_sim, W_a_to_W_m_rw_sim, Y_a_to_Y_sim];
    end
        % Calculate Loss (percentage difference squared)
        Loss_mat = ((soln_mat - Actual_MM)./Actual_MM).^2;

        % Summarize by using weighted average
        avg_loss = mean(Weight.*Loss_mat,"all");

        % Store the loss.
        Search_results_Mat(i,7) = avg_loss;

    i = i+1;
end  

%% See what's the best combination of parameters.
% get location
[minloss, index] = min(Search_results_Mat(:,7));

% get the parameter values.
Best_para = Search_results_Mat(index,:);

% Make it into a nice table
Best_para_tab = array2table(Best_para, "VariableNames", {'Z_m_0', 'g_z_m', 'Z_t_0', 'g_z_t', 'H_bar_0', 'g_H_bar', 'Loss'});

%% Solve the model with the calibrated parameters.
Best_Z_m_0 = Best_para(1);
Best_g_z_m = Best_para(2);
Best_Z_t_0 = Best_para(3);
Best_g_z_t = Best_para(4);
Best_H_bar_0 = Best_para(5);
Best_g_H_bar = Best_para(6);

% Create exogenous variables (Z_mt, Z_tt, H_bar)
Z_m_series = [Best_Z_m_0, Best_Z_m_0 * (1+Best_g_z_m)];
Z_t_series = [Best_Z_t_0, Best_Z_t_0 * (1+Best_g_z_t)];
H_bar_series = [Best_H_bar_0, Best_H_bar_0*(1+Best_g_H_bar)];

% Matrix to store results in each period (2*3)
soln_mat_calibrated = zeros(2,3);

% Solve the model for the two period.
for t = 1:2
    Z_lt = Z_l_series(t);
    Z_tt = Z_t_series(t);
    Z_mt = Z_m_series(t);

    Hbar_t = H_bar_series(t);
    P_t = P_series(t);
    T_t = T_series(t);
    
    % Choose a reasonable starting point
    x0 = [0.5; Z_lt; Z_mt];

    % Find the equilibrium
    x_star = fsolve(@(x)Final_Model_Function(x, Z_lt, Z_tt, Z_mt, Mu, Hbar_t, SD, P_t, T_t), x0);

    % Calculate W_m * avg h in sector m
    h_in_m = (Hbar_t +  ...
    SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)));

    Real_world_W_m = x_star(3)*h_in_m;

    % Calculate Y_a_t, Y_m_t
    Y_a_t = ((Z_lt*x_star(1))^((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu))...
    ^(Mu/(Mu - 1));
    Y_m_t = Z_mt * h_in_m;

    % Calculate the simulated moments
    L_a_sim = x_star(1);
    W_a_to_W_m_rw_sim = x_star(2)/Real_world_W_m;
    Y_a_to_Y_sim = Y_a_t/(Y_a_t + Y_m_t);

    % Store in the matrix.
    soln_mat_calibrated(t,:) = [L_a_sim, W_a_to_W_m_rw_sim, Y_a_to_Y_sim];
end
